import os
import re

from datetime import datetime
from flask import Blueprint, request, jsonify, send_from_directory,current_app,Response
from werkzeug.utils import secure_filename
from bson import ObjectId
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.schema import Document
from extensions import mongo
from utils import token_required
from langchain_ollama import OllamaEmbeddings  # import embeddings
from utils.faiss_index import delete_from_index, update_faiss_index,CachedEmbeddingModel
from data_ingestion import (search_arxiv, search_pubmed,
                            search_openalex, merge_and_deduplicate)  # fonctions de recherche
from langchain_community.vectorstores import FAISS
import requests



from utils.paths import get_user_folder, get_user_index_folder
from utils.extract_text_from_pdf import extract_pdf_content

from pdf2image import convert_from_path

from .multi_functions_routes import ollama_query
from deep_translator import GoogleTranslator
from utils.traitement_document import extract_text_from_pdf,extract_file_content,create_thumbnail,find_best_fuzzy_match,translate_text,allowed_file,calculate_file_hash

bp_documents = Blueprint("documents", __name__)
BASE_DIR = os.path.dirname(__file__)
UPLOAD_FOLDER = "documents"
THUMBNAIL_FOLDER = "thumbnails"

EXTERNAL_UPLOAD_FOLDER = "./external_uploads"



@bp_documents.route("/upload", methods=["POST"])
@token_required
def upload_documents(user):
    data = request.form or {}
    conversation_id = data.get("conversation_id")
    if not conversation_id:
        return jsonify({"error": "conversation_id manquant"}), 400

    if "files" not in request.files:
        return jsonify({"error": "Aucun fichier"}), 400

    files = request.files.getlist("files")
    saved_files = []

    user_folder = os.path.join(UPLOAD_FOLDER, str(user["_id"]))
    os.makedirs(user_folder, exist_ok=True)

    new_documents = []

    for f in files:
        if f and allowed_file(f.filename):
            filename = secure_filename(f.filename)
            filepath = os.path.join(user_folder, filename)
            f.save(filepath)
            chunks = extract_file_content(filepath)

            if chunks:
                user_thumb_folder = os.path.join(THUMBNAIL_FOLDER, str(user["_id"]))
                os.makedirs(user_thumb_folder, exist_ok=True)
                thumb_filename = f"{os.path.splitext(filename)[0]}.png"
                thumb_path = os.path.join(user_thumb_folder, thumb_filename)
                create_thumbnail(filepath, thumb_path)
            file_hash = calculate_file_hash(filepath)

            # Vérifier doublon
            if mongo.db.documents.find_one({"user_id": user["_id"], "conversation_id": conversation_id, "file_hash": file_hash}):
                print(f"[INFO] Fichier déjà uploadé : {filename}")
                continue

            file_id = ObjectId()
            for i, chunk in enumerate(chunks):
                doc_entry = {
                    "filename": filename,
                    "path": filepath,
                    "user_id": user["_id"],
                    "conversation_id": conversation_id,
                    "content": chunk,
                    "chunk_index": i,
                    "total_chunks": len(chunks),
                    "file_type": os.path.splitext(filename)[1][1:].upper(),
                    "upload_date": datetime.utcnow(),
                    "thumbnail": thumb_filename,
                    "is_master": i == 0,
                    "file_id": file_id,
                    "file_hash": file_hash,
                }
                inserted = mongo.db.documents.insert_one(doc_entry)
                if i == 0:
                    saved_files.append(doc_entry)

                new_documents.append({
                    "id": str(inserted.inserted_id),
                    "content": chunk,
                    "metadata": {
                        "filename": filename,
                        "chunk_index": i,
                        "conversation_id": conversation_id,
                        "file_type": doc_entry["file_type"],
                        "file_id": str(file_id),
                        "is_master": i == 0
                    }
                })

                if new_documents:
                    index_folder = os.path.join(user_folder, "conversations", str(conversation_id))
                    os.makedirs(index_folder, exist_ok=True)
                    update_faiss_index(index_folder, new_documents, user_id=user["_id"])

                return jsonify({"files": saved_files})


import os

def sanitize_filename(name):
    # Remplace les caractères interdits sous Windows et limite la longueur
    name = re.sub(r'[\\/*?:"<>|]', "_", name)
    return name[:240]


@bp_documents.route("/download_external_document", methods=["POST"])
@token_required
def download_external_document(user):
    def arxiv_to_pdf_url(url):
        if "arxiv.org/abs/" in url:
            return url.replace("arxiv.org/abs/", "arxiv.org/pdf/") + ".pdf"
        return url

    data = request.get_json()
    file_url = data.get("file_url")
    filename = data.get("filename")
    source = data.get("source", "unknown")

    if not file_url or not filename:
        return jsonify({"error": "URL ou nom de fichier manquant"}), 400

    file_url = arxiv_to_pdf_url(file_url)
    filename = sanitize_filename(filename)

    if not filename.lower().endswith(".pdf"):
        filename += ".pdf"

    try:
        response = requests.get(file_url, allow_redirects=True, timeout=10)
        if response.status_code != 200:
            return jsonify({"error": f"Erreur téléchargement : {response.status_code}"}), 400

        content_type = response.headers.get("Content-Type", "")
        if not content_type.startswith("application/pdf"):
            return jsonify({"error": "Le fichier téléchargé n'est pas un PDF valide."}), 400

        user_folder = os.path.join(EXTERNAL_UPLOAD_FOLDER, str(user["_id"]))
        os.makedirs(user_folder, exist_ok=True)
        filepath = os.path.join(user_folder, filename)

        with open(filepath, "wb") as f:
            f.write(response.content)
        print(f"[Téléchargement] Fichier enregistré à : {filepath}")

        content_data = extract_pdf_content(filepath)
        text_content = content_data["text_sections"][:100_000]

        text_splitter = RecursiveCharacterTextSplitter(chunk_size=1500, chunk_overlap=50)
        chunks = text_splitter.split_text(text_content)[:200]
        print(f"[Chunking] Nombre de chunks générés : {len(chunks)}")

        documents = [
            Document(page_content=chunk, metadata={"source": f"{filename}_chunk_{i}"})
            for i, chunk in enumerate(chunks)
        ]

        index_path = os.path.join(user_folder, "faiss_index")
        embedding = CachedEmbeddingModel(user["_id"])

        try:
            store = FAISS.load_local(index_path, embeddings=embedding, allow_dangerous_deserialization=True)
            store.add_documents(documents)
            print("[FAISS] Index existant mis à jour.")
        except Exception:
            store = FAISS.from_documents(documents, embedding=embedding)
            print("[FAISS] Nouvel index créé.")

        store.save_local(index_path)
        print("[FAISS] Index sauvegardé.")

        original_url = data.get("original_url", file_url)

        doc_entry = {
            "_id": original_url,
            "title": filename,
            "path": filepath,
            "source": source,
            "content": text_content,
            "tables_count": content_data["tables_count"],
            "figures_count": content_data["figures_count"],
            "upload_date": datetime.utcnow(),
            "indexed": True,
            "url": file_url,
        }

        mongo = current_app.mongo
        mongo.db.documents_externes.replace_one(
            {"_id": original_url},
            doc_entry,
            upsert=True
        )

        inserted_doc = mongo.db.documents_externes.find_one({"_id": original_url})
        print("[MongoDB] Document inséré ou mis à jour :", inserted_doc)

        return jsonify({"message": "Téléchargement et indexation OK", "document": doc_entry})

    except Exception as e:
        print(f"[Erreur] Exception lors de l’indexation : {e}")
        return jsonify({"error": f"Erreur serveur : {str(e)}"}), 500


@bp_documents.route("/document_chunks", methods=["POST"])
@token_required
def get_document_chunks(user):
    try:
        data = request.get_json()
        document_id = data.get("document_id")
        language = data.get("language", "en")  # Récupère la langue, défaut 'en'

        if not document_id:
            return jsonify({"error": "document_id requis"}), 400

        # Conversion URL arXiv
        if "arxiv.org/abs/" in document_id:
            document_id = document_id.replace("arxiv.org/abs/", "arxiv.org/pdf/") + ".pdf"

        # Recherche en MongoDB
        doc = mongo.db.documents_externes.find_one({"_id": document_id})

        # Téléchargement + extraction si absent
        if not doc:
            filepath = os.path.join(EXTERNAL_UPLOAD_FOLDER, str(user["_id"]), os.path.basename(document_id))
            os.makedirs(os.path.dirname(filepath), exist_ok=True)

            if not os.path.exists(filepath):
                response = requests.get(document_id, stream=True, timeout=15)
                response.raise_for_status()

                content_type = response.headers.get("Content-Type", "")
                if not content_type.startswith("application/pdf"):
                    return jsonify({"error": "Le fichier distant n'est pas un PDF valide."}), 400

                with open(filepath, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

            content_data = extract_pdf_content(filepath)
            text_content = content_data["text_sections"][:100_000]

            doc = {
                "_id": document_id,
                "title": os.path.basename(filepath),
                "path": filepath,
                "source": "unknown",
                "content": text_content,
                "tables_count": content_data["tables_count"],
                "figures_count": content_data["figures_count"],
                "upload_date": datetime.utcnow(),
                "indexed": False,
                "url": document_id
            }
            mongo.db.documents_externes.insert_one(doc)

        # Traduction du contenu si nécessaire
        content_to_use = doc.get("content", "")
        if language != "en":
            try:
                content_to_use = translate_text(content_to_use, language)
            except Exception as e:
                print(f"[⚠️ Traduction contenu échouée] : {e}")

        # Indexation FAISS (sans traduction, on indexe le texte original)
        if not doc.get("indexed", False):
            text_splitter = RecursiveCharacterTextSplitter(chunk_size=1500, chunk_overlap=50)
            chunks = text_splitter.split_text(doc["content"])[:200]

            documents = [
                Document(page_content=chunk, metadata={"source": f"{doc['title']}_chunk_{i}"})
                for i, chunk in enumerate(chunks)
            ]

            user_folder = os.path.join(EXTERNAL_UPLOAD_FOLDER, str(user["_id"]))
            os.makedirs(user_folder, exist_ok=True)
            index_path = os.path.join(user_folder, "faiss_index")

            embedding = CachedEmbeddingModel(user["_id"])

            try:
                store = FAISS.load_local(index_path, embeddings=embedding, allow_dangerous_deserialization=True)
                store.add_documents(documents)
            except Exception:
                store = FAISS.from_documents(documents, embedding=embedding)

            store.save_local(index_path)

            mongo.db.documents_externes.update_one(
                {"_id": document_id},
                {"$set": {"indexed": True}}
            )

        # Si on a déjà highlight_chunks en base, renvoyer
        if "highlight_chunks" in doc:
            full_text = " ".join(doc.get("content", "").splitlines())[:20000]
            valid_chunks = [chunk for chunk in doc["highlight_chunks"] if chunk.get("text") in full_text]

            # Traduction des extraits avant renvoi
            if language != "en":
                for chunk in valid_chunks:
                    try:
                        chunk["text"] = translate_text(chunk["text"], language)
                        if chunk.get("label"):
                            chunk["label"] = translate_text(chunk["label"], language).upper()
                    except Exception as e:
                        print(f"[⚠️ Traduction highlight_chunks échouée] : {e}")

            return jsonify({
                "chunks": valid_chunks,
                "full_text": content_to_use[:20000]  # texte potentiellement traduit
            })

        # Construction prompt avec contenu traduit (ou original si en='en')
        prompt = f"""
        Tu es un expert en lecture de documents scientifiques. Voici le contenu partiel d’un document.
        Identifie les passages les plus pertinents concernant l’objectif, les méthodes, les résultats et la conclusion.

        Pour chaque passage trouvé, ajoute un titre en majuscules (ex: OBJECTIF, METHODE, RESULTATS, CONCLUSION).
        Donne-moi une liste de 3 à 5 extraits pertinents dans ce format :

        OBJECTIF:
        [passage]

        METHODE:
        [passage]

        ...

        Document :
        {content_to_use[:5000]}
        """

        response = ollama_query(prompt)
        raw_chunks = [chunk.strip() for chunk in response.strip().split("\n\n") if chunk.strip()]
        if raw_chunks and raw_chunks[0].lower().startswith("here are"):
            raw_chunks = raw_chunks[1:]

        full_text = " ".join(doc.get("content", "").splitlines())[:20000]

        chunks = [chunk.replace('\n', ' ').strip() for chunk in raw_chunks]

        valid_chunks_with_pos = []
        for chunk in chunks:
            label_match = re.match(r"^(OBJECTIF|METHODE|RESULTATS|CONCLUSION):", chunk, re.IGNORECASE)
            label = label_match.group(1).capitalize() if label_match else None
            clean_chunk = re.sub(r"^(OBJECTIF|METHODE|RESULTATS|CONCLUSION):", "", chunk, flags=re.IGNORECASE).strip()

            matched_text, start, end = find_best_fuzzy_match(full_text, clean_chunk, min_ratio=0.5)
            if matched_text:
                valid_chunks_with_pos.append({
                    "text": matched_text,
                    "start": start,
                    "end": end,
                    "label": label
                })
            else:
                print("No match found.")

        # Traduction des extraits trouvés avant sauvegarde et envoi
        if language != "en":
            for chunk in valid_chunks_with_pos:
                try:
                    chunk["text"] = translate_text(chunk["text"], language)
                    if chunk.get("label"):
                        chunk["label"] = translate_text(chunk["label"], language).upper()
                except Exception as e:
                    print(f"[⚠️ Traduction des extraits échouée] : {e}")

        mongo.db.documents_externes.update_one(
            {"_id": document_id},
            {"$set": {"highlight_chunks": valid_chunks_with_pos}}
        )

        return jsonify({
            "chunks": valid_chunks_with_pos,
            "full_text": content_to_use[:20000]
        })

    except Exception as e:
        print(f"[Erreur] Exception dans /document_chunks : {e}")
        return jsonify({"error": str(e)}), 500



@bp_documents.route('/api/favorites', methods=['POST'])
@token_required
def add_favorite(user):
    data = request.get_json()
    document_id = data.get('documentId')
    document_data = data.get('documentData')  # optionnel, pour infos du doc

    if not document_id:
        return jsonify({'error': 'documentId manquant'}), 400

    # Recherche document dans la base
    doc = mongo.db.documents_externes.find_one({"_id": document_id})

    if not doc:
        # Si document non trouvé, on l'ajoute à la base (avec infos minimales)
        # Si tu n'as pas toutes les données, prends au moins ce que tu peux
        if not document_data:
            return jsonify({'error': 'Document non trouvé et pas assez d\'infos pour l\'ajouter'}), 400

        # Conversion _id string en ObjectId
        document_data["_id"] = document_id
        mongo.db.documents_externes.insert_one(document_data)

    # Ajout aux favoris de l'utilisateur
    user_id = user["_id"]

    mongo.db.users_profile.update_one(
        {"_id": user_id},
        {"$addToSet": {"favorites": document_id}},
        upsert=True
    )

    return jsonify({'message': 'Document ajouté aux favoris'})


@bp_documents.route('/favorites/<document_id>', methods=['DELETE'])
@token_required
def remove_favorite(user, document_id):
    user_id = user["_id"]
    mongo.db.users.update_one(
        {"_id": user_id},
        {"$pull": {"favorites": document_id}}
    )
    return jsonify({'message': 'Document retiré des favoris'})


@bp_documents.route('/favorites', methods=['GET'])
@token_required
def get_favorites(user):
    user_id = user["_id"]

    # Récupérer la liste des _id favoris de l'utilisateur
    user_data = mongo.db.users_profile.find_one({"_id": user_id}, {"favorites": 1})
    favorites_ids = user_data.get("favorites", [])

    # Chercher les documents dans documents_externes
    favorites_docs = list(
        mongo.db.documents_externes.find({"_id": {"$in": favorites_ids}})
    )

    # Pour renvoyer au frontend, on peut transformer ObjectId en string si besoin
    for doc in favorites_docs:
        doc["_id"] = str(doc["_id"])

    return jsonify({"favorites": favorites_docs})


@bp_documents.route('/thumbnails/<user_id>/<filename>')
def serve_thumbnail(user_id, filename):
    thumb_path = os.path.join(THUMBNAIL_FOLDER, user_id)
    return send_from_directory(thumb_path, filename)


@bp_documents.route('/user_documents/<conversation_id>', methods=['GET'])
@token_required
def get_user_documents(user, conversation_id):
    docs = list(mongo.db.documents.find({"conversation_id": conversation_id, "user_id": user["_id"]}))
    for doc in docs:
        doc["_id"] = str(doc["_id"])
        if isinstance(doc["user_id"], ObjectId):
            doc["user_id"] = str(doc["user_id"])
    return jsonify({"files": docs})

@bp_documents.route("/documents", methods=["GET"])
@token_required
def get_user_all_documents(user):
    docs = list(mongo.db.documents.find({"user_id": user["_id"]}))
    
    formatted = []
    for doc in docs:
        upload_date = doc.get("upload_date")
        # Convertir en ISO string si upload_date existe
        upload_date_str = upload_date.isoformat() if upload_date else None
        
        formatted.append({
            "id": str(doc["_id"]),
            "title": doc.get("title", ""),
            "filename": doc.get("filename", ""),
            "uploaded_at": upload_date_str
        })

    return jsonify(formatted)


@bp_documents.route("/delete_document/<doc_id>", methods=["DELETE"])
@token_required
def delete_document(user, doc_id):
    try:
        object_id = ObjectId(doc_id)
    except Exception:
        return jsonify({"error": "ID de document invalide"}), 400

    doc = mongo.db.documents.find_one({"_id": object_id, "user_id": user["_id"]})
    if not doc:
        return jsonify({"error": "Document introuvable"}), 404

    file_id = doc.get("file_id")

    if file_id:
        # Récupérer tous les chunks liés
        chunks = list(mongo.db.documents.find({
            "file_id": file_id,
            "user_id": user["_id"]
        }))
        chunk_ids = [str(chunk["_id"]) for chunk in chunks]

        # Supprimer les documents (chunks) dans la base
        result = mongo.db.documents.delete_many({
            "file_id": file_id,
            "user_id": user["_id"]
        })

        # Supprimer le cache des embeddings associés
        try:
            mongo.db.embeddings_cache.delete_many({
                "chunk_id": {"$in": chunk_ids},
                "user_id": user["_id"]
            })
        except Exception as e:
            print(f"Erreur suppression cache embeddings : {e}")

        # Supprimer fichiers physiques et thumbnail si master
        if doc.get("is_master"):
            try:
                if os.path.exists(doc["path"]):
                    os.remove(doc["path"])
            except Exception as e:
                print(f"Erreur suppression fichier : {e}")

            try:
                thumb_path = os.path.join(
                    THUMBNAIL_FOLDER,
                    str(user["_id"]),
                    doc["thumbnail"]
                )
                if os.path.exists(thumb_path):
                    os.remove(thumb_path)
            except Exception as e:
                print(f"Erreur suppression thumbnail : {e}")

            # Supprimer de l'index FAISS
            try:
                delete_from_index(user_id=user["_id"], doc_ids=chunk_ids)
            except Exception as e:
                print(f"Erreur suppression FAISS : {e}")

        return jsonify({
            "message": "Document, chunks, indexation et cache supprimés",
            "deleted_count": result.deleted_count
        }), 200

    # Fallback pour les documents anciens sans file_id
    mongo.db.documents.delete_one({"_id": object_id})
    mongo.db.embeddings_cache.delete_many({"chunk_id": str(object_id), "user_id": user["_id"]})

    return jsonify({"message": "Document supprimé (ancien format + cache)"}), 200





@bp_documents.route("/proxy_pdf")
def proxy_pdf():
    url = request.args.get("url")
    if not url:
        return "Missing url parameter", 400
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200 or "pdf" not in r.headers.get("Content-Type", ""):
            return "Not a valid PDF file", 400
        return Response(r.content, mimetype="application/pdf")
    except Exception as e:
        return str(e), 500