import os
from datetime import datetime
from flask import Blueprint, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
from bson import ObjectId
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.schema import Document
from extensions import mongo
from utils import token_required
from langchain_ollama import OllamaEmbeddings  # import embeddings
from utils.faiss_index import delete_from_index  # fonction suppression index FAISS
from data_ingestion import (search_arxiv, search_pubmed,
                            search_openalex, merge_and_deduplicate)  # fonctions de recherche
from langchain_community.vectorstores import FAISS
import PyPDF2
import pdfplumber
from PIL import Image
bp_documents = Blueprint("documents", __name__)

UPLOAD_FOLDER = "uploads"
THUMBNAIL_FOLDER = "thumbnails"
ALLOWED_EXTENSIONS = {"pdf", "txt", "docx", "png", "jpg"}


def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS

def create_thumbnail(filepath, save_path, size=(200, 200)):
    ext = filepath.rsplit('.', 1)[1].lower()
    try:
        if ext == 'pdf':
            # Extraire la première page du PDF comme image
            with pdfplumber.open(filepath) as pdf:
                page = pdf.pages[0]
                pil_image = page.to_image(resolution=150).original
                pil_image.thumbnail(size)
                pil_image.save(save_path)
        elif ext in ['jpg', 'jpeg', 'png', 'gif']:
            img = Image.open(filepath)
            img.thumbnail(size)
            img.save(save_path)
        else:
            # Pas de thumbnail pour autres fichiers
            pass
    except Exception as e:
        print(f"Erreur création thumbnail pour {filepath} : {e}")


def extract_pdf_content(filepath):
    text = ""
    try:
        with pdfplumber.open(filepath) as pdf:
            for page in pdf.pages:
                page_text = page.extract_text()
                if page_text:
                    # Nettoyage minimal : retirer espaces en trop
                    cleaned_text = page_text.strip()
                    text += cleaned_text + "\n"
    except Exception as e:
        print(f"Erreur extraction PDF {filepath} : {e}")
    
    if not text.strip():
        # Cas où pdfplumber ne renvoie rien : essayer PyPDF2
        try:
            with open(filepath, "rb") as f:
                reader = PyPDF2.PdfReader(f)
                for page in reader.pages:
                    page_text = page.extract_text()
                    if page_text:
                        text += page_text.strip() + "\n"
        except Exception as e2:
            print(f"PyPDF2 aussi a échoué sur {filepath} : {e2}")

    return text.strip()

@bp_documents.route("/upload", methods=["POST"])
@token_required
def upload_documents(user):
    conversation_id = request.form.get("conversation_id")
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

            content = ""
            ext = filename.rsplit(".", 1)[1].lower()
            if ext == "pdf":
                content = extract_pdf_content(filepath)
            else:
                try:
                    with open(filepath, "r", encoding="utf-8") as ff:
                        content = ff.read()
                except Exception as e:
                    print(f"Erreur lecture fichier {filename}: {e}")

            if content:
                text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=100)
                chunks = text_splitter.split_text(content)
                for i, chunk in enumerate(chunks):
                    new_documents.append(Document(page_content=chunk, metadata={"source": f"{filename}_chunk_{i}"}))

            user_thumb_folder = os.path.join(THUMBNAIL_FOLDER, str(user["_id"]))
            os.makedirs(user_thumb_folder, exist_ok=True)
            thumb_path = os.path.join(user_thumb_folder, f"{filename}.png")
            create_thumbnail(filepath, thumb_path)

            doc_entry = {
                "filename": filename,
                "path": filepath,
                "thumbnail": f"{filename}.png",
                "user_id": user["_id"],
                "conversation_id": conversation_id,
                "content": content,
                "upload_date": datetime.utcnow()
            }
            inserted = mongo.db.documents.insert_one(doc_entry)
            doc_entry["_id"] = str(inserted.inserted_id)
            saved_files.append(doc_entry)

            # Ajouter à la liste pour indexation FAISS (avec contenu complet)
            new_doc = Document(page_content=content, metadata={"source": filepath})
            new_documents.append(new_doc)

    if not new_documents:
        return jsonify({"error": "Aucun fichier valide à indexer"}), 400

    embedding = OllamaEmbeddings(model="nomic-embed-text")
    index_path = os.path.join(user_folder, "faiss_index")

    try:
        store = FAISS.load_local(index_path, embedding=embedding, allow_dangerous_deserialization=True)
        store.add_documents(new_documents)
    except Exception:
        store = FAISS.from_documents(new_documents, embedding=embedding)

    store.save_local(index_path)

    return jsonify({"files": saved_files})


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


@bp_documents.route("/delete_document/<doc_id>", methods=["DELETE"])
@token_required
def delete_document(user, doc_id):
    doc = mongo.db.documents.find_one({"_id": ObjectId(doc_id), "user_id": user["_id"]})
    if not doc:
        return jsonify({"error": "Document introuvable"}), 404

    if os.path.exists(doc["path"]):
        os.remove(doc["path"])

    thumb_path = os.path.join(THUMBNAIL_FOLDER, str(user["_id"]), doc["thumbnail"])
    if os.path.exists(thumb_path):
        os.remove(thumb_path)

    success = delete_from_index(doc["path"], user["_id"])

    mongo.db.documents.delete_one({"_id": ObjectId(doc_id)})

    return jsonify({"message": "Document supprimé", "faiss_deleted": success})


@bp_documents.route('/fetch_documents', methods=['POST'])
@token_required
def fetch_documents(user):
    data = request.get_json()
    query = data.get('query')
    source = data.get('source')
    max_results = data.get('max_results', 5)

    documents = []
    if source == 'all':
        all_docs = []
        all_docs.extend(search_arxiv(query, max_results))
        all_docs.extend(search_pubmed(query, max_results))
        # all_docs.extend(search_semantic_scholar(query, max_results))
        all_docs.extend(search_openalex(query, max_results))
        documents = merge_and_deduplicate(all_docs)
    elif source == 'arxiv':
        documents = search_arxiv(query, max_results)
    elif source == 'pubmed':
        documents = search_pubmed(query, max_results)
    # elif source == 'semantic_scholar':
    #     documents = search_semantic_scholar(query, max_results)
    elif source == 'openalex':
        documents = search_openalex(query, max_results)
    else:
        # gérer source inconnue si besoin
        documents = []

    documents = sorted(documents, key=lambda d: d.get('title', '').lower())
    return jsonify({'documents': documents})
