import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import re
import time
import traceback
import eventlet

eventlet.monkey_patch()

from datetime import datetime, timedelta
import pdfplumber
import numpy as np
from bson import ObjectId
import ollama
import requests
import PyPDF2
import jwt as pyjwt
from flask import Flask, request, jsonify, send_from_directory, Blueprint, url_for
from flask_cors import CORS, cross_origin
from flask_socketio import SocketIO, emit, disconnect
from werkzeug.utils import secure_filename
from flask_jwt_extended import JWTManager
from flask_bcrypt import Bcrypt
from flask_pymongo import PyMongo
from bson.objectid import ObjectId
from langdetect import detect, LangDetectException
from langchain_community.vectorstores import FAISS
from langchain.chains import RetrievalQA
from langchain_ollama import OllamaLLM
from langchain.schema import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_ollama import OllamaEmbeddings
import faiss
from functools import wraps
from data_ingestion import search_arxiv, search_pubmed,search_openalex,search_semantic_scholar,merge_and_deduplicate  # Assurez-vous que ce fichier est bien présent

from utils.faiss_index import delete_from_index
from PIL import Image

# --- Config Flask ---
app = Flask(__name__)
app.config.update({
    "SECRET_KEY": "votre_cle_super_secrete",
    "JWT_SECRET_KEY": "votre_cle_super_secrete",
    "MONGO_URI": "mongodb://localhost:27017/PFE",
})

BASE_DIR = os.path.dirname(__file__)
UPLOAD_FOLDER = os.path.join(BASE_DIR, "documents")
THUMBNAIL_FOLDER = os.path.join(BASE_DIR, "thumbnails")
DOWNLOAD_FOLDER = os.path.join(BASE_DIR, "downloads")
client = ollama.Client()
ALLOWED_EXTENSIONS = {"txt", "pdf"}
MODEL_NAME = "llama3"

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(DOWNLOAD_FOLDER, exist_ok=True)

# --- Extensions ---
CORS(app, supports_credentials=True)
socketio = SocketIO(app, cors_allowed_origins=["http://localhost:3000"], async_mode="eventlet")
jwt_mgr = JWTManager(app)
mongo = PyMongo(app)
bcrypt = Bcrypt(app)

auth_bp = Blueprint("auth", __name__)
protected_bp = Blueprint("protected", __name__)
connected_users = {}  # sid -> user_id mapping
faiss_index = None
documents = [] 

# --- Helpers JWT pour Socket.IO ---
def decode_jwt(token: str):
    try:
        return pyjwt.decode(token, app.config['SECRET_KEY'], algorithms=['HS256'])
    except (pyjwt.ExpiredSignatureError, pyjwt.InvalidTokenError):
        return None

# --- Utils ---
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def remove_think_blocks(text):
    return re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL)

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


# --- Décorateur token_required pour routes REST ---
def token_required(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        if request.method == 'OPTIONS':
            return '', 200 
        auth = request.headers.get('Authorization', None)
        if not auth or not auth.startswith('Bearer '):
            return jsonify({'message': 'Token manquant'}), 401
        token = auth.split()[1]
        payload = decode_jwt(token)
        if not payload or 'user_id' not in payload:
            return jsonify({'message': 'Token invalide ou expiré'}), 401
        user = mongo.db.users.find_one({'_id': ObjectId(payload['user_id'])})
        if not user:
            return jsonify({'message': 'Utilisateur introuvable'}), 404
        return f(user, *args, **kwargs)
    return wrapper



def embed_text(text):
    return np.random.rand(768).astype("float32")

def index_document(text, doc_id):
    global faiss_index, documents
    vec = embed_text(text)
    if faiss_index is None:
        faiss_index = faiss.IndexFlatL2(vec.shape[0])
    faiss_index.add(np.expand_dims(vec, axis=0))
    documents.append({"page_content": text, "doc_id": doc_id})

def query_faiss_index(question_embedding, top_k=5):
    if faiss_index is None or faiss_index.ntotal == 0:
        return []

    distances, indices = faiss_index.search(np.expand_dims(question_embedding, axis=0), top_k)

    results = []
    for dist, idx in zip(distances[0], indices[0]):
        if idx < len(documents):
            results.append((documents[idx], dist))
    return results

# --- Routes Auth REST ---
@auth_bp.route("/register", methods=["POST"])
@cross_origin(origin="localhost", headers=["Content-Type", "Authorization"])
def register():
    data = request.json or {}
    if not (data.get("email") and data.get("password")):
        return jsonify({"error": "Email et mot de passe requis"}), 400
    if mongo.db.users.find_one({"email": data["email"]}):
        return jsonify({"error": "Email déjà pris"}), 409
    hashed_password = bcrypt.generate_password_hash(data["password"]).decode()
    uid = mongo.db.users.insert_one({
        "email": data["email"],
        "name": data.get("name"),
        "password": hashed_password
    }).inserted_id
    return jsonify({"message": "Inscription réussie", "user_id": str(uid)}), 201

@auth_bp.route("/login", methods=["POST"])
@cross_origin(origin="localhost", headers=["Content-Type","Authorization"])
def login():
    d=request.json or {}
    u=mongo.db.users.find_one({"email":d.get("email")})
    if not u or not bcrypt.check_password_hash(u["password"],d.get("password","")):
        return jsonify({"msg":"Identifiants invalides"}),401
    tok=pyjwt.encode({
        "user_id":str(u["_id"]),
        "exp":datetime.utcnow()+timedelta(hours=12)
    }, app.config["SECRET_KEY"], algorithm="HS256")
    return jsonify({"token":tok,"user":{"id":str(u["_id"]),"email":u["email"],"name":u.get("name")}})


@protected_bp.route("/me", methods=["GET"])
@token_required
def me(user):
    return jsonify({
        "id": str(user["_id"]),
        "email": user["email"],
        "name": user.get("name")
    })

app.register_blueprint(auth_bp, url_prefix="/auth")
app.register_blueprint(protected_bp, url_prefix="/protected")

# --- SocketIO handlers ---
@socketio.on('connect')
def on_connect(auth):
    token = auth.get('token') if auth else None
    if not token:
        print("Connexion refusée : pas de token")
        return disconnect()
    payload = decode_jwt(token)
    if not payload or 'user_id' not in payload:
        print("Connexion refusée : token invalide")
        return disconnect()
    connected_users[request.sid] = payload['user_id']
    print(f"Utilisateur connecté : {payload['user_id']} (SID: {request.sid})")

@socketio.on('disconnect')
def on_disconnect():
    user_id = connected_users.pop(request.sid, None)
    print(f"Déconnexion: {user_id} (SID: {request.sid})")

@socketio.on("chat_message")
def handle_message(data):
    model = (data.get("model") or MODEL_NAME).strip()
    uid = connected_users.get(request.sid)
    now = datetime.utcnow().isoformat() + "Z"
    if not uid:
        emit("auth_failed", {"message": "❌ Auth requise"})
        return disconnect()

    user_input = (data.get("user_input") or "").strip()
    conv_id = data.get("conversation_id")
    if not user_input:
        return emit("stream_response", {"token": "⚠️ Message vide"})

    try:
        # ==== 🔍 Étape 1 : Chargement FAISS ====
        user_folder = os.path.join(UPLOAD_FOLDER, str(uid))
        index_folder = os.path.join(user_folder, "faiss_index")
        index_file = os.path.join(index_folder, "index.faiss")

        if not os.path.exists(index_file):
            emit("stream_response", {"token": "⚠️ Aucun index disponible pour cet utilisateur"})
            raise FileNotFoundError(f"Index non trouvé pour l'utilisateur {uid}")

        embedding = OllamaEmbeddings(model="nomic-embed-text")
        store = FAISS.load_local(index_folder, embedding, allow_dangerous_deserialization=True)
        print(f"[chat_message] 📁 Index FAISS chargé avec {len(store.docstore._dict)} documents.")

        # ==== 🔍 Étape 2 : Embedding de la question et recherche FAISS ====
        question_embedding = embedding.embed_query(user_input)
        docs_similaires = store.similarity_search_with_score_by_vector(question_embedding, top_k=5)
        print(f"[chat_message] 🔎 Résultats FAISS bruts (top 5) :")

        for i, (doc, dist) in enumerate(docs_similaires):
            print(f"  → Doc #{i+1} | Distance : {dist:.4f} | Extrait : {doc.page_content[:100]}...")

        # ==== 🔍 Étape 3 : Filtrage des documents pertinents ====
        threshold = 0.7
        relevant = [
            doc for doc, dist in docs_similaires
            if dist < threshold and len(doc.page_content.strip()) > 100
        ]
        print(f"[chat_message] ✅ {len(relevant)} documents pertinents trouvés avec seuil {threshold}")

          # Détection de langue sécurisée avec fallback
        try:
            detected_lang = detect(user_input)
        except Exception:
            detected_lang = "fr"  # Langue par défaut


        # ==== 🔁 Étape 4 : Si documents pertinents → RAG ====
        if relevant:
            print("[chat_message] 🤖 RAG activé : génération via RetrievalQA")

            # ⚠️ Message temporaire vers le frontend (indique que c'est lent mais normal)
            emit("stream_response", {"token": "⏳ Je réfléchis à ta question, un instant..."})

            system_prompt = f"Tu es un assistant IA utile. Réponds toujours en {detected_lang}."
            # Limiter à 2 documents pertinents max
            retriever = store.as_retriever(
                search_type="similarity",
                search_kwargs={"k": 2}
            )

            # 🧠 Création de la chaîne QA avec mode map_reduce (meilleur pour gros contextes)
            qa_chain = RetrievalQA.from_chain_type(
                llm=OllamaLLM(model=model,system=system_prompt,
        temperature=0.3),
                chain_type="map_reduce",  # ou "refine"
                retriever=retriever
            )

            # 🚀 Lancement
            result = qa_chain.invoke(user_input)
            response = remove_think_blocks(result["result"].strip())

            emit("stream_response", {"token": response})

            # 💾 Sauvegarde conversation
            if conv_id:
                print(f"[chat_message] 💾 Sauvegarde (RAG) dans conversation {conv_id}")
                token = data.get("token")
                headers = {"Authorization": f"Bearer {token}"}
                requests.post(
                    "http://localhost:5000/save_chat",
                    json={
                        "conversation_id": conv_id,
                        "messages": [
                            {
                                "role": "user",
                                "content": user_input,
                                "timestamp": now
                            },
                            {
                                "role": "assistant",
                                "content": buffer,
                                "timestamp": datetime.utcnow().isoformat() + "Z"
                            }
                        ]
                    },
                    headers=headers
                )
            return
        else:
            print("[chat_message] ⚠️ Aucun document suffisamment pertinent - fallback LLM activé")

    except Exception as e:
        print(f"[RAG] ❌ Erreur: {e}")
        emit("stream_response", {"token": f"⚠️ Erreur RAG: {e}"})

    # ==== 🧠 Fallback vers Ollama seul ====
    try:
        detected_lang = detect(user_input)
        system_prompt = f"Tu es un assistant IA. Réponds en {detected_lang}."
        print(f"[chat_message] 🔁 Fallback - génération directe via Ollama ({model})")

        stream = ollama.chat(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_input}
            ],
            stream=True
        )
        buffer = ""
        for chunk in stream:
            token = chunk["message"]["content"]
            buffer += token
        emit("stream_response", {"token": buffer})

        if conv_id:
            print(f"[chat_message] 💾 Sauvegarde (Fallback) dans conversation {conv_id}")
            token = data.get("token")
            headers = {"Authorization": f"Bearer {token}"}
            requests.post(
    "http://localhost:5000/save_chat",
    json={
        "conversation_id": conv_id,
        "messages": [
            {
                "role": "user",
                "content": user_input,
                "timestamp": now
            },
            {
                "role": "assistant",
                "content": buffer,
                "timestamp": datetime.utcnow().isoformat() + "Z"
            }
        ]
    },
    headers=headers
)

    except Exception as e:
        import traceback
        print(f"[LLM Fallback] ⚠️ Exception: {e}\n{traceback.format_exc()}")
        emit("stream_response", {"token": f"⚠️ Erreur LLM: {e}"})


    
@app.route("/save_chat", methods=["POST"])
@token_required
def save_chat(user):
    data = request.json or {}
    conv_id = data.get("conversation_id")
    messages = data.get("messages")

    if not (conv_id and messages and isinstance(messages, list)):
        return jsonify({"error": "Données manquantes ou incorrectes"}), 400

    conv = mongo.db.conversations.find_one({"_id": ObjectId(conv_id), "user_id": user["_id"]})
    if not conv:
        return jsonify({"error": "Conversation introuvable"}), 404

    # Ajouter chaque message individuellement dans la liste
    for msg in messages:
        mongo.db.conversations.update_one(
            {"_id": ObjectId(conv_id)},
            {"$push": {"messages": {
                "role": msg["role"],
                "content": msg["content"],
                "timestamp": msg.get("timestamp", datetime.utcnow())
            }}}
        )

    # Si c'est le premier message, générer un titre automatiquement
    if conv["title"] == "Nouvelle conversation" and len(conv.get("messages", [])) == 0:
        first_user_msg = next((m for m in messages if m["role"] == "user"), None)
        if first_user_msg:
            generated_title = first_user_msg["content"].strip().split("?")[0][:50]
            mongo.db.conversations.update_one(
                {"_id": ObjectId(conv_id)},
                {"$set": {"title": generated_title}}
            )

    return jsonify({"status": "ok"})

@app.route("/upload", methods=["POST"])
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

    # Initialiser la liste des nouveaux documents à indexer
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

            # Création thumbnail
            user_thumb_folder = os.path.join(THUMBNAIL_FOLDER, str(user["_id"]))
            os.makedirs(user_thumb_folder, exist_ok=True)
            thumb_path = os.path.join(user_thumb_folder, f"{filename}.png")
            create_thumbnail(filepath, thumb_path)

            # Enregistrement en base MongoDB (avec contenu complet)
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

            # Ajouter à la liste pour indexation
            new_doc = Document(page_content=content, metadata={"source": filepath})
            new_documents.append(new_doc)

    if not new_documents:
        return jsonify({"error": "Aucun fichier valide à indexer"}), 400

    # Chargement ou création de l'index FAISS
    embedding = OllamaEmbeddings(model="nomic-embed-text")

    index_path = os.path.join(user_folder, "faiss_index")

    try:
        store = FAISS.load_local(index_path, embedding=embedding, allow_dangerous_deserialization=True)
        store.add_documents(new_documents)
    except Exception:
        store = FAISS.from_documents(new_documents, embedding=embedding)

    # Sauvegarde de l’index FAISS
    store.save_local(index_path)

    return jsonify({"files": saved_files})

# --- Route d'accès aux miniatures ---
@app.route('/thumbnails/<user_id>/<filename>')
def serve_thumbnail(user_id, filename):
    thumb_path = os.path.join(THUMBNAIL_FOLDER, user_id)
    return send_from_directory(thumb_path, filename)


# --- Endpoint GET documents par conversation_id ---
@app.route('/user_documents/<conversation_id>', methods=['GET'])
@token_required
def get_user_documents(user, conversation_id):
    docs = list(mongo.db.documents.find({"conversation_id": conversation_id, "user_id": user["_id"]}))
    # Convertir ObjectId en string
    for doc in docs:
        doc["_id"] = str(doc["_id"])
        if isinstance(doc["user_id"], ObjectId):
            doc["user_id"] = str(doc["user_id"])
    return jsonify({"files": docs})

# --- Endpoint GET messages par conversation_id ---
@app.route("/conversation/<conversation_id>", methods=["GET"])
@token_required
def get_conversation_messages(user, conversation_id):
    conv = mongo.db.conversations.find_one({"_id": ObjectId(conversation_id), "user_id": user["_id"]})
    if not conv:
        return jsonify({"error": "Conversation non trouvée"}), 404

    # On renvoie les messages (supposé dans conv["messages"])
    return jsonify({"messages": conv.get("messages", [])})


# @app.route('/user_documents', methods=['GET'])
# @token_required
# def user_documents():
#     user_id = connected_users.get(request.sid)  # récupère l'ID ou identifiant de l'utilisateur depuis le token

#     # Dossier où sont stockés les fichiers de l'utilisateur
#     user_upload_folder = os.path.join(app.config['UPLOAD_FOLDER'], str(user_id))

#     files_list = []
#     if os.path.exists(user_upload_folder):
#         for filename in os.listdir(user_upload_folder):
#             filepath = os.path.join(user_upload_folder, filename)

#             # Ici tu peux adapter si tu as des thumbnails
#             thumbnail_path = None  
#             # Par exemple, si tu as un dossier thumbnails ou tu as un fichier thumbnail nommé filename + ".png"
#             thumbnail_candidate = os.path.join(user_upload_folder, "thumbnails", filename + ".png")
#             if os.path.exists(thumbnail_candidate):
#                 thumbnail_path = f"/uploads/{user_id}/thumbnails/{filename}.png"  # ou ton URL accessible

#             files_list.append({
#                 "filename": filename,
#                 "thumbnail": thumbnail_path
#             })

#     return jsonify(files=files_list)

@app.route("/delete_document/<doc_id>", methods=["DELETE"])
@token_required
def delete_document(user, doc_id):
    doc = mongo.db.documents.find_one({"_id": ObjectId(doc_id), "user_id": user["_id"]})
    if not doc:
        return jsonify({"error": "Document introuvable"}), 404

    # Supprimer le fichier
    if os.path.exists(doc["path"]):
        os.remove(doc["path"])

    # Supprimer la miniature
    thumb_path = os.path.join(THUMBNAIL_FOLDER, str(user["_id"]), doc["thumbnail"])
    if os.path.exists(thumb_path):
        os.remove(thumb_path)

    # Supprimer de l'index FAISS
    success = delete_from_index(doc["path"], user["_id"])

    # Supprimer de MongoDB
    mongo.db.documents.delete_one({"_id": ObjectId(doc_id)})

    return jsonify({"message": "Document supprimé", "faiss_deleted": success})

@app.route("/conversations",methods=["POST"])
@token_required
def create_conv(user):
    t=(request.json or {}).get("title","Nouvelle conversation")
    res=mongo.db.conversations.insert_one({"user_id":user["_id"],"title":t,"messages":[],"created_at":datetime.utcnow()})
    return jsonify({"conversation_id":str(res.inserted_id)}),201

@app.route("/conversations",methods=["GET"])
@token_required
def list_conv(user):
    convs=mongo.db.conversations.find({"user_id":user["_id"]})
    return jsonify([{"id":str(c["_id"]),"title":c["title"],"created_at":c["created_at"]} for c in convs])

@app.route("/conversations/<cid>",methods=["GET"])
@token_required
def get_conv(user,cid):
    conv=mongo.db.conversations.find_one({"_id":ObjectId(cid),"user_id":user["_id"]})
    if not conv:return jsonify({"error":"Introuvable"}),404
    return jsonify({"id":str(conv["_id"]),"title":conv["title"],"messages":conv["messages"]})

@app.route("/conversations/<cid>",methods=["PUT"])
@token_required
def rename_conv(user,cid):
    nt=(request.json or {}).get("title")
    if not nt:return jsonify({"error":"title manquant"}),400
    res=mongo.db.conversations.update_one({"_id":ObjectId(cid),"user_id":user["_id"]},{"$set":{"title":nt}})
    if res.matched_count==0:return jsonify({"error":"Introuvable"}),404
    return jsonify({"msg":"Title updated"})

@app.route("/conversations/<cid>", methods=["DELETE"])
@token_required
def delete_conv(user, cid):
    res = mongo.db.conversations.delete_one({"_id": ObjectId(cid), "user_id": user["_id"]})
    if res.deleted_count == 0:
        return jsonify({"error": "Introuvable"}), 404
    return jsonify({"msg": "Conversation supprimée"})


@app.route('/fetch_documents', methods=['POST'])
@token_required
def fetch_documents(user):
    data = request.get_json()
    query = data.get('query')
    source = data.get('source')
    max_results = data.get('max_results', 5)

    documents = []
    if source == 'all':
        all_docs = []
        all_docs.append(search_arxiv(query, max_results))
        all_docs.append(search_pubmed(query, max_results))
        all_docs.append(search_semantic_scholar(query, max_results))
        all_docs.append(search_openalex(query, max_results))
        documents = merge_and_deduplicate(all_docs)
        # etc.
    elif source == 'arxiv':
        documents = search_arxiv(query, max_results)
    elif source == 'pubmed':
        documents = search_pubmed(query, max_results)
    elif source == 'sematic_scholar':
        documents = search_pubmed(query, max_results)
    elif source == 'OpenAlex':
        documents = search_pubmed(query, max_results)
    else:
        # Optionnel: gérer source inconnue
        pass
    
    documents = sorted(documents, key=lambda d: d.get('title', '').lower())
    return jsonify({'documents': documents})

if __name__ == "__main__":
        socketio.run(
        app,
        host='0.0.0.0',
        port=5000,
        debug=True,
        use_reloader=False
    ) 