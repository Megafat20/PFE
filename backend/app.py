import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import re
import eventlet
eventlet.monkey_patch()

import requests
from datetime import datetime, timedelta
import numpy as np
import ollama
import jwt as pyjwt
from flask import Flask, request, jsonify, Blueprint
from flask_cors import CORS, cross_origin
from flask_socketio import SocketIO, emit, disconnect
from flask_jwt_extended import JWTManager,decode_token, get_jwt_identity
from flask_bcrypt import Bcrypt
from werkzeug.utils import secure_filename
from langdetect import detect
from langchain_community.vectorstores import FAISS
from langchain.chains import RetrievalQA
from langchain_ollama import OllamaLLM, OllamaEmbeddings
from langchain.schema import Document
import faiss
import langid
# Modules locaux
from data_ingestion import search_arxiv, search_pubmed, search_openalex, merge_and_deduplicate
from utils import mongo, token_required
from utils.faiss_index import delete_from_index
from routes import chat_bp, bp_documents, bp_conversations,bp_multi
from utils.auth import auth_bp
from utils.auth import init_auth
from utils.protected import protected_bp
# --- Config Flask ---
app = Flask(__name__)
app.config.update({
    "SECRET_KEY": "votre_cle_super_secrete",
    "JWT_SECRET_KEY": "votre_cle_super_secrete",
    "MONGO_URI": "mongodb://localhost:27017/PFE",
})

# --- Dossiers ---
BASE_DIR = os.path.dirname(__file__)
UPLOAD_FOLDER = os.path.join(BASE_DIR, "documents")
THUMBNAIL_FOLDER = os.path.join(BASE_DIR, "thumbnails")
DOWNLOAD_FOLDER = os.path.join(BASE_DIR, "downloads")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(DOWNLOAD_FOLDER, exist_ok=True)

ALLOWED_EXTENSIONS = {"txt", "pdf"}
MODEL_NAME = "llama3"
client = ollama.Client()

# --- Extensions ---
CORS(app, supports_credentials=True)
socketio = SocketIO(app, cors_allowed_origins=["http://localhost:3000"], async_mode="eventlet")
jwt_mgr = JWTManager(app)
mongo.init_app(app)
bcrypt = Bcrypt(app)
init_auth(app)
# --- Blueprints ---
app.register_blueprint(chat_bp)
app.register_blueprint(bp_documents)
app.register_blueprint(bp_conversations)
app.register_blueprint(protected_bp, url_prefix="/user")
app.register_blueprint(auth_bp, url_prefix='/auth')
app.register_blueprint(bp_multi, url_prefix='/multi')
# --- Variables globales ---
connected_users = {}
faiss_index = None
documents = []

# --- Utils ---
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def remove_think_blocks(text):
    return re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL)

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
    return [(documents[idx], dist) for dist, idx in zip(distances[0], indices[0]) if idx < len(documents)]

def load_faiss_index_for_user(uid):
    user_folder = os.path.join(UPLOAD_FOLDER, str(uid))
    index_folder = os.path.join(user_folder, "faiss_index")
    index_file = os.path.join(index_folder, "index.faiss")
    if not os.path.exists(index_file):
        raise FileNotFoundError("Index FAISS non trouvé")
    embedding = OllamaEmbeddings(model="nomic-embed-text")
    return FAISS.load_local(index_folder, embedding, allow_dangerous_deserialization=True)

def decode_jwt(token):
    try:
        return pyjwt.decode(token, app.config["SECRET_KEY"], algorithms=["HS256"])
    except Exception:
        return None



# --- Socket.IO handlers ---
@socketio.on('connect')
def on_connect(auth):
    token = auth.get('token') if auth else None
    if not token:
        return disconnect()
    payload = decode_jwt(token)
    if not payload or 'user_id' not in payload:
        return disconnect()
    connected_users[request.sid] = payload['user_id']
    print(f"[connect] ✅ Utilisateur connecté : {payload['user_id']}")

@socketio.on('disconnect')
def on_disconnect():
    user_id = connected_users.pop(request.sid, None)
    print(f"[disconnect] 🔌 Déconnexion : {user_id}")
@socketio.on("chat_message")
def handle_message(data):
    model_name = (data.get("model") or MODEL_NAME).strip()
    uid = connected_users.get(request.sid)
    now = datetime.utcnow().isoformat() + "Z"

    if not uid:
        emit("auth_failed", {"message": "❌ Auth requise"})
        return disconnect()

    user_input = (data.get("user_input") or "").strip()
    conv_id = data.get("conversation_id")

    if not user_input:
        return emit("stream_response", {"token": "⚠️ Message vide"})

    lang_name = "français"  # fallback de base
    buffer = ""

    try:
        # ==== 🔍 Étape 1 : Chargement FAISS ====
        user_folder = os.path.join(UPLOAD_FOLDER, str(uid))
        index_folder = os.path.join(user_folder, "faiss_index")
        index_file = os.path.join(index_folder, "index.faiss")

        if not os.path.exists(index_file):
            emit("stream_response", {"token": "⚠️ Aucun index disponible pour cet utilisateur"})
            raise FileNotFoundError(f"Index non trouvé pour l'utilisateur {uid}")

        embedding = OllamaEmbeddings(model="nomic-embed-text")
        store = load_faiss_index_for_user(uid)
        question_embedding = embedding.embed_query(user_input)
        docs_similaires = store.similarity_search_with_score_by_vector(question_embedding, top_k=5)

        print(f"[chat_message] 🔎 Résultats FAISS bruts (top 5) :")
        for i, (doc, dist) in enumerate(docs_similaires):
            print(f"  → Doc #{i+1} | Distance : {float(dist):.4f} | Extrait : {doc.page_content[:100]}...")

        # ==== 🔍 Étape 2 : Filtrage des documents pertinents ====
        threshold = 0.7
        relevant = [
            doc for doc, dist in docs_similaires
            if float(dist) < threshold and len(doc.page_content.strip()) > 100
        ]
        print(f"[chat_message] ✅ {len(relevant)} documents pertinents trouvés avec seuil {threshold}")

        # ==== 🌐 Détection de langue ====
        try:
            detected_lang, conf = langid.classify(user_input)
            print(f"Langue détectée: {detected_lang} (confiance: {conf:.2f})")
        except Exception:
            detected_lang = "fr"

        lang_map = {
            "en": "anglais",
            "fr": "français",
            "es": "espagnol",
            "al": "allemand",
        }
        lang_name = lang_map.get(detected_lang, "français")

        # ==== 🔁 Étape 3 : RAG ====
        if relevant:
            print("[chat_message] 🤖 RAG activé : génération via RetrievalQA")
            emit("stream_response", {"token": "⏳ Je réfléchis à ta question, un instant..."})

            system_prompt = f"Tu es un assistant IA utile. Réponds toujours en {lang_name}."
            retriever = store.as_retriever(search_type="similarity", search_kwargs={"k": 2})

            qa_chain = RetrievalQA.from_chain_type(
                llm=OllamaLLM(model=model_name, system=system_prompt, temperature=0.3),
                chain_type="map_reduce",
                retriever=retriever
            )

            result = qa_chain.invoke(user_input)
            response = remove_think_blocks(result["result"].strip())
            buffer = response
            emit("stream_response", {"token": response})

            if conv_id:
                print(f"[chat_message] 💾 Sauvegarde (RAG) dans conversation {conv_id}")
                token = data.get("token")
                headers = {"Authorization": f"Bearer {token}"}
                requests.post(
                    "http://localhost:5000/save_chat",
                    json={
                        "conversation_id": conv_id,
                        "messages": [
                            {"role": "user", "content": user_input, "timestamp": now},
                            {"role": "assistant", "content": buffer, "timestamp": datetime.utcnow().isoformat() + "Z"}
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

    # ==== 🧠 Fallback vers Ollama ====
    try:
        system_prompt = f"Tu es un assistant IA. Réponds en {lang_name}."
        print(f"[chat_message] 🔁 Fallback - génération directe via ({model_name})")
        stream = ollama.chat(
            model=model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_input}
            ],
            stream=True
        )
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
                        {"role": "user", "content": user_input, "timestamp": now},
                        {"role": "assistant", "content": buffer, "timestamp": datetime.utcnow().isoformat() + "Z"}
                    ]
                },
                headers=headers
            )

    except Exception as e:
        import traceback
        print(f"[LLM Fallback] ⚠️ Exception: {e}\n{traceback.format_exc()}")
        emit("stream_response", {"token": f"⚠️ Erreur LLM: {e}"})


    
if __name__ == "__main__":
        socketio.run(
        app,
        host='0.0.0.0',
        port=5000,
        debug=True,
        use_reloader=False
    ) 