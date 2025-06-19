import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import re
import eventlet
eventlet.monkey_patch()
import traceback
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
import logging
from data_ingestion import search_arxiv, search_pubmed, search_openalex, merge_and_deduplicate
from utils import mongo, token_required
from utils.faiss_index import delete_from_index,update_faiss_index,CachedEmbeddingModel
from routes import chat_bp, bp_documents, bp_conversations,bp_multi,bp_rating,bp_protected
from routes.conversations_routes import get_conversation_history
from utils.auth import auth_bp
from utils.auth import init_auth
from utils.paths import  get_user_index_folder
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from bson import ObjectId
from rank_bm25 import BM25Okapi

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


# --- Extensions ---
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}}, supports_credentials=True)

socketio = SocketIO(app, cors_allowed_origins=["http://localhost:3000"],async_mode="eventlet")
jwt_mgr = JWTManager(app)
mongo.init_app(app)
bcrypt = Bcrypt(app)
init_auth(app)
# --- Blueprints ---
app.register_blueprint(chat_bp)
app.register_blueprint(bp_documents)
app.register_blueprint(bp_conversations)
app.register_blueprint(bp_protected, url_prefix="/protected")
app.register_blueprint(auth_bp, url_prefix='/auth')
app.register_blueprint(bp_multi, url_prefix='/multi')
app.register_blueprint(bp_rating, url_prefix='/rating') 
# --- Variables globales ---
connected_users = {}
faiss_index = None
documents = []

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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

def decode_jwt(token):
    try:
        return pyjwt.decode(token, app.config["SECRET_KEY"], algorithms=["HS256"])
    except Exception:
        return None
def cosine_similarity(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def find_validated_answer_similar(question, threshold=0.8):
    embedding_model = OllamaEmbeddings(model="nomic-embed-text")
    question_embedding = embedding_model.embed_query(question)

    docs = list(mongo.db.validated_answers.find({"validated": True, "embedding": {"$exists": True}}))

    best_match = None
    best_score = -1

    for doc in docs:
        try:
            doc_embedding = doc["embedding"]
            # Optionnel : convertir en numpy array si besoin
            # doc_embedding = np.array(doc_embedding)

            score = cosine_similarity(question_embedding, doc_embedding)
            if score > best_score:
                best_score = score
                best_match = doc
        except Exception as e:
            print(f"Erreur sur doc {doc.get('_id')}: {e}")
            continue

    if best_score >= threshold and best_match:
        return {"found": True, "answer": best_match["answer"], "similarity": best_score}
    else:
        return {"found": False, "similarity": best_score}

def preprocess_text(text):
    """Nettoyage de texte : tokenisation + suppression des stopwords"""
    tokens = word_tokenize(text.lower())
    stop_words = set(stopwords.words('english'))
    return [word for word in tokens if word.isalnum() and word not in stop_words]

def hybrid_search(query, store,conv_id,uid, top_k=5):
    
    if store is None:
        return "", []
    """Recherche hybride : FAISS + BM25 avec logs pour debug"""
    embedding_model = CachedEmbeddingModel(uid)
    query_vector = embedding_model.embed_query(query)

    # Recherche sémantique (FAISS)
    semantic_results = store.similarity_search_with_score_by_vector(query_vector, top_k=top_k)
    # Recherche lexicale (BM25)
    documents = list(mongo.db.documents.find({
        "user_id": ObjectId(uid),
        "conversation_id": conv_id
    }))
    if not documents:
        print("❌ Aucun document trouvé en base MongoDB.")
        return "⚠️ Aucun document disponible pour la recherche.", []

    contents = [doc['content'] for doc in documents]
    tokenized_corpus = [preprocess_text(text) for text in contents]
    bm25 = BM25Okapi(tokenized_corpus)
    bm25_scores = bm25.get_scores(preprocess_text(query))

    # Top résultats BM25
    top_bm25 = np.argsort(bm25_scores)[::-1][:top_k]
    bm25_results = [{
        "content": documents[i]['content'],
        "source": documents[i].get('filename', 'Document'),
        "score": bm25_scores[i],
        "type": "bm25"
    } for i in top_bm25]

    results = []

    # ✅ MODIF : Ne filtre pas par score FAISS (uniquement longueur)
    for doc, score in semantic_results:
        if len(doc.page_content.strip()) > 100:
            results.append({
                "content": doc.page_content,
                "source": doc.metadata.get("filename", "Document"),
                "score": score,
                "type": "semantic"
            })

    # ✅ MODIF : Tolérance plus grande pour BM25
    for r in bm25_results:
        if r['score'] > 0.2 and not any(r['content'] == res['content'] for res in results):
            results.append(r)

    if not results:
        print("❌ Aucun résultat trouvé par FAISS ni BM25.")
        return "⚠️ Aucun résultat pertinent trouvé dans vos documents.", []

    # ✅ Tri mixte selon le type
    results.sort(key=lambda r: r['score'] if r['type'] == 'semantic' else -r['score'])

    # Génération du contexte
    context_str = "DOCUMENT CONTEXT (Use only if directly relevant to the question):\n"
    sources_used = []
    for i, r in enumerate(results[:3]):
        context_str += f"\n--- Source {i+1} ({r['type']}) ---\n{r['content'][:1000]}\n"
        sources_used.append({"source": r["source"], "type": r["type"]})

    print(f"✅ CONTEXT GÉNÉRÉ POUR LA QUESTION : {query}\n{context_str}")
    return context_str, sources_used


def load_faiss_index_for_user(uid, conv_id):
    try:
        user_folder = os.path.join(UPLOAD_FOLDER, str(uid))
        index_folder = os.path.join(user_folder, "conversations", str(conv_id))
        index_file = os.path.join(index_folder, "index.faiss")

        if not os.path.exists(index_file):
            raise FileNotFoundError(f"Index non trouvé pour l'utilisateur {uid}")

        embedding_model = CachedEmbeddingModel(uid)
        return FAISS.load_local(index_folder, embedding_model, allow_dangerous_deserialization=True)
    except Exception as e:
        logger.error(f"FAISS load error: {str(e)}")
        return None



# --- Socket.IO handlers ---
@socketio.on('connect')
def on_connect(auth):
    print("🔌 Tentative de connexion socket.io...")
    print("Auth reçu :", auth)

    token = auth.get('token') if auth else None
    if not token:
        print("❌ Token manquant")
        return disconnect()

    payload = decode_jwt(token)
    if not payload or 'user_id' not in payload:
        print("❌ Token invalide")
        return disconnect()

    connected_users[request.sid] = payload['user_id']
    print(f"[connect] ✅ Utilisateur connecté : {payload['user_id']}")

@socketio.on('disconnect')
def on_disconnect():
    user_id = connected_users.pop(request.sid, None)
    print(f"[disconnect] 🔌 Déconnexion : {user_id}")
    
@socketio.on("chat_message")
def handle_message(data):
    import traceback
    model_name = (data.get("model") or MODEL_NAME).strip()
    uid = connected_users.get(request.sid)
    now = datetime.utcnow().isoformat() + "Z"
    user_input = (data.get("user_input") or "").strip()
    conv_id = data.get("conversation_id")
    force_llm = data.get("force_llm", False)  # Flag pour forcer génération directe
   
    if not uid:
        emit("auth_failed", {"message": "❌ Auth requise"})
        return disconnect()

    if not user_input:
        emit("stream_response", {"token": "⚠️ Message vide"})
        return

    # Détection langue avec fallback
    try:
        detected_lang, conf = langid.classify(user_input)
        if conf < 0:
            detected_lang = "fr"
    except Exception:
        detected_lang = "fr"

    lang_map = {
        "en": "english",
        "fr": "french",
        "es": "spanish",
        "de": "german",
    }
    lang_name = lang_map.get(detected_lang, "french")

    if lang_name == "english":
        prompt_suffix = (
            "When providing code examples, ensure they are complete, clear, "
            "and well formatted using the appropriate Markdown code blocks "
            "(```python`, ```html`, etc.). Use realistic and concrete examples."
        )
    else:
        prompt_suffix = (
            "Lorsque tu donnes des exemples de code, assure-toi qu'ils soient complets, "
            "clairs, et bien formatés en utilisant les blocs Markdown appropriés "
            "(```python`, ```html`, etc.). Utilise des exemples concrets et réalistes."
        )

    PERSONALITY_TEMPLATES = {
    "formelle": "Tu es un assistant professionnel et rigoureux, tu donnes des réponses précises et bien structurées.",
    "amicale": "Tu es un assistant chaleureux et sympathique, tu expliques les choses simplement avec un ton bienveillant.",
    "concise": "Tu es un assistant qui va droit au but. Tu réponds de façon brève et claire."
    }
    CONTEXT_TEMPLATES = {
    "recherche scientifique": "Tu aides à analyser, comprendre et résumer des documents scientifiques.",
    "juridique": "Tu agis comme un assistant juridique, tu aides à comprendre les lois, jurisprudences et contrats.",
    "général": "Tu es un assistant polyvalent pour toutes sortes de tâches."
    }

    personality = PERSONALITY_TEMPLATES.get(data.get("personality", "formelle"))
    context = CONTEXT_TEMPLATES.get(data.get("context", "général"))
    # Fonction utilitaire pour formater l'historique
    def format_history(history):
        formatted = ""
        for msg in history:
            role = "Utilisateur" if msg["role"] == "user" else "Assistant"
            formatted += f"{role} : {msg['content']}\n"
        return formatted

    try:
        # Recherche réponse validée similaire (optionnel)
        result = find_validated_answer_similar(user_input)
        if result["found"]:
            print(f"[handle_message] Réponse validée trouvée avec similarité {result['similarity']:.2f}")
            answer_str = str(result["answer"])
            for token in answer_str.split():
                emit("stream_response", {"token": token + " "}, room=request.sid, namespace='/')
                socketio.sleep(0.01)  # pour lisser le flux

            socketio.emit("stream_end", room=request.sid, namespace='/')
            return

    except Exception as e:
        print(f"[handle_message] Erreur recherche validée: {e}")

    # Instruction commune pour gérer continuité ou nouvelle question
    continuity_instruction = (
        "Voici un extrait de conversation précédente entre l'utilisateur et toi. "
        "Si la nouvelle question est indépendante ou sans lien clair, réponds directement sans utiliser l'historique. "
        "Mais si la nouvelle question dépend du contexte précédent, utilise-le pour construire ta réponse.\n"
        "N'indique pas cette instruction dans ta réponse."
    )

    # --- Génération directe forcée ---
    if force_llm:
        try:
            history = get_conversation_history(conv_id, k=5) if conv_id else []
            formatted_history = format_history(history)
            system_prompt = (
                f"{personality}\n{context}\n\n"
                f"Tu es un assistant IA. Tu réponds toujours en {lang_name}.\n{prompt_suffix}\n\n"
                f"[INSTRUCTION IMPORTANTE]\n{continuity_instruction}\n[FIN INSTRUCTION]\n\n"
                f"[HISTORIQUE DE LA CONVERSATION]\n{formatted_history}[FIN HISTORIQUE]\n"
            )
            print(f"[chat_message] 🔁 Génération directe (force_llm) via ({model_name})")
            stream = ollama.chat(
                model=model_name,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_input},
                ],
                stream=True,
            )
            buffer = ""
            for chunk in stream:
                token = chunk["message"]["content"]
                buffer += token
                emit("stream_response", {"token": token})
            if conv_id:
                token = data.get("token")
                headers = {"Authorization": f"Bearer {token}"}
                requests.post(
                    "http://localhost:5000/save_chat",
                    json={
                        "conversation_id": conv_id,
                        "messages": [
                            {"role": "user", "content": user_input, "timestamp": now},
                            {"role": "assistant", "content": buffer, "timestamp": datetime.utcnow().isoformat() + "Z"},
                        ],
                    },
                    headers=headers,
                )
            return
        except Exception as e:
            print(f"[LLM force_llm] ⚠️ Exception: {e}\n{traceback.format_exc()}")
            emit("stream_response", {"token": f"⚠️ Erreur LLM: {e}"})
            return

    # --- Logique RAG classique ---
    try:
        store = load_faiss_index_for_user(uid, conv_id)
        if store is None:
            raise ValueError("NoIndex")

        context_str, combined_docs = hybrid_search(user_input, store, conv_id, uid, top_k=5)
        print("[chat_message] CONTEXTE généré :", context_str[:200])

        history = get_conversation_history(conv_id, k=5) if conv_id else []
        formatted_history = format_history(history)

        if combined_docs:
            system_prompt = (
                f"{personality}\n{context}\n\n"
                f"Tu es un assistant IA utile et intelligent. "
                f"Tu réponds toujours en {lang_name}. {prompt_suffix}\n\n"
                f"[INSTRUCTION IMPORTANTE]\n{continuity_instruction}\n[FIN INSTRUCTION]\n\n"
                f"[CONTEXTE DOCUMENTAIRE]\n{context_str}\n[FIN CONTEXTE]\n\n"
                f"[CONVERSATION PRÉCÉDENTE]\n{formatted_history}[FIN CONVERSATION]\n\n"
                f"Nouvelle question : {user_input}\n"
                f"Réponse :"
            )

            stream = ollama.chat(
                model=model_name,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_input},
                ],
                stream=True,
            )
            buffer = ""
            for chunk in stream:
                token = chunk["message"]["content"]
                buffer += token
                emit("stream_response", {"token": token})

            if conv_id:
                token = data.get("token")
                headers = {"Authorization": f"Bearer {token}"}
                requests.post(
                    "http://localhost:5000/save_chat",
                    json={
                        "conversation_id": conv_id,
                        "messages": [
                            {"role": "user", "content": user_input, "timestamp": now},
                            {"role": "assistant", "content": buffer, "timestamp": datetime.utcnow().isoformat() + "Z"},
                        ],
                    },
                    headers=headers,
                )
            return
        else:
            print("[chat_message] ⚠️ Aucun document pertinent trouvé - fallback LLM activé")

    except Exception as e:
        if str(e) == "NoIndex":
            # Ne rien envoyer dans stream_response, juste passer au fallback silencieusement
            pass
        else:
            # Pour les autres erreurs, afficher un message
            emit("stream_response", {"token": f"⚠️ Erreur RAG: {e}"})
    try:
        def log_conversation_history(history):
            formatted_history = ""
            for i, msg in enumerate(history):
                role = "Utilisateur" if msg["role"] == "user" else "Assistant"
                content_preview = msg["content"][:100].replace("\n", " ")  # limite à 100 caractères
                formatted_history += f"{i+1}. {role} : {content_preview}\n"
            print("[DEBUG] Historique formaté pour prompt :\n" + formatted_history)
            return formatted_history

        history = get_conversation_history(conv_id, k=5)
        print(f"[DEBUG] Raw history for conv_id={conv_id}: {history}")
        formatted_history = log_conversation_history(history)
        system_prompt = (
            f"{personality}\n{context}\n\n"
            f"Tu es un assistant IA. Tu réponds toujours en {lang_name}.\n{prompt_suffix}\n\n"
            f"[INSTRUCTION IMPORTANTE]\n{continuity_instruction}\n[FIN INSTRUCTION]\n\n"
            f"[HISTORIQUE DE LA CONVERSATION]\n{formatted_history}[FIN HISTORIQUE]\n"
        )
        print(f"[chat_message] 🔁 Fallback - génération directe via ({model_name})")
        stream = ollama.chat(
            model=model_name,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_input},
            ],
            stream=True,
        )
        buffer = ""
        for chunk in stream:
            token = chunk["message"]["content"]
            buffer += token
            emit("stream_response", {"token": token})

        if conv_id:
            token = data.get("token")
            headers = {"Authorization": f"Bearer {token}"}
            requests.post(
                "http://localhost:5000/save_chat",
                json={
                    "conversation_id": conv_id,
                    "messages": [
                        {"role": "user", "content": user_input, "timestamp": now},
                        {"role": "assistant", "content": buffer, "timestamp": datetime.utcnow().isoformat() + "Z"},
                    ],
                },
                headers=headers,
            )

    except Exception as e:
        print(f"[LLM Fallback] ⚠️ Exception: {e}\n{traceback.format_exc()}")
        emit("stream_response", {"token": f"⚠️ Erreur LLM: {e}"})
        socketio.emit("stream_end", room=request.sid)






    
if __name__ == "__main__":
        socketio.run(
        app,
        host='0.0.0.0',
        port=5000,
        debug=True,
        use_reloader=False
    ) 