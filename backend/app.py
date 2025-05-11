import os
import re
import time
import traceback

from datetime import datetime, timedelta
import pdfplumber
import ollama
from flask import Flask, request, jsonify, send_from_directory, Blueprint, current_app
from flask_cors import CORS, cross_origin
from flask_socketio import SocketIO, emit
from werkzeug.utils import secure_filename
import requests
from data_ingestion import fetch_arxiv_papers, fetch_pubmed_articles
from auth import init_auth, token_required

from langchain_community.vectorstores import FAISS
from langchain.chains import RetrievalQA
from langchain_ollama import OllamaLLM
from langchain.schema import Document
from langchain.embeddings.base import Embeddings
from langchain.text_splitter import RecursiveCharacterTextSplitter

from flask_bcrypt import Bcrypt
from bson.objectid import ObjectId
from flask_pymongo import PyMongo
import jwt
from functools import wraps

UPLOAD_FOLDER = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'documents')
INDEX_FOLDER = os.path.join(UPLOAD_FOLDER, 'faiss_index')
DOWNLOAD_FOLDER = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'downloads')
ALLOWED_EXTENSIONS = {"txt", "pdf"}
MODEL_NAME = "deepseek-r1:7b"

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["DOWNLOAD_FOLDER"] = DOWNLOAD_FOLDER
app.config["SECRET_KEY"] = "votre_cle_super_secrete"
app.config["MONGO_URI"] = "mongodb://localhost:27017/PFE"

CORS(app, supports_credentials=True)
socketio = SocketIO(app, cors_allowed_origins=["http://localhost:3000"], async_mode='threading')

mongo = PyMongo(app)
bcrypt = Bcrypt(app)
auth_bp = Blueprint('auth', __name__)
protected_bp = Blueprint("protected", __name__)

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(DOWNLOAD_FOLDER, exist_ok=True)

# --- Utilities ---
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def remove_think_blocks(text):
    return re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL)

def extract_pdf_content(filepath):
    try:
        with pdfplumber.open(filepath) as pdf:
            return "\n".join(page.extract_text() or '' for page in pdf.pages)
    except Exception as e:
        raise Exception(f"Erreur lecture PDF : {e}")

def is_valid_email(email):
    email_regex = r'^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$'
    return re.match(email_regex, email) is not None

def is_valid_password(password):
    password_regex = r'^(?=.*[A-Za-z])(?=.*\d)'
    return re.match(password_regex, password) is not None

# --- Ollama Embeddings ---
class OllamaEmbeddings(Embeddings):
    def __init__(self, model: str = "nomic-embed-text"):
        self.model = model

    def embed_documents(self, texts: list[str]) -> list[list[float]]:
        return [ollama.embeddings(model=self.model, prompt=text)["embedding"] for text in texts]

    def embed_query(self, text: str) -> list[float]:
        return ollama.embeddings(model=self.model, prompt=text)["embedding"]

# --- Token Decorator ---
def token_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        token = None
        if 'Authorization' in request.headers:
            token = request.headers['Authorization'].split(" ")[1]

        if not token:
            return jsonify({'message': 'Token manquant !'}), 401

        try:
            data = jwt.decode(token, current_app.config["SECRET_KEY"], algorithms=["HS256"])
            current_user = mongo.db.users.find_one({"_id": ObjectId(data['user_id'])})
        except jwt.ExpiredSignatureError:
            return jsonify({'message': 'Le token a expiré !'}), 401
        except jwt.InvalidTokenError:
            return jsonify({'message': 'Token invalide !'}), 401

        return f(current_user, *args, **kwargs)
    return decorated_function

@protected_bp.route("/me", methods=["GET"])
@token_required
def get_profile(current_user):
    return jsonify({
        "name":current_user["name"],
        "email": current_user["email"],
        "id": str(current_user["_id"])
    })
    
# --- Auth ---
@auth_bp.route("/register", methods=["POST"])
@cross_origin(origin='localhost', headers=['Content-Type', 'Authorization'])
def register():
    data = request.json
    name = data.get("name")
    email = data.get("email")
    password = data.get("password")
    print(f"Nom: {name}, Email: {email}, Mot de passe: {password}")
    if not email or not password:
        return jsonify({"error": "L'email et le mot de passe sont requis"}), 400
    if not is_valid_email(email):
        return jsonify({"error": "Email invalide"}), 400
    if not is_valid_password(password):
        print("Mot de passe invalide")
        return jsonify({"error": "Mot de passe invalide"}), 400 
    if mongo.db.users.find_one({"email": email}):
        return jsonify({"error": "Email déjà utilisé"}), 409

    hashed_pw = bcrypt.generate_password_hash(password).decode('utf-8')
    user_id = mongo.db.users.insert_one({"name":name,"email": email, "password": hashed_pw}).inserted_id

    return jsonify({"message": "Inscription réussie", "user_id": str(user_id)}), 201

@auth_bp.route("/login", methods=["POST"])
@cross_origin(origin='localhost', headers=['Content-Type', 'Authorization'])
def login():
    data = request.json
    email = data.get("email")
    name = data.get("name")
    password = data.get("password")

    user = mongo.db.users.find_one({"email": email})
    if user and bcrypt.check_password_hash(user["password"], password):
        token = jwt.encode({
            "user_id": str(user["_id"]),
            "exp": datetime.utcnow() + timedelta(hours=12)
        }, current_app.config["SECRET_KEY"], algorithm="HS256")
        print("Token généré:", token)

        # Renvoi du token et des informations utilisateur
        return jsonify({
            "token": token,
            "user": {
                "id": str(user["_id"]),
                "email": user["email"], 
                "name": user["name"]  # Ajoutez d'autres informations utilisateur si nécessaire
            }
        })
    else :
        return jsonify({"msg": "Identifiants invalides"}), 401





# --- Upload & indexation ---
@app.route("/upload", methods=["POST"])
@token_required
def upload_file(current_user):
    file = request.files.get("file")

    if not file or file.filename == "":
        return jsonify({"error": "Aucun fichier sélectionné"}), 400
    if not allowed_file(file.filename):
        return jsonify({"error": "Type de fichier non autorisé"}), 400

    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
    file.save(filepath)

    try:
        if filename.endswith(".pdf"):
            content = extract_pdf_content(filepath)
        else:
            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

        splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
        documents = splitter.create_documents([content])

        embeddings = OllamaEmbeddings()
        vector_store = FAISS.from_documents(documents, embeddings, index_type="IndexFlatL2")
        user_index_folder = os.path.join(INDEX_FOLDER, str(current_user['_id']))
        os.makedirs(user_index_folder, exist_ok=True)
        vector_store.save_local(user_index_folder)

        mongo.db.documents.insert_one({
            "user_id": current_user['_id'],
            "filename": filename,
            "filepath": filepath,
            "uploaded_at": datetime.utcnow()
        })

        return jsonify({"message": "Indexation réussie."}), 200

    except Exception as e:
        print("Erreur lors de l'upload :", traceback.format_exc())
        return jsonify({"error": str(e)}), 500



# --- Téléchargements ---
@app.route("/uploads/<filename>")
@token_required
def uploaded_file(current_user, filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

@app.route("/downloads/<filename>")
@token_required
def download_document(current_user, filename):
    file_path = os.path.join(app.config['DOWNLOAD_FOLDER'], filename)
    if os.path.exists(file_path):
        return send_from_directory(app.config['DOWNLOAD_FOLDER'], filename)
    else:
        return jsonify({"error": "Fichier non trouvé"}), 404

# --- Documents scientifiques ---
@app.route('/fetch_documents', methods=['POST'])
@token_required
def fetch_documents(current_user):
    data = request.json
    query = data.get("query", "machine learning")
    max_results = int(data.get("max_results", 5))
    source = data.get("source", "arxiv")

    try:
        if source == "arxiv":
            documents = fetch_arxiv_papers(query, max_results)
        elif source == "pubmed":
            documents = fetch_pubmed_articles(query, max_results)
        elif source == "all":
            documents = fetch_arxiv_papers(query, max_results) + fetch_pubmed_articles(query, max_results)
        else:
            documents = []
        return jsonify({"documents": documents})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

# --- Run ---
app.register_blueprint(auth_bp, url_prefix="/auth")
app.register_blueprint(protected_bp, url_prefix="/protected")

@app.route("/save_chat", methods=["POST"])
@token_required
def save_chat(current_user):
    data = request.json
    mongo.db.history.insert_one({
        "user_id": current_user["_id"],
        "question": data.get("question"),
        "response": data.get("response"),
        "timestamp": datetime.utcnow()
    })
    return jsonify({"message": "Historique sauvegardé"})
@app.route("/api/history", methods=["GET"])
@token_required
def get_history(current_user):
    history_cursor = mongo.db.history.find({"user_id": current_user["_id"]}).sort("timestamp", -1)
    history = []
    for item in history_cursor:
        history.append({
            "question": item.get("question"),
            "response": item.get("response"),
            "timestamp": item.get("timestamp").isoformat() if item.get("timestamp") else None
        })
    return jsonify({"history": history})


# Créer une nouvelle conversation
@app.route('/conversations', methods=['POST'])
@token_required
def create_conversation(current_user):
    data = request.get_json()
    title = data.get('title', 'Nouvelle conversation')
    
    # Créer une nouvelle conversation
    new_conv = {
        'user_id': current_user['_id'],
        'title': title,
        'messages': [],
        'created_at': datetime.utcnow()
    }
    
    # Insérer la conversation dans la base de données
    result = mongo.db.conversations.insert_one(new_conv)
    return jsonify({"conversation_id": str(result.inserted_id)}), 201

# Ajouter un message à une conversation
@app.route('/conversations/<conversation_id>/messages', methods=['POST'])
@token_required
def add_message(current_user, conversation_id):
    data = request.get_json()
    
    # Vérifier les données nécessaires
    if "role" not in data or "content" not in data:
        return jsonify({"error": "Données manquantes (role et content)"}), 400

    message = {
        "role": data["role"],
        "content": data["content"],
        "timestamp": datetime.utcnow()
    }
    
    # Ajouter le message à la conversation
    mongo.db.conversations.update_one(
        {"_id": ObjectId(conversation_id), "user_id": current_user["_id"]},
        {"$push": {"messages": message}}
    )
    return jsonify({"msg": "Message ajouté"}), 200

# Récupérer toutes les conversations de l'utilisateur
@app.route('/conversations', methods=['GET'])
@token_required
def get_conversations(current_user):
    conversations = mongo.db.conversations.find({"user_id": current_user["_id"]})
    
    # Renvoi des conversations avec pagination optionnelle (par exemple)
    return jsonify([{
        "id": str(conv["_id"]),
        "title": conv["title"],
        "created_at": conv["created_at"]
    } for conv in conversations]), 200

# Récupérer les messages d'une conversation spécifique
@app.route('/conversations/<conversation_id>', methods=['GET'])
@token_required
def get_conversation(current_user, conversation_id):
    conv = mongo.db.conversations.find_one({
        "_id": ObjectId(conversation_id),
        "user_id": current_user["_id"]
    })
    
    if not conv:
        return jsonify({"error": "Conversation introuvable"}), 404
    
    return jsonify({
        "id": str(conv["_id"]),
        "title": conv["title"],
        "messages": conv["messages"]
    }), 200

@socketio.on("message")
def handle_message(data):
    user_input = data.get("user_input", "")
    print(f"Received user input: {user_input}")

    if not user_input:
        return

    try:
        if os.path.exists(INDEX_FOLDER):
            vector_store = FAISS.load_local(
                INDEX_FOLDER,
                embeddings=OllamaEmbeddings(),
                allow_dangerous_deserialization=True
            )

            if vector_store.index.ntotal > 0:
                retriever = vector_store.as_retriever(search_type="similarity", search_kwargs={"k": 3})
                docs_with_scores = vector_store.similarity_search_with_score(user_input, k=3)

                relevant_docs = [
                    doc for doc, score in docs_with_scores
                    if score >= 0.75 and len(doc.page_content.strip()) > 100
                ]

                if relevant_docs:
                    qa_chain = RetrievalQA.from_chain_type(
                        llm=OllamaLLM(model=MODEL_NAME),
                        chain_type="stuff",
                        retriever=retriever
                    )
                    result = qa_chain.invoke(user_input)
                    response = remove_think_blocks(result["result"].strip())
                    save_chat_data = {
                        "question": user_input,
                        "response": response
                    }
                    requests.post("http://localhost:5000/save_chat", json=save_chat_data, headers={
                        'Authorization': f"Bearer {data['token']}"
                    })
                    for token in response.split():
                        time.sleep(0.03)
                        emit("stream_response", {"token": token})
                    return

    except Exception as e:
        emit("stream_response", {"token": f"⚠️ Erreur index : {e}"})

    # Fallback direct LLM
    emit("stream_response", {"token": "🤖 Réponse générée sans document indexé.\n"})
    try:
        buffer = ""
        stream = ollama.chat(model=MODEL_NAME, messages=[{"role": "user", "content": user_input}], stream=True)
        for chunk in stream:
            token = chunk["message"]["content"]
        time.sleep(0.1)
        buffer += token 
        response = buffer
        save_chat_data = {
                        "question": user_input,
                        "response": response
                    }
        token = data.get("token")
        if not token:
            emit("stream_response", {"token": "⚠️ Erreur : jeton d'authentification manquant."})
            return

        requests.post("http://localhost:5000/save_chat", json=save_chat_data, headers={
            'Authorization': f"Bearer {token}"
        })
        emit("stream_response", {"token": token})
    except Exception as e:
        emit("stream_response", {"token": f"⚠️ Erreur LLM : {e}"})
        return jsonify({"error": f"Erreur LLM : {e}"}), 500
if __name__ == '__main__':
    socketio.run(app, debug=True)
