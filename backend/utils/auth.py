import os
from flask import Blueprint, request, jsonify, current_app
from extensions import mongo
from flask_bcrypt import Bcrypt
import jwt
from datetime import datetime, timedelta
from functools import wraps
from bson.objectid import ObjectId
from langchain_community.vectorstores import FAISS
from langchain.schema import Document
from langchain_ollama import OllamaEmbeddings
from .jwt_utils import decode_jwt


auth_bp = Blueprint('auth', __name__)
bcrypt = Bcrypt()


@auth_bp.route("/register", methods=["POST"])
def register():
    mongo = current_app.mongo
    data = request.json
    email = data.get("email")
    name = data.get("name")
    password = data.get("password")

    if not email or not password:
        return jsonify({"error": "Champs requis"}), 400

    if mongo.db.users.find_one({"email": email}):
        return jsonify({"error": "Utilisateur existe déjà"}), 409

    hashed_pw = bcrypt.generate_password_hash(password).decode('utf-8')
    result = mongo.db.users.insert_one({
        "name": name,
        "email": email,
        "password": hashed_pw
    })
    user_id = str(result.inserted_id)
    init_user_index(user_id)
    return jsonify({"message": "Inscription réussie"}), 201

@auth_bp.route("/login", methods=["POST"])
def login():
    mongo = current_app.mongo
    data = request.json
    email = data.get("email")
    password = data.get("password")

    user = mongo.db.users.find_one({"email": email})
    if user and bcrypt.check_password_hash(user["password"], password):
        token = jwt.encode({
            "user_id": str(user["_id"]),
            "exp": datetime.utcnow() + timedelta(hours=12)
        }, current_app.config["SECRET_KEY"], algorithm="HS256")
        print("Token généré:", token)
        return jsonify({"token": token})
    return jsonify({"error": "Identifiants invalides"}), 401



def token_required(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        
        if request.method == "OPTIONS":
            # simple réponse pour le preflight
            return '', 200
        print("token_required called")
        auth = request.headers.get('Authorization', None)
        print("Authorization header:", auth)
        if not auth or not auth.startswith('Bearer '):
            print("Token manquant")
            return jsonify({'message': 'Token manquant'}), 401
        token = auth.split()[1]
        payload = decode_jwt(token)
        print("Payload decoded:", payload)
        if not payload or 'user_id' not in payload:
            print("Token invalide ou expiré")
            return jsonify({'message': 'Token invalide ou expiré'}), 401
        user = mongo.db.users.find_one({'_id': ObjectId(payload['user_id'])})
        if not user:
            print("Utilisateur introuvable")
            return jsonify({'message': 'Utilisateur introuvable'}), 404
        return f(user, *args, **kwargs)
    return wrapper

@auth_bp.route("/update_profile", methods=["PUT", "OPTIONS"])
@token_required
def update_profile(current_user):
    if request.method == "OPTIONS":
        # simple réponse pour le preflight
        return '', 200
    mongo = current_app.mongo
    data = request.json

    # Champs que l'utilisateur peut mettre à jour
    allowed_fields = ["name", "persona", "context"]
    update_fields = {field: data[field] for field in allowed_fields if field in data}

    if not update_fields:
        return jsonify({"error": "Aucun champ à mettre à jour"}), 400

    mongo.db.users_profile.update_one(
        {"_id": current_user["_id"]},
        {"$set": update_fields}
    )

    return jsonify({"message": "Profil mis à jour avec succès", "updated_fields": update_fields}), 200


def init_auth(app):
    bcrypt.init_app(app)
    mongo.init_app(app)
    app.mongo = mongo



def init_user_index(user_id):
    user_folder = os.path.join("documents", str(user_id))
    index_folder = os.path.join(user_folder, "faiss_index")
    os.makedirs(index_folder, exist_ok=True)

    embedding = OllamaEmbeddings(model="nomic-embed-text")

    # Index FAISS vide avec un dummy document supprimé ensuite
    dummy_doc = [Document(page_content="Initialisation", metadata={"doc_id": "init"})]
    index = FAISS.from_documents(dummy_doc, embedding)
    index.index.reset()
    index.save_local(index_folder)

    print(f"[init_user_index] ✅ Index créé pour l'utilisateur {user_id}")