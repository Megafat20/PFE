from flask import Blueprint, request, jsonify, current_app
from extensions import mongo
from flask_bcrypt import Bcrypt
import jwt
from datetime import datetime, timedelta
from functools import wraps
from bson.objectid import ObjectId

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
    mongo.db.users.insert_one({
        "name": name,
        "email": email,
        "password": hashed_pw
    })
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


def init_auth(app):
    bcrypt.init_app(app)
    mongo.init_app(app)
    app.mongo = mongo
