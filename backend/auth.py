from flask import Blueprint, request, jsonify, current_app
from flask_pymongo import PyMongo
from flask_bcrypt import Bcrypt
import jwt
from datetime import datetime, timedelta
from functools import wraps
from bson.objectid import ObjectId

auth_bp = Blueprint('auth', __name__)
bcrypt = Bcrypt()

def init_auth(app):
    bcrypt.init_app(app)
    mongo = PyMongo(app)
    app.mongo = mongo
    app.register_blueprint(auth_bp, url_prefix="/auth")

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
    mongo.db.users.insert_one({"name": name, "email": email, "password": hashed_pw})
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
    def decorated(*args, **kwargs):
        mongo = current_app.mongo
        token = request.headers.get("Authorization", "").replace("Bearer ", "")
        if not token:
            return jsonify({"error": "Token requis"}), 403
        try:
            data = jwt.decode(token, current_app.config["SECRET_KEY"], algorithms=["HS256"])
            user = mongo.db.users.find_one({"_id": ObjectId(data["user_id"])})
            if not user:
                return jsonify({"error": "Utilisateur non trouvé"}), 403
        except Exception as e:
            return jsonify({"error": f"Token invalide : {e}"}), 403
        return f(user, *args, **kwargs)
    return decorated
