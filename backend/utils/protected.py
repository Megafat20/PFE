from flask import Blueprint, jsonify
from .auth import token_required

protected_bp = Blueprint("protected", __name__)

@protected_bp.route("/me", methods=["GET"])
@token_required
def me(user):
    return jsonify({
        "id": str(user["_id"]),
        "email": user["email"],
        "name": user.get("name")
    })
