from flask import Blueprint, jsonify
from utils.auth import token_required

bp_protected = Blueprint("protected", __name__)

@bp_protected.route("/me", methods=["GET"])
@token_required
def me(user):
    return jsonify({
        "id": str(user["_id"]),
        "email": user["email"],
        "name": user.get("name")
    })
