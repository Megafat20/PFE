from flask import Blueprint, request, jsonify
from bson import ObjectId
from datetime import datetime
from extensions import mongo
from utils import token_required

chat_bp = Blueprint('chat', __name__)

@chat_bp.route("/save_chat", methods=["POST"])
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

    for msg in messages:
        mongo.db.conversations.update_one(
            {"_id": ObjectId(conv_id)},
            {"$push": {"messages": {
                "role": msg["role"],
                "content": msg["content"],
                "timestamp": msg.get("timestamp", datetime.utcnow())
            }}}
        )

    if conv["title"] == "Nouvelle conversation" and len(conv.get("messages", [])) == 0:
        first_user_msg = next((m for m in messages if m["role"] == "user"), None)
        if first_user_msg:
            generated_title = first_user_msg["content"].strip().split("?")[0][:50]
            mongo.db.conversations.update_one(
                {"_id": ObjectId(conv_id)},
                {"$set": {"title": generated_title}}
            )

    return jsonify({"status": "ok"})