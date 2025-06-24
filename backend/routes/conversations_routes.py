from flask import Blueprint, request, jsonify
from bson import ObjectId
from datetime import datetime
from extensions import mongo
from utils import token_required

bp_conversations = Blueprint('conversations', __name__)

@bp_conversations.route("/conversations", methods=["POST"])
@token_required
def create_conv(user):
    title = (request.json or {}).get("title", "Nouvelle conversation")
    res = mongo.db.conversations.insert_one({
        "user_id": user["_id"],
        "title": title,
        "messages": [],
        "created_at": datetime.utcnow()
    })
    return jsonify({"conversation_id": str(res.inserted_id)}), 201

@bp_conversations.route("/conversations", methods=["GET"])
@token_required
def list_conv(user):
    convs = mongo.db.conversations.find({"user_id": user["_id"]})
    result = []
    for c in convs:
        created_at = c.get("created_at")
        created_at_str = created_at.isoformat() if created_at else None
        result.append({
            "id": str(c["_id"]),
            "title": c["title"],
            "created_at": created_at_str
        })
    return jsonify(result)



@bp_conversations.route("/conversations/search", methods=["GET"])
@token_required
def search_conversation_messages(user):
    query = request.args.get("q", "").strip().lower()
    if not query:
        return jsonify({"error": "Paramètre de recherche manquant"}), 400

    results = []
    conversations = mongo.db.conversations.find({"user_id": user["_id"]})
    
    for conv in conversations:
        conv_id = str(conv["_id"])
        title = conv["title"]
        for msg in conv.get("messages", []):
            if query in msg.get("content", "").lower():
                results.append({
                    "conversation_id": conv_id,
                    "conversation_title": title,
                    "role": msg.get("role", ""),
                    "content": msg.get("content", ""),
                    "timestamp": msg.get("timestamp")
                })

    return jsonify(results)

@bp_conversations.route("/conversations/<cid>", methods=["GET"])
@token_required
def get_conv(user, cid):
    conv = mongo.db.conversations.find_one({"_id": ObjectId(cid), "user_id": user["_id"]})
    if not conv:
        return jsonify({"error": "Introuvable"}), 404

    # Optionnel : nettoyer les messages si besoin (par exemple supprimer _id des messages)
    messages = conv.get("messages", [])
    for m in messages:
        if "_id" in m:
            del m["_id"]

    return jsonify({
        "id": str(conv["_id"]),
        "title": conv["title"],
        "messages": messages
    })

@bp_conversations.route("/conversations/<cid>", methods=["PUT"])
@token_required
def rename_conv(user, cid):
    new_title = (request.json or {}).get("title")
    if not new_title:
        return jsonify({"error": "title manquant"}), 400
    res = mongo.db.conversations.update_one(
        {"_id": ObjectId(cid), "user_id": user["_id"]},
        {"$set": {"title": new_title}}
    )
    if res.matched_count == 0:
        return jsonify({"error": "Introuvable"}), 404
    return jsonify({"msg": "Title updated"})

@bp_conversations.route("/conversations/<cid>", methods=["DELETE"])
@token_required
def delete_conv(user, cid):
    res = mongo.db.conversations.delete_one({"_id": ObjectId(cid), "user_id": user["_id"]})
    if res.deleted_count == 0:
        return jsonify({"error": "Introuvable"}), 404
    return jsonify({"msg": "Conversation supprimée"})


def get_conversation_history(conversation_id, k=5):
    
    try:
        conv_oid = ObjectId(conversation_id)
    except Exception:
        conv_oid = conversation_id 
    """
    Récupère les k derniers échanges depuis la conversation Mongo.
    """
    convo = mongo.db.conversations.find_one({"_id": conv_oid})
    if not convo or "messages" not in convo:
        return []

    all_messages = convo["messages"]
    # Prendre les k derniers échanges = 2k derniers messages
    history = all_messages[-(k * 2):]
    return history

