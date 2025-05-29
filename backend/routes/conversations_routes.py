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
    return jsonify([{"id": str(c["_id"]), "title": c["title"], "created_at": c["created_at"]} for c in convs])

@bp_conversations.route("/conversations/<cid>", methods=["GET"])
@token_required
def get_conv(user, cid):
    conv = mongo.db.conversations.find_one({"_id": ObjectId(cid), "user_id": user["_id"]})
    if not conv:
        return jsonify({"error": "Introuvable"}), 404
    return jsonify({"id": str(conv["_id"]), "title": conv["title"], "messages": conv["messages"]})

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
