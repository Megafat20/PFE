from extensions import mongo
from utils import token_required
from flask import Blueprint, jsonify,request
from bson import ObjectId

bp_notifications = Blueprint("notifications", __name__)


@bp_notifications.route("/enable_notifications", methods=["POST"])
@token_required
def enable_notifications(user):
    data = request.get_json()
    enabled = data.get("enabled", False)
    interests = data.get("interests", [])

    if not isinstance(enabled, bool):
        return jsonify({"error": "enabled doit être un booléen"}), 400
    if not isinstance(interests, list):
        return jsonify({"error": "interests doit être une liste"}), 400

    update_data = {
        "notifications_enabled": enabled,
        "interests": interests if enabled else [],
    }
    mongo.db.users_profile.update_one({"_id": user["_id"]}, {"$set": update_data})
    return jsonify({"success": True})

@bp_notifications.route("/notifications", methods=["GET"])
@token_required
def get_notifications(user):
    notifications = list(mongo.db.notifications.find({"user_id": user["_id"], "read": False}))
    for n in notifications:
        n["_id"] = str(n["_id"])
        n["user_id"] = str(n["user_id"])
        n["created_at"] = n["created_at"].isoformat()
    return jsonify(notifications)

@bp_notifications.route("/notifications/<notif_id>/read", methods=["POST"])
@token_required
def mark_notification_read(user, notif_id):
    notif = mongo.db.notifications.find_one({"_id": ObjectId(notif_id), "user_id": user["_id"]})
    if not notif:
        return jsonify({"error": "Notification introuvable"}), 404
    mongo.db.notifications.update_one({"_id": ObjectId(notif_id)}, {"$set": {"read": True}})
    return jsonify({"success": True})