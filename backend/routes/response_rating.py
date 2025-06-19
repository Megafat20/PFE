from flask import Blueprint, request, jsonify

from utils import token_required
from datetime import datetime
from extensions import mongo
from utils.faiss_index import get_cached_embedding
import numpy as np
from bson import ObjectId



bp_rating = Blueprint('rating', __name__)


@bp_rating.route("/validate_answer", methods=["POST"])
@token_required
def validate_answer(user):
    try:
        data = request.json
        conv_id = data.get("conversation_id")
        question = data.get("question")
        answer = data.get("answer")
        rating = data.get("rating")
        validated = data.get("validated", False)

        if not conv_id or rating is None or answer is None or question is None:
            return jsonify({"error": "Données manquantes"}), 400

        try:
            conv_id_obj = ObjectId(conv_id)
        except Exception:
            return jsonify({"error": "ID de conversation invalide"}), 400

        # 🔍 Vérifie que la conversation existe bien
        conv = mongo.db.conversations.find_one({"_id": conv_id_obj, "user_id": user["_id"]})
        if not conv:
            return jsonify({"error": "Conversation introuvable"}), 404

        # ✅ Mise à jour du message dans la conversation
        messages = conv.get("messages", [])
        updated = False
        for i, msg in enumerate(messages):
            if msg.get("content") == answer or msg.get("text") == answer:
                messages[i]["rating"] = rating
                messages[i]["validated"] = validated
                messages[i]["question"] = question
                updated = True
                break

        if not updated:
            return jsonify({"error": "Message introuvable"}), 404

        mongo.db.conversations.update_one(
            {"_id": conv_id_obj, "user_id": user["_id"]},
            {"$set": {"messages": messages}}
        )

        # ✅ Si note suffisante, stocker dans validated_answers avec vectorisation via cache
        if validated and rating >= 4:
            question_vector = get_cached_embedding(
                text=question,
                uid=user["_id"],
                embedding_type="query"
            )

            mongo.db.validated_answers.insert_one({
                "user_id": user["_id"],
                "conversation_id": conv_id_obj,
                "question": question,
                "embedding": question_vector,
                "answer": answer,
                "rating": rating,
                "validated": True,
                "timestamp": datetime.utcnow()
            })

        return jsonify({"msg": "Note enregistrée et réponse validée."})

    except Exception as e:
        import traceback
        traceback_str = traceback.format_exc()
        print(traceback_str)
        return jsonify({"error": "Erreur serveur interne", "details": str(e)}), 500

