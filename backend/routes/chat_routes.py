from flask import Blueprint, request, jsonify
from bson import ObjectId
from datetime import datetime
from extensions import mongo
from utils import token_required
from langchain.vectorstores import FAISS
from langchain.embeddings import OllamaEmbeddings
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
import ollama
import os


chat_bp = Blueprint('chat', __name__)

EXTERNAL_UPLOAD_FOLDER = "./external_uploads"

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


@chat_bp.route("/chat_document", methods=["POST"])
@token_required
def chat_document(user):
    data = request.get_json()
    document_id = data.get("document_id")
    question = data.get("question", "").strip()

    if not document_id or not question:
        return jsonify({"error": "document_id et question requis"}), 400

    user_folder = os.path.join(EXTERNAL_UPLOAD_FOLDER, str(user["_id"]))
    index_path = os.path.join(user_folder, "faiss_index")

    if not os.path.exists(index_path):
        return jsonify({"error": "Index FAISS non trouvé pour cet utilisateur"}), 404

    embedding = OllamaEmbeddings(model="nomic-embed-text")

    try:
        store = FAISS.load_local(index_path, embeddings=embedding, allow_dangerous_deserialization=True)
    except Exception as e:
        return jsonify({"error": f"Impossible de charger l'index FAISS: {str(e)}"}), 500

    docs = store.similarity_search(question, k=5)

    context = "\n\n".join([doc.page_content for doc in docs])
    max_context_length = 4000
    if len(context) > max_context_length:
        context = context[:max_context_length] + "..."

    prompt = f"""Tu es un assistant IA qui répond uniquement à partir du contexte fourni ci-dessous.
Contexte:
{context}

Question:
{question}

Réponds précisément, en français si la question est en français."""
   
    try:
        response = ollama.chat(
            model="llama3",
            messages=[{"role": "user", "content": prompt}]
        )
        answer = response.get("choices", [{}])[0].get("message", {}).get("content", "Pas de réponse")
    except Exception as e:
        return jsonify({"error": f"Erreur lors de l'appel à Ollama : {str(e)}"}), 500

    sources = []
    for doc in docs:
        sources.append({
            "content": doc.page_content[:500],
            "metadata": doc.metadata,
        })

    return jsonify({"answer": answer, "sources": sources})
