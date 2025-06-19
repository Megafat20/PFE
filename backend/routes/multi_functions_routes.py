# routes/multi_functions.py
from flask import Blueprint, request, jsonify
import ollama
from utils import token_required
from data_ingestion import search_arxiv,search_openalex,search_pubmed,merge_and_deduplicate
from transformers import pipeline,AutoTokenizer ,MarianMTModel, MarianTokenizer
from utils.extract_text_from_pdf import extract_pdf_content ,summarize_text,generate_table_of_contents,extract_key_passages,chunk_text
from extensions import mongo
from datetime import datetime
from bson import ObjectId
import textwrap
import tempfile
import os

import traceback

bp_multi = Blueprint("multi_functions", __name__)

def ollama_query(prompt):
    try:
        response = ollama.chat(model="llama3", messages=[{"role": "user", "content": prompt}])
        print("OLLAMA RESPONSE:", response)  # debug

        # Vérification robuste du format de la réponse
        if "text" in response:
            return response["text"]
        elif "message" in response and "content" in response["message"]:
            return response["message"]["content"]
        else:
            raise ValueError(f"Réponse inattendue de Ollama : {response}")
    except Exception as e:
        # Loguer ou gérer plus proprement selon besoin
        raise RuntimeError(f"Erreur lors de l'appel à Ollama : {str(e)}")

def create_prompt(template: str, **kwargs) -> str:
    return template.format(**kwargs)





summarizer = pipeline("summarization", model="facebook/bart-large-cnn")
tokenizer = AutoTokenizer.from_pretrained("facebook/bart-large-cnn")


@bp_multi.route('/upload_summarize', methods=['POST', 'OPTIONS'])
@token_required
def upload_summarize(user):
    if request.method == 'OPTIONS':
        return '', 200

    try:
        if 'file' not in request.files:
            return jsonify({"error": "Aucun fichier fourni"}), 400

        file = request.files['file']
        if not file.filename or not file.filename.endswith('.pdf'):
            return jsonify({"error": "Fichier PDF invalide"}), 400

        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
            file.save(tmp.name)
            tmp_path = tmp.name

        extracted = extract_pdf_content(tmp_path)
        os.remove(tmp_path)

        text = extracted.get("text_sections", "")
        if not isinstance(text, str) or not text.strip():
            return jsonify({"error": "Le document est vide ou illisible"}), 400

        summary = summarize_text(text)
        toc = generate_table_of_contents(text)
        highlights = extract_key_passages(text)

        mongo.db.summary_documents.insert_one({
            "user_id": ObjectId(user["id"]),
            "filename": file.filename,
            "content": text,
            "summary": summary,
            "toc": toc,
            "highlights": highlights
        })

        return jsonify({
            "summary": summary,
            "toc": toc,
            "highlights": highlights
        })

    except Exception as e:
        return jsonify({"error": f"Erreur lors de l'analyse : {str(e)}"}), 500




@bp_multi.route('/summary', methods=['POST'])
@token_required
def summary(user):
    try:
        data = request.get_json()
        if data is None:
            return jsonify({"error": "Le corps de la requête n'est pas un JSON valide"}), 400

        text = data.get("documents", "")
        if not isinstance(text, str):
            return jsonify({"error": "Le champ 'documents' doit être une chaîne de caractères"}), 400
        if not text.strip():
            return jsonify({"error": "Aucun texte fourni"}), 400

        # Nettoyage et découpage
        text = " ".join(text.split())
        chunks = chunk_text(text, max_tokens=1024)

        summaries = []
        for chunk in chunks:
            try:
                chunk_tokens = tokenizer.encode(chunk, add_special_tokens=True)
                max_input_len = 1024 - 150
                if len(chunk_tokens) > max_input_len:
                    chunk_tokens = chunk_tokens[:max_input_len]
                    chunk = tokenizer.decode(chunk_tokens, skip_special_tokens=True)

                max_len = min(150, max(30, len(chunk_tokens) // 2))
                summary = summarizer(chunk, max_length=max_len, min_length=30, do_sample=False)
                summaries.append(summary[0]['summary_text'])
            except Exception as e:
                print("Erreur chunk résumé:", e)
                continue

        final_summary = " ".join(summaries)
        toc = generate_table_of_contents(text)
        highlights = extract_key_passages(text)

        return jsonify({
            "summary": final_summary,
            "toc": toc,
            "highlights": highlights
        })

    except Exception as e:
        print("❌ Erreur backend résumé:", e)
        return jsonify({"error": "Erreur serveur : " + str(e)}), 500
    
def chunk_text(text, max_length=1024):
    """Diviser le texte en morceaux gérables pour le modèle."""
    return textwrap.wrap(text, max_length)

@bp_multi.route("/summary/<doc_id>", methods=["GET"])
def summarize_document(doc_id):
    try:
        document = mongo.db.documentsExterne.find_one({"_id": ObjectId(doc_id)})

        if not document:
            return jsonify({"error": "Document non trouvé"}), 404

        content = document.get("content", "")
        if not content.strip():
            return jsonify({"error": "Document vide"}), 400

        chunks = chunk_text(content)
        summaries = []

        for chunk in chunks:
            summary = summarizer(chunk, max_length=150, min_length=40, do_sample=False)[0]['summary_text']
            summaries.append(summary)

        full_summary = " ".join(summaries)

        # Mise à jour du document avec le résumé
        mongo.db.documentsExterne.update_one(
            {"_id": ObjectId(doc_id)},
            {"$set": {"summary": full_summary}}
        )

        return jsonify({
            "summary": full_summary,
            "chunks": len(chunks)
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500
        
@bp_multi.route('/multi/summarize', methods=['POST'])
@token_required
def summarize():
    data = request.json
    text = data.get('text', '')

    if not text:
        return jsonify({"error": "Aucun texte fourni"}), 400

    try:
        # HuggingFace summarization
        summary_list = summarizer(text, max_length=150, min_length=40, do_sample=False)
        summary_text = summary_list[0]['summary_text']

        return jsonify({"summary": summary_text})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@bp_multi.route('/translate', methods=['POST'])
@token_required
def translate_text(user):
    data = request.json
    text = data.get("text", "")
    target_language = data.get("language", "")

    if not text or not target_language:
        return jsonify({"error": "Texte ou langue cible manquant"}), 400

    prompt = f"Traduis le texte suivant en {target_language} (langue ) :\n\n{text}"
    try:
        translation = ollama_query(prompt)
        return jsonify({"translation": translation})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

import requests




@bp_multi.route('/translate_texte', methods=['POST'])
def translate_texte():
    data = request.json
    if not data or 'text' not in data:
        return jsonify({'error': 'Missing "text" field in JSON body'}), 400

    texts = data['text']
    target_lang = data.get('language', 'fr')
    source_lang = data.get('source_language', 'fr')  # ici par défaut 'fr' si tu veux

    if isinstance(texts, str):
        texts = [texts]

    # Si source == target, on ne traduit pas
    if source_lang == target_lang:
        return jsonify({'translations': texts})

    translated_texts = []
    for text in texts:
        try:
            response = requests.post(
                url="http://localhost:5001/translate",
                json={
                    "q": text,
                    "source": source_lang,
                    "target": target_lang,
                    "format": "text"
                },
                headers={"Content-Type": "application/json"}
            )
            if response.status_code == 200:
                translated = response.json().get("translatedText")
                translated_texts.append(translated)
            else:
                print("Erreur LibreTranslate:", response.text)
                translated_texts.append(text)  # fallback
        except Exception as e:
            print(f"Exception pendant la traduction : {str(e)}")
            translated_texts.append(text)

    return jsonify({'translations': translated_texts})


@bp_multi.route('/generate', methods=['POST'])
@token_required
def generate_intro(user):
    data = request.json
    topic = data.get("topic", "")
    if not topic:
        return jsonify({"error": "Pas de sujet fourni"}), 400

    template = "Écris une introduction académique sur le sujet : {topic}"
    prompt = create_prompt(template, topic=topic)
    try:
        result = ollama_query(prompt)
        return jsonify({"introduction": result})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp_multi.route('/rapport', methods=['POST'])
@token_required
def rapport(user):
    data = request.json
    info = data.get('info', '')
    if not info:
        return jsonify({'error': 'Pas d\'information fournie'}), 400
    prompt = f"Rédige un rapport détaillé à partir des informations suivantes : {info}"
    try:
        rapport = ollama_query(prompt)
        return jsonify({'rapport': rapport})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp_multi.route('/paraphrase', methods=['POST'])
@token_required
def paraphrase(user):
    data = request.json
    texte = data.get('texte', '').strip()
    style = data.get('style', 'fluent').strip().lower()

    if not texte:
        return jsonify({'error': 'Pas de texte fourni'}), 400

    # Dictionnaire de styles → instructions personnalisées
    style_prompts = {
        "fluent": "Reformule ce texte de manière fluide et naturelle :",
        "formal": "Reformule ce texte dans un style formel :",
        "creative": "Reformule ce texte de manière créative et originale :",
        "academic": "Reformule ce texte dans un style académique clair et rigoureux :",
        "professional": "Reformule ce texte dans un style professionnel :",
    }

    instruction = style_prompts.get(style, style_prompts["fluent"])
    prompt = f"{instruction} {texte}"

    try:
        reformulation = ollama_query(prompt)
        return jsonify({'paraphrase': reformulation})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# --- Fonctions supplémentaires ---

@bp_multi.route('/extract_key_info', methods=['POST'])
@token_required
def extract_key_info(user):
    data = request.json
    query = data.get('query', '')
    if not query:
        return jsonify({'error': 'Pas de requête fournie'}), 400

    prompt = f"""
    Tu es un assistant IA scientifique. À partir de la requête ci-dessous,
    extrait les citations, références et informations clés.

    Requête : {query}
    """
    try:
        result = ollama_query(prompt)
        return jsonify({'key_info': result})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp_multi.route('/correct_grammar', methods=['POST'])
@token_required
def correct_grammar(user):
    data = request.json
    text = data.get('text', '')
    if not text:
        return jsonify({'error': 'Pas de texte fourni'}), 400

    prompt = f"Améliore la grammaire et le style scientifique du passage suivant :\n{text}"
    try:
        corrected = ollama_query(prompt)
        return jsonify({'corrected_text': corrected})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp_multi.route('/analyze_trends', methods=['POST'])
@token_required
def analyze_trends(user):
    data = request.json
    topic = data.get('topic', '')
    if not topic:
        return jsonify({'error': 'Pas de sujet fourni'}), 400

    prompt = f"Analyse les tendances récentes dans les documents scientifiques et brevets concernant : {topic}"
    try:
        analysis = ollama_query(prompt)
        return jsonify({'trend_analysis': analysis})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
@bp_multi.route('/search', methods=['POST'])
@token_required
def handle_document_search(user):
    data = request.get_json()  # ✅ récupération correcte
    print("Received start_document_search:", data)

    query = data.get("query", "").strip()
    max_results = int(data.get("max_results", 3))

    if not query:
        return jsonify({"error": "Aucune requête de recherche fournie"}), 400

    try:
        # Recherches
        arxiv_docs = search_arxiv(query, max_results)
        pubmed_docs = search_pubmed(query, max_results)
        openalex_docs = search_openalex(query, max_results)

        # Fusion et suppression des doublons
        merged_docs = merge_and_deduplicate(arxiv_docs, pubmed_docs, openalex_docs)
        review_summary = generate_literature_summary(merged_docs, query)
        return jsonify({
                "documents": merged_docs,
                "summary": review_summary  # 👈 ajouter la synthèse ici
            }), 200

    except Exception as e:
        print("Erreur lors de la recherche :", str(e))
        return jsonify({"search_error": str(e)}), 500
    
def generate_literature_summary(docs, query):
    from textwrap import shorten
    import requests

    # Préparer le prompt
    entries = ""
    for doc in docs[:5]:
        entries += f"Title: {doc.get('title')}\n"
        entries += f"Summary: {shorten(doc.get('summary', ''), width=500)}\n\n"

    prompt = (
        f"Tu es un expert scientifique. Fais une revue de littérature sur : '{query}'.\n\n"
        f"Voici quelques articles trouvés :\n\n{entries}\n"
        "Élabore une synthèse des idées principales, des approches, des tendances et des limites."
    )

    # Appel à Ollama (DeepSeek ou autre)
    try:
        response = requests.post(
            "http://localhost:11434/api/generate",
            json={
                "model": "deepseek-coder:latest",  # ou deepseek-llm
                "prompt": prompt,
                "stream": False
            }
        )
        result = response.json()
        return result.get("response", "")
    except Exception as e:
        print("Erreur LLM:", e)
        return "Synthèse non disponible (erreur lors de l'appel au modèle)."
    
def save_user_search_history(user_id, query):
    mongo.db.search_history.update_one(
        {"user_id": user_id},
        {
            "$push": {
                "history": {
                    "$each": [{"query": query, "timestamp": datetime.utcnow()}],
                    "$slice": -10  # garder les 10 dernières recherches
                }
            }
        },
        upsert=True
    )
def get_user_search_history(user_id):
    doc = mongo.db.search_history.find_one({"user_id": user_id})
    if doc and "history" in doc:
        return sorted(doc["history"], key=lambda x: x["timestamp"], reverse=True)
    return []
@bp_multi.route('/history', methods=['GET'])
@token_required
def get_search_history(user):
    history = get_user_search_history(user["id"])
    return jsonify({"history": history}), 200