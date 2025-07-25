# routes/multi_functions.py
from flask import Blueprint, request, jsonify ,Response,send_file
import ollama
from utils import token_required
from data_ingestion import search_arxiv,search_openalex,search_pubmed,merge_and_deduplicate,search_core
from transformers import pipeline,AutoTokenizer ,MarianMTModel, MarianTokenizer
from utils.extract_text_from_pdf import extract_pdf_content ,summarize_text,generate_table_of_contents,extract_key_passages,chunk_text
from extensions import mongo
from datetime import datetime
from bson import ObjectId
import textwrap
import tempfile
import os
from fpdf import FPDF
import traceback
import requests
import json
from deep_translator import GoogleTranslator

bp_multi = Blueprint("multi_functions", __name__)

def ollama_query(prompt, model="llama3", stream=False):
    url = "http://localhost:11434/api/chat"
    headers = {"Content-Type": "application/json"}
    payload = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "stream": stream,
    }

    if not stream:
        response = requests.post(url, json=payload)
        response.raise_for_status()
        data = response.json()
        return data.get("message", {}).get("content", "") or data.get("text", "")

    # Mode streaming : retourne un générateur de texte
    def stream_response():
        with requests.post(url, json=payload, stream=True) as response:
            response.raise_for_status()
            for line in response.iter_lines():
                if line:
                    try:
                        clean = line.decode("utf-8").replace("data: ", "")
                        if clean.strip() == "[DONE]":
                            break
                        chunk = json.loads(clean)
                        yield chunk.get("message", {}).get("content", "")
                    except Exception as e:
                        print(f"[Streaming error] {e}")

    return stream_response()

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
            "user_id":user["_id"],
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


@bp_multi.route("/generate_report", methods=['POST', 'OPTIONS'])
@token_required
def generate_report(user):
    try:
        data = request.get_json()
        objective = data.get("objective", "")
        content_blocks = data.get("content_blocks", [])
        language = data.get("language", "fr")
        style = data.get("style", "scientifique")


        prompt = (
        f"Tu es un expert en rédaction {style}. En te basant uniquement sur l'objectif suivant :\n\n"
        f"{objective}\n\n"
        f"Génère un rapport complet et structuré en {language} comportant les parties appropriées "
        f"(comme l’introduction, la méthodologie, les résultats, la discussion, la conclusion...). "
        f"Utilise un style {style}, clair et professionnel, avec des titres et sous-titres adaptés."
    )


        def generate():
            for chunk in ollama_query(prompt, stream=True):  # 🧠 ton modèle doit supporter `stream=True`
                yield chunk

        return Response(generate(), mimetype="text/plain")

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@bp_multi.route("/generate_report_pdf", methods=["POST"])
@token_required
def generate_report_pdf(user):
    try:
        data = request.get_json()
        content = data.get("content", "")
        filename = data.get("filename", "rapport")

        pdf = FPDF()
        pdf.add_page()
        pdf.set_auto_page_break(auto=True, margin=15)
        pdf.set_font("Arial", size=12)

        for line in content.split("\n"):
            pdf.multi_cell(0, 10, line)

        filepath = f"/tmp/{filename}.pdf"
        pdf.output(filepath)

        return send_file(filepath, as_attachment=True, download_name=f"{filename}.pdf")

    except Exception as e:
        return jsonify({"error": str(e)}), 500

from langdetect import detect

@bp_multi.route('/paraphrase', methods=['POST'])
@token_required
def paraphrase(user):
    data = request.json
    texte = data.get('texte', '').strip()
    style = data.get('style', 'fluent').strip().lower()

    if not texte:
        return jsonify({'error': 'Pas de texte fourni'}), 400

    # Détecter la langue d'origine du texte
    try:
        detected_lang = detect(texte)
    except:
        detected_lang = "fr"  # fallback

    # Instructions personnalisées
    style_prompts = {
        "fluent": "Reformule ce texte de manière fluide et naturelle",
        "formal": "Reformule ce texte dans un style formel",
        "creative": "Reformule ce texte de manière créative et originale",
        "academic": "Reformule ce texte dans un style académique clair et rigoureux",
        "professional": "Reformule ce texte dans un style professionnel",
    }

    instruction = style_prompts.get(style, style_prompts["fluent"])

    # Ajoute la consigne explicite sur la langue
    prompt = (
        f"{instruction}.\n"
        f"Réponds dans la même langue que le texte d’origine ({detected_lang}).\n\n"
        f"Texte : {texte}"
    )

    try:
        reformulation = ollama_query(prompt)
        return jsonify({'paraphrase': reformulation})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


# --- Fonctions supplémentaires ---

@bp_multi.route("/ai_writer", methods=["POST"])
@token_required
def ai_writer(user):
    data = request.get_json()
    prompt = data.get("prompt", "").strip()
    if not prompt:
        return jsonify({"error": "Le prompt est requis"}), 400

    try:
        response = ollama.chat(
            model="llama3",
            messages=[{"role": "user", "content": prompt}],
        )
        answer = response.message.content if hasattr(response, "message") else "Pas de réponse"
        return jsonify({"text": answer})
    except Exception as e:
        return jsonify({"error": f"Erreur lors de l'appel à Ollama : {str(e)}"}), 500

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


THEMATIC_KEYWORDS = {
    "Biologie": ["biologie", "gène", "cellule", "ADN", "protéine", "organisme", "microbe", "bio"],
    "Mathématiques": ["math", "algèbre", "probabilité", "statistique", "géométrie", "équation", "matrice", "nombre"],
    "Informatique": ["IA", "intelligence artificielle", "machine learning", "deep learning", "réseau de neurones", "algorithme", "informatique", "computer"],
    "Physique": ["physique", "quantique", "relativité", "mécanique", "particule", "énergie", "champ"],
    "Chimie": ["chimie", "molécule", "réaction", "solvant", "atome", "composé"],
    # Tu peux enrichir avec d'autres thèmes
}

def classify_by_theme(doc):
    content = f"{doc.get('title', '')} {doc.get('summary', '')}".lower()
    for theme, keywords in THEMATIC_KEYWORDS.items():
        if any(kw in content for kw in keywords):
            return theme
    return "Autres"

def group_documents_by_theme(docs):
    grouped = {}
    for doc in docs:
        theme = classify_by_theme(doc)
        grouped.setdefault(theme, []).append(doc)
    return grouped


@bp_multi.route('/search', methods=['POST'])
@token_required
def handle_document_search(user):
    data = request.get_json()  # ✅ récupération correcte
    print("Received start_document_search:", data)
    language = data.get("language", "en")
    query = data.get("query", "").strip()
    max_results = int(data.get("max_results", 3))

    if not query:
        return jsonify({"error": "Aucune requête de recherche fournie"}), 400

    try:
        # Recherches
        arxiv_docs = search_arxiv(query, max_results, language)
        pubmed_docs = search_pubmed(query, max_results, language)
        openalex_docs = search_openalex(query, max_results, language)
        core_docs = search_core(query, max_results, language)
        merged_docs = merge_and_deduplicate(arxiv_docs, pubmed_docs, openalex_docs,core_docs)
        grouped = group_documents_by_theme(merged_docs)
        grouped_documents = [
            {"theme": theme, "documents": docs}
            for theme, docs in grouped.items()
        ]

        review_summary = generate_literature_summary(merged_docs, query , language)
        return jsonify({
            "grouped_documents": grouped_documents,
            "summary": review_summary
        }), 200

    except Exception as e:
        print("Erreur lors de la recherche :", str(e))
        return jsonify({"search_error": str(e)}), 500
    

def translate_text(text, target_lang):
    try:
        return GoogleTranslator(source='auto', target=target_lang).translate(text)
    except Exception:
        return text  # fallback

def generate_literature_summary(docs, query , language):
    from textwrap import shorten
    import requests

    # Préparer le prompt
    entries = ""
    for doc in docs[:5]:
        entries += f"Title: {doc.get('title')}\n"
        entries += f"Summary: {shorten(doc.get('summary', ''), width=500)}\n\n"

    prompt = (
        f"You are a scientific expert. Write a literature review about: '{query}'.\n\n"
        f"Here are a few articles:\n\n{entries}\n"
        "Write a structured and concise synthesis of the main ideas, approaches, trends, and limitations."
    )

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
        summary = result.get("response", "")

        # ✅ Traduction automatique en français
        summary_fr = translate_text(summary, language)
        return summary_fr

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
    

@bp_multi.route('/history', methods=['POST'])
@token_required
def save_user_search_history(user):
    data = request.get_json()
    query = data.get("query", "").strip()
    if not query:
        return jsonify({"error": "Requête vide"}), 400
    history = save_user_search_history(user["_id"], query)
    return jsonify({"history": history}), 200   