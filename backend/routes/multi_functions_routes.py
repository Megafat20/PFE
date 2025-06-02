# routes/multi_functions.py
from flask import Blueprint, request, jsonify
import ollama
from utils import token_required
from data_ingestion import search_arxiv,search_openalex,search_pubmed,merge_and_deduplicate
from transformers import pipeline

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

@bp_multi.route('/summary', methods=['POST'])
@token_required
def summary(user):
    data = request.json
    documents = data.get("documents", "")
    if not documents:
        return jsonify({"error": "Aucun document fourni"}), 400

    prompt = f"Résume le texte suivant de façon claire et concise:\n{documents}"
    try:
        result = ollama_query(prompt)
        return jsonify({"summary": result})
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


# @bp_multi.route('/translate_bulk', methods=['POST'])
# @token_required
# def translate_bulk():
#     data = request.get_json()
#     texts = data.get("texts", [])
#     language = data.get("language", "fr")

#     if not texts or not isinstance(texts, list):
#         return jsonify({"error": "Liste de textes invalide"}), 400

#     try:
#         translations = translator.translate(texts, target_language=language)
#         translated_texts = [item["translatedText"] for item in translations]
#         return jsonify({"translations": translated_texts})
#     except Exception as e:
#         return jsonify({"error": str(e)}), 500


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
    texte = data.get('texte', '')
    if not texte:
        return jsonify({'error': 'Pas de texte fourni'}), 400
    prompt = f"Reformule ce texte avec d'autres mots tout en gardant le sens : {texte}"
    try:
        reformulation = ollama_query(prompt)
        return jsonify({'paraphrase': reformulation})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


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