# data_ingestion.py
import os
import requests
from flask import Flask, request, jsonify
import arxiv
from flask_socketio import SocketIO
from Bio import Entrez
import time
import xml.etree.ElementTree as ET
from deep_translator import GoogleTranslator

def translate_text(text, target_lang):
    try:
        return GoogleTranslator(source='auto', target=target_lang).translate(text)
    except Exception:
        return text
app = Flask(__name__)
socketio = SocketIO(app, cors_allowed_origins="*")  # Autoriser CORS si frontend sur un autre port
Entrez.email = "superfataou13@gmail.com"

UNPAYWALL_EMAIL = "superfataou13@gmail.com"  # Remplace par ton email validé sur Unpaywall

def search_arxiv(query="machine learning", max_results=3 , language="en"):
    search = arxiv.Search(query=query, max_results=max_results)
    results = []
    for result in search.results():
        title = result.title
        summary = result.summary
        if language != "en":
            title = translate_text(title, language)
            summary = translate_text(summary, language)
        doc = {
            "title": title,
            "url": result.entry_id,
            "source": "arxiv",
            "summary": summary
        }
        socketio.emit('document', doc)
        results.append(doc)
    return results

def search_pubmed(query="machine learning", max_results=3, language="en"):
    try:
        if language != "en":
            query = translate_text(query, "en")
    except Exception as e:
        print(f"[⚠️ Traduction échouée] Query non traduite : {e}")

    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    results = []

    for id_ in ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=id_, rettype="xml", retmode="text")
        article_data = fetch_handle.read()

        if isinstance(article_data, bytes):
            article_data = article_data.decode("utf-8", errors="replace")
        root = ET.fromstring(article_data)

        # Résumé
        abstract = ""
        for abstract_text in root.findall(".//AbstractText"):
            if abstract_text.text:
                abstract += abstract_text.text.strip() + " "

        # Titre
        title = ""
        title_elem = root.find(".//ArticleTitle")
        if title_elem is not None and title_elem.text:
            title = title_elem.text.strip()
        else:
            title = f"PubMed Article {id_}"

        # Traduction si besoin
        if language != "en":
            title = translate_text(title, language)
            if abstract:
                abstract = translate_text(abstract, language)

        doc = {
            "title": title,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{id_}",
            "source": "pubmed",
            "summary": abstract.strip() if abstract else "Résumé non disponible"
        }

        socketio.emit('document', doc)
        results.append(doc)

    return results

def reconstruct_abstract(inverted_index):
    if not isinstance(inverted_index, dict):
        return "Résumé non disponible"

    word_positions = []
    for word, positions in inverted_index.items():
        for pos in positions:
            word_positions.append((pos, word))

    sorted_words = [word for _, word in sorted(word_positions)]
    return " ".join(sorted_words)


def search_openalex(query="machine learning", max_results=3, language="en"):
    try:
        if language != "en":
            query = translate_text(query, "en")
    except Exception as e:
        print(f"[⚠️ Traduction échouée] Query non traduite : {e}")

    url = f"https://api.openalex.org/works?search={query}&per-page={max_results}"
    response = requests.get(url)
    results = []

    if response.status_code != 200:
        print("Erreur OpenAlex:", response.text)
        return results

    data = response.json()
    for item in data.get('results', []):
        title = item.get("title", "Titre non disponible")
        url = item.get("id", "#")
        authors = [auth["author"]["display_name"] for auth in item.get("authorships", [])]
        abstract_data = item.get("abstract_inverted_index")

        summary = reconstruct_abstract(abstract_data) if abstract_data else "Résumé non disponible"
        primary_location = item.get("primary_location") or {}
        source = primary_location.get("source") or {}
        pdf_url = source.get("url") or url

        # Traduction
        if language != "en":
            try:
                title = translate_text(title,language)
                if summary and summary != "Résumé non disponible":
                    summary = translate_text(summary, language)
            except Exception as e:
                print(f"[⚠️ Traduction titre/résumé OpenAlex échouée] : {e}")

        doc = {
            "title": title,
            "authors": authors,
            "summary": summary,
            "pdf_url": pdf_url,
            "url": url,
            "source": "openalex"
        }

        socketio.emit('document', doc)
        results.append(doc)

    return results

def search_core(query="machine learning", max_results=3, language="en"):
    results = []
    base_url = "https://api.core.ac.uk/v3/search/works"
    try:
        if language != "en":
            query = translate_text(query, "en")
    except Exception as e:
        print(f"[⚠️ Traduction échouée] Query non traduite : {e}")

    params = {
        "q": query,
        "page": 1,
        "pageSize": max_results,
        "fields": "title,authors,abstract,openAccessUrl,url",
    }
    headers = {
        "Accept": "application/json",
        # "Authorization": "Bearer TON_API_KEY",
    }

    try:
        response = requests.get(base_url, params=params, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()

        for item in data.get("results", []):
            title = item.get("title") or "Titre non disponible"
            authors = [a.get("name") for a in item.get("authors", []) if a.get("name")] or []
            summary = item.get("abstract") or "Résumé non disponible"
            pdf_url = item.get("openAccessUrl") or item.get("url") or "#"

            # Traduction
            if language != "en":
                try:
                    title =translate_text(title,language).text
                    if summary and summary != "Résumé non disponible":
                        summary = translate_text(summary,language).text
                except Exception as e:
                    print(f"[⚠️ Traduction titre/résumé CORE échouée] : {e}")

            doc = {
                "title": title,
                "authors": authors,
                "summary": summary,
                "pdf_url": pdf_url,
                "url": item.get("url") or "#",
                "source": "core",
            }

            socketio.emit('document', doc)
            results.append(doc)

    except Exception as e:
        print(f"Erreur CORE API: {e}")

    return results




def merge_and_deduplicate(*doc_sources):
    seen_urls = set()
    seen_titles = set()
    merged = []

    for docs in doc_sources:
        for doc in docs:
            url = doc.get('url', '').lower()
            title = doc.get('title', '').lower()
            if url in seen_urls or title in seen_titles:
                continue
            seen_urls.add(url)
            seen_titles.add(title)
            merged.append(doc)

    merged.sort(key=lambda d: d.get('score', 0), reverse=True)
    return merged


