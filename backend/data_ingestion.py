# data_ingestion.py
import os
import requests
from flask import Flask, request, jsonify
import arxiv
from flask_socketio import SocketIO
from Bio import Entrez
import time
import xml.etree.ElementTree as ET

app = Flask(__name__)
socketio = SocketIO(app, cors_allowed_origins="*")  # Autoriser CORS si frontend sur un autre port
Entrez.email = "superfataou13@gmail.com"

UNPAYWALL_EMAIL = "superfataou13@gmail.com"  # Remplace par ton email validé sur Unpaywall

def search_arxiv(query="machine learning", max_results=3, save_dir="./data/arxiv"):
    os.makedirs(save_dir, exist_ok=True)
    search = arxiv.Search(query=query, max_results=max_results)
    results = []
    for result in search.results():
        paper_path = os.path.join(save_dir, f"{result.get_short_id()}.pdf")
        if not os.path.exists(paper_path):
            result.download_pdf(filename=paper_path)
        doc = {
            "title": result.title,
            "url": result.entry_id,
            "file_path": paper_path,
            "source": "arxiv",
            "summary": result.summary
        }
        socketio.emit('document', doc)
        results.append(doc)
    return results

def search_pubmed(query="machine learning", max_results=3, save_dir="./data/pubmed"):
    os.makedirs(save_dir, exist_ok=True)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    results = []

    for id_ in ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=id_, rettype="xml", retmode="text")
        article_data = fetch_handle.read()
        # 💡 Correction ici
        if isinstance(article_data, bytes):
            article_data = article_data.decode("utf-8", errors="replace")

        file_path = os.path.join(save_dir, f"{id_}.xml")
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(article_data)
            
        root = ET.fromstring(article_data)
        
        
        abstract = ""
        for abstract_text in root.findall(".//AbstractText"):
            if abstract_text.text:
                abstract += abstract_text.text.strip() + " "
                
                  # 🏷️ Titre
        title = ""
        title_elem = root.find(".//ArticleTitle")
        if title_elem is not None and title_elem.text:
            title = title_elem.text.strip()
        else:
            title = f"PubMed Article {id_}"
        doc = {
             "title": title,
            "file_path": file_path,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{id_}",
            "source": "pubmed",
            "summary": abstract.strip() if abstract else "Résumé non disponible"
        }

        socketio.emit('document', doc)
        results.append(doc)

    return results

# Fonction Semantic Scholar (commentée - à activer si besoin)
# def search_semantic_scholar(query="machine learning", max_results=2):
#     url = "https://api.semanticscholar.org/graph/v1/paper/search"
#     params = {
#         "query": query,
#         "limit": max_results,
#         "fields": "title,url,abstract,openAccessPdf"
#     }
#     results = []
#     retries = 3
#     backoff = 5  # secondes à attendre en cas de 429
# 
#     for attempt in range(retries):
#         try:
#             response = requests.get(url, params=params)
#             if response.status_code == 429:
#                 print(f"Rate limit hit, waiting {backoff}s before retry ({attempt+1}/{retries})")
#                 time.sleep(backoff)
#                 continue  # réessayer
#             response.raise_for_status()  # lance une exception pour autre erreur HTTP
#             data = response.json()
#             papers = data.get("data", [])
#             for paper in papers:
#                 open_access_pdf = paper.get("openAccessPdf")
#                 pdf_url = open_access_pdf.get("url") if open_access_pdf else None
#                 if pdf_url:
#                     doc = {
#                         "title": paper.get("title", "No title"),
#                         "url": paper.get("url", ""),
#                         "file_path": pdf_url,
#                         "summary": "",  # à compléter si besoin
#                         "source": "semantic_scholar"
#                     }
#                     socketio.emit('document', doc)
#                     results.append(doc)
#             break  # succès, sortir de la boucle retry
#         except requests.RequestException as e:
#             print(f"Erreur lors de la requête Semantic Scholar : {e}")
#             time.sleep(backoff)
#     else:
#         print("Échec après plusieurs tentatives.")
#     return results

def reconstruct_abstract(inverted_index):
    if not isinstance(inverted_index, dict):
        return "Résumé non disponible"

    word_positions = []
    for word, positions in inverted_index.items():
        for pos in positions:
            word_positions.append((pos, word))

    sorted_words = [word for _, word in sorted(word_positions)]
    return " ".join(sorted_words)


def search_openalex(query="machine learning", max_results=3):
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


