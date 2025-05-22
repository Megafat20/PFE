# data_ingestion.py
import os
import requests
from flask import Flask, request, jsonify
import arxiv
from flask_socketio import SocketIO
from Bio import Entrez

app = Flask(__name__)
socketio = SocketIO(app)
Entrez.email = "superfataou13@gmail.com"

UNPAYWALL_EMAIL = "superfataou13@gmail.com"  # Remplace par ton email validé sur Unpaywall

def search_arxiv(query="machine learning", max_results=5, save_dir="./data/arxiv"):
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

def search_pubmed(query="machine learning", max_results=5, save_dir="./data/pubmed"):
    os.makedirs(save_dir, exist_ok=True)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    results = []
    for id_ in ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=id_, rettype="xml", retmode="text")
        article_data = fetch_handle.read()
        file_path = os.path.join(save_dir, f"{id_}.xml")
        with open(file_path, "wb") as f:  # <-- mode binaire
            f.write(article_data)
        doc = {
            "title": f"PubMed Article {id_}",
            "file_path": file_path,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{id_}",
            "source": "pubmed",
            "summary": "Résumé non disponible pour PubMed."
        }
        socketio.emit('document', doc)
        results.append(doc)
    return results


def search_semantic_scholar(query="machine learning", max_results=5):
    url = f"https://api.semanticscholar.org/graph/v1/paper/search"
    params = {
        "query": query,
        "limit": max_results,
        "fields": "title,url,abstract,openAccessPdf"
    }
    results = []
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()  # déclenche une exception pour codes erreur HTTP
        data = response.json()
        papers = data.get("data", [])
        for paper in papers:
            open_access_pdf = paper.get("openAccessPdf")
            pdf_url = open_access_pdf.get("url") if open_access_pdf else None
            if pdf_url:
                doc = {
                    "title": paper.get("title", "No title"),
                    "url": paper.get("url", ""),
                    "file_path": pdf_url,
                    "summary": paper.get("abstract", ""),
                    "source": "semantic_scholar"
                }
                socketio.emit('document', doc)
                results.append(doc)
    except requests.RequestException as e:
        print(f"Erreur lors de la requête Semantic Scholar : {e}")
    except Exception as e:
        print(f"Erreur inattendue : {e}")
    return results

def search_openalex(query="machine learning", max_results=5):
    base_url = "https://api.openalex.org/works"
    response = requests.get(base_url, params={"search": query, "per_page": max_results})
    results = []
    if response.status_code == 200:
        for item in response.json().get("results", []):
            doi = item.get("doi")
            if doi:
                unpaywall_url = f"https://api.unpaywall.org/v2/{doi}?email={UNPAYWALL_EMAIL}"
                up_resp = requests.get(unpaywall_url)
                if up_resp.status_code == 200:
                    data = up_resp.json()
                    best_oa_location = data.get("best_oa_location")
                    pdf_url = best_oa_location.get("url_for_pdf") if best_oa_location else None
                    if pdf_url:
                        abstract = item.get("abstract_inverted_index", "")
                        if not isinstance(abstract, str):
                            abstract = ""
                        doc = {
                            "title": item.get("title", "No title"),
                            "url": f"https://doi.org/{doi}",
                            "file_path": pdf_url,
                            "summary": abstract,
                            "source": "OpenAlex"
                        }
                        socketio.emit('document', doc)
                        results.append(doc)
    return results



def merge_and_deduplicate(doc_lists):
    seen_urls = set()
    seen_titles = set()
    merged = []

    for docs in doc_lists:
        for doc in docs:
            url = doc.get('url', '').lower()
            title = doc.get('title', '').lower()
            if url in seen_urls or title in seen_titles:
                continue  # On ignore les doublons
            seen_urls.add(url)
            seen_titles.add(title)
            merged.append(doc)

    # Trie par score (ou 0 si absent) décroissant
    merged.sort(key=lambda d: d.get('score', 0), reverse=True)
    return merged

@app.route('/get_documents', methods=['GET'])
def get_documents():
    query = request.args.get('query', default="machine learning", type=str)
    max_results = int(request.args.get('max_results', 5))

    documents = []
    documents += search_arxiv(query, max_results)
    documents += search_pubmed(query, max_results)
    documents += search_semantic_scholar(query, max_results)
    documents += search_openalex(query, max_results)

    return jsonify({"documents": documents})
