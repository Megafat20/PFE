# data_ingestion.py
import os
from flask import Flask, request, jsonify
import arxiv
from flask_socketio import SocketIO, emit
from Bio import Entrez

app = Flask(__name__)
socketio = SocketIO(app)
Entrez.email = "superfataou13@gmail.com"  # Remplace par ton email réel

def fetch_arxiv_papers(query="machine learning", max_results=5, save_dir="./data/arxiv"):
    os.makedirs(save_dir, exist_ok=True)
    search = arxiv.Search(query=query, max_results=max_results)
    results = []
    for result in search.results():
        paper_path = os.path.join(save_dir, f"{result.get_short_id()}.pdf")
        if not os.path.exists(paper_path):
            result.download_pdf(filename=paper_path)
        
        # Emit message en temps réel au frontend avec SocketIO
        socketio.emit('document', {
            "title": result.title,
            "url": result.entry_id,
            "file_path": paper_path
        })
        results.append({
            "title": result.title,
            "url": result.entry_id,
            "file_path": paper_path,
            "summary": result.summary
        })
    return results

def fetch_pubmed_articles(query="machine learning", max_results=5, save_dir="./data/pubmed"):
    os.makedirs(save_dir, exist_ok=True)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]
    summary = "Résumé non disponible pour PubMed."
    for id_ in ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=id_, rettype="xml", retmode="text")
        article_data = fetch_handle.read()
        file_path = os.path.join(save_dir, f"{id_}.xml")
        with open(file_path, "w") as f:
            f.write(article_data)
        
        socketio.emit('document', {
            "title": f"PubMed Article {id_}",
            "file_path": file_path,
            "summary": summary
        })

        return [{
            "title": f"PubMed Article {id_}",
            "file_path": file_path,
            "summary": summary,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{id_}"
        } for id_ in ids]
    
    return [{"id": id_, "file_path": os.path.join(save_dir, f"{id_}.xml")} for id_ in ids]

@app.route('/get_documents', methods=['GET'])
def get_documents():
    # Combine tous les fichiers arxiv et pubmed
    arxiv_files = fetch_arxiv_papers(query="machine learning", max_results=5)
    pubmed_files = fetch_pubmed_articles(query="machine learning", max_results=5)
    all_documents = arxiv_files + pubmed_files

    return jsonify({"documents": all_documents})