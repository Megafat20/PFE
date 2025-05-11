import os
from pathlib import Path
import torch
import arxiv
from Bio import Entrez
import faiss
import numpy as np
from pymongo import MongoClient
from llama_index.core import VectorStoreIndex, SimpleDirectoryReader, StorageContext
from llama_index.core.node_parser import SentenceSplitter
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.core import Settings
from llama_index.llms.ollama import Ollama
from llama_index.core.query_engine import RetrieverQueryEngine
from llama_index.core.retrievers import VectorIndexRetriever, QueryFusionRetriever
from llama_index.retrievers.bm25 import BM25Retriever
from llama_index.core.postprocessor import SimilarityPostprocessor
from llama_index.core.schema import QueryBundle

# 1. Setup LLM - Use quantized model for faster inference
Settings.llm = Ollama(model="mistral:7b-instruct", request_timeout=3000.0)

# 2. Setup Entrez for PubMed
Entrez.email = "yayahissein30@gmail.com"  # Required by NCBI




# Connect to MongoDB
client = MongoClient("mongodb://localhost:27017/")
db = client["research_db"]
collection = db["documents"]

# Setup embeddings
embedding_model = HuggingFaceEmbedding(model_name="all-MiniLM-L6-v2")

# Function to build FAISS index
def build_faiss_index(nodes):
    # Convert documents to embeddings
    embeddings = np.array([embedding_model.embed(node.text) for node in nodes])
    
    # Create FAISS index
    dimension = embeddings.shape[1]  # Dimensionality of the embedding
    index = faiss.IndexFlatL2(dimension)  # Use L2 distance metric
    index.add(embeddings)  # Add embeddings to the index
    
    return index

# Function to store document metadata in MongoDB
def store_metadata_in_mongo(documents, embeddings):
    for doc, embedding in zip(documents, embeddings):
        collection.insert_one({
            "title": doc.title,
            "author": doc.author,
            "embedding": embedding.tolist()
        })

# Preprocessing and building index with FAISS and MongoDB
def process_documents_for_faiss(documents):
    nodes = preprocess_documents(documents)
    embeddings = np.array([embedding_model.embed(node.text) for node in nodes])
    
    # Build FAISS index
    faiss_index = build_faiss_index(nodes)
    
    # Store metadata in MongoDB
    store_metadata_in_mongo(documents, embeddings)
    
    return faiss_index
# 3. Fetch Arxiv papers
def fetch_arxiv_papers(query="machine learning", max_results=2, save_dir="./data/arxiv"):
    os.makedirs(save_dir, exist_ok=True)
    search = arxiv.Search(query=query, max_results=max_results)
    for result in search.results():
        paper_path = os.path.join(save_dir, f"{result.get_short_id()}.pdf")
        if not os.path.exists(paper_path):
            result.download_pdf(filename=paper_path)
    print(f"Fetched {max_results} Arxiv papers to {save_dir}")

# 4. Fetch PubMed abstracts
def fetch_pubmed_articles(query="machine learning", max_results=2, save_dir="./data/pubmed"):
    os.makedirs(save_dir, exist_ok=True)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]

    for pubmed_id in ids:
        abstract_path = os.path.join(save_dir, f"pubmed_{pubmed_id}.txt")
        if not os.path.exists(abstract_path):
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="abstract", retmode="text")
            abstract = handle.read()
            with open(abstract_path, "w", encoding="utf-8") as f:
                f.write(abstract)
    print(f"Fetched {len(ids)} PubMed abstracts to {save_dir}")

# 5. Load documents
def load_documents(data_dir: str):
    reader = SimpleDirectoryReader(
        input_dir=data_dir,
        recursive=True,
        required_exts=[".pdf", ".txt", ".pptx", ".docx", ".html", ".md"]
    )
    documents = reader.load_data()
    print(f"Loaded {len(documents)} documents")
    return documents

# 6. Preprocess into chunks/nodes
def preprocess_documents(documents):
    splitter = SentenceSplitter(
        chunk_size=1024,
        chunk_overlap=128,
        paragraph_separator="\n\n",
        secondary_chunking_regex="[^,.;。]+[,.;。]?"
    )
    nodes = splitter.get_nodes_from_documents(documents)
    print(f"Split into {len(nodes)} nodes/chunks")
    return nodes

# 7. Setup embeddings
def setup_embeddings(model_name: str = "all-MiniLM-L6-v2"):
    return HuggingFaceEmbedding(
        model_name=model_name,
        device="cuda" if torch.cuda.is_available() else "cpu",
    )

# 8. Build vector index
def build_vector_index(nodes, persist_dir: str = "./storage"):
    Settings.embed_model = setup_embeddings()
    Settings.chunk_size = 1024
    Settings.context_window = 4096

    storage_context = StorageContext.from_defaults()
    index = VectorStoreIndex(nodes, storage_context=storage_context)
    index.storage_context.persist(persist_dir=persist_dir)
    print(f"Index persisted to {persist_dir}")
    return index

# 9. Configuration de la recherche hybride
def create_hybrid_query_engine(index, similarity_top_k=5, bm25_top_k=5):
    # Recherche sémantique
    vector_retriever = VectorIndexRetriever(index=index, similarity_top_k=similarity_top_k)
    
    # Recherche lexicale (BM25)
    bm25_retriever = BM25Retriever.from_defaults(docstore=index.docstore, similarity_top_k=bm25_top_k)

    # Fusion des résultats
    retriever = QueryFusionRetriever(retrievers=[bm25_retriever, vector_retriever], mode="reciprocal_rerank")
    
    # Post-traitement (filtrage par score de similarité)
    postprocessors = [SimilarityPostprocessor(similarity_cutoff=0.7)]

    # Création de l'engine de requêtes
    query_engine = RetrieverQueryEngine.from_args(
        retriever=retriever,
        postprocessors=postprocessors,
    )

    return query_engine

# 10. Traitement de la requête utilisateur
def preprocess_query(user_query: str):
    return QueryBundle(user_query.strip().lower())

def search(query_engine, user_query: str):
    query_bundle = preprocess_query(user_query)
    response = query_engine.query(query_bundle)
    return response

# === Main pipeline ===
if __name__ == "__main__":
    DATA_DIR = "./data"
    PERSIST_DIR = "./storage"

    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
    Path(PERSIST_DIR).mkdir(parents=True, exist_ok=True)

    # Fetch data
    fetch_arxiv_papers(query="large language models", max_results=2)
    fetch_pubmed_articles(query="deep learning for healthcare", max_results=2)

    # Load documents and build index
    if not os.path.exists(os.path.join(PERSIST_DIR, "vector_store.json")):
        documents = load_documents(DATA_DIR)
        if documents:
            nodes = preprocess_documents(documents)
            index = build_vector_index(nodes, PERSIST_DIR)
        else:
            print("No documents found in the data directory.")
    else:
        print("Index already exists. Loading from disk.")
        index = VectorStoreIndex.from_persist_dir(PERSIST_DIR)

    # Configure hybrid search
    query_engine = create_hybrid_query_engine(index)

    # Test query
    response = search(query_engine, "Summarize recent advances in LLMs for healthcare")
    print("\nTest Query Response:\n", response)
