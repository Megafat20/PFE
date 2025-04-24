# data_ingestion.py

import os
from pathlib import Path
import torch
import arxiv
from Bio import Entrez
from llama_index.core import VectorStoreIndex, SimpleDirectoryReader, StorageContext
from llama_index.core.node_parser import SentenceSplitter
from llama_index.embeddings.huggingface import HuggingFaceEmbedding
from llama_index.core import Settings
from llama_index.llms.ollama import Ollama

Settings.llm = Ollama(model="deepseek-r1:latest", request_timeout=3000.0)

# Set up Entrez for PubMed
Entrez.email = "yayahissein30@gmail.com"  # Replace with your email (required by NCBI)


def fetch_arxiv_papers(query="machine learning", max_results=5, save_dir="./data/arxiv"):
    os.makedirs(save_dir, exist_ok=True)
    search = arxiv.Search(query=query, max_results=max_results)
    for result in search.results():
        paper_path = os.path.join(save_dir, f"{result.get_short_id()}.pdf")
        result.download_pdf(filename=paper_path)
    print(f"Fetched {max_results} Arxiv papers to {save_dir}")


def fetch_pubmed_articles(query="machine learning", max_results=5, save_dir="./data/pubmed"):
    os.makedirs(save_dir, exist_ok=True)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    ids = record["IdList"]

    for i, pubmed_id in enumerate(ids):
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="abstract", retmode="text")
        abstract = handle.read()
        with open(os.path.join(save_dir, f"pubmed_{pubmed_id}.txt"), "w", encoding="utf-8") as f:
            f.write(abstract)
    print(f"Fetched {len(ids)} PubMed abstracts to {save_dir}")


def load_documents(data_dir: str):
    reader = SimpleDirectoryReader(
        input_dir=data_dir,
        recursive=True,
        required_exts=[".pdf", ".txt", ".pptx", ".docx", ".html", ".md"]
    )
    documents = reader.load_data()
    print(f"Loaded {len(documents)} documents")
    return documents


def preprocess_documents(documents):
    splitter = SentenceSplitter(
        chunk_size=512,
        chunk_overlap=128,
        paragraph_separator="\n\n",
        secondary_chunking_regex="[^,.;。]+[,.;。]?"
    )
    nodes = splitter.get_nodes_from_documents(documents)
    print(f"Split into {len(nodes)} nodes/chunks")
    return nodes


def setup_embeddings(model_name: str = "BAAI/bge-small-en"):
    return HuggingFaceEmbedding(
        model_name=model_name,
        device="cuda" if torch.cuda.is_available() else "cpu",
    )


def build_vector_index(nodes, persist_dir: str = "./storage"):
    Settings.embed_model = setup_embeddings()
    Settings.chunk_size = 512
    Settings.context_window = 4096

    # Build new storage context (DO NOT load from disk)
    storage_context = StorageContext.from_defaults()

    # Build and persist index
    index = VectorStoreIndex(nodes, storage_context=storage_context)
    index.storage_context.persist(persist_dir=persist_dir)
    print(f"Index persisted to {persist_dir}")
    return index



if __name__ == "__main__":
    DATA_DIR = "./data"
    PERSIST_DIR = "./storage"

    # Create directories
    Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
    Path(PERSIST_DIR).mkdir(parents=True, exist_ok=True)

    # Fetch external data
    fetch_arxiv_papers(query="large language models", max_results=5)
    fetch_pubmed_articles(query="deep learning for healthcare", max_results=5)

    # Load, process, and build index
    documents = load_documents(DATA_DIR)
    if documents:
        nodes = preprocess_documents(documents)
        index = build_vector_index(nodes, PERSIST_DIR)

        # Run a test query
        query_engine = index.as_query_engine()
        response = query_engine.query("Summarize recent advances in LLMs for healthcare")
        print("\nTest Query Response:\n", response)
    else:
        print("No documents found in the data directory.")
