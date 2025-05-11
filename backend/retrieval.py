import os
from llama_index.core.query_engine import RetrieverQueryEngine
from llama_index.core.retrievers import VectorIndexRetriever, QueryFusionRetriever
from llama_index.retrievers.bm25 import BM25Retriever
from llama_index.core import Settings
from llama_index.core.postprocessor import SimilarityPostprocessor
from llama_index.core.query_pipeline import QueryPipeline
from llama_index.core.schema import QueryBundle
from llama_index.core.query_engine import CustomQueryEngine




# 1. Configuration de la recherche hybride
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

# 2. Traitement de la requête utilisateur (pré-traitement simple ici, à étendre avec NLP)
def preprocess_query(user_query: str):
    # Ici, tu peux ajouter de la tokenization, suppression de stopwords, expansion de requête...
    return QueryBundle(user_query.strip().lower())

# 3. Fonction de recherche principale
def search(query_engine, user_query: str):
    query_bundle = preprocess_query(user_query)
    response = query_engine.query(query_bundle)
    return response
