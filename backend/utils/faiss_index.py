import os
from langchain.vectorstores import FAISS
from langchain_ollama import OllamaEmbeddings
from langchain.schema import Document
import logging
from langchain.embeddings.base import Embeddings
logger = logging.getLogger(__name__)


from bson import ObjectId
from datetime import datetime
from extensions import mongo


class CachedEmbeddingModel(Embeddings):
    def __init__(self, user_id):
        self.user_id = user_id
        self.base_model = OllamaEmbeddings(model="nomic-embed-text")

    def embed_query(self, text: str) -> list[float]:
        return get_cached_embedding(text, self.user_id, embedding_type="query")

    def embed_documents(self, texts: list[str]) -> list[list[float]]:
        # Tu peux aussi cacher ici si besoin, ou simplement faire en boucle
        return [self.embed_query(t) for t in texts]

class CachedFAISS(FAISS):
    def __init__(self, user_id, *args, **kwargs):
        embedding = CachedEmbeddingModel(user_id)
        super().__init__(*args, embedding=embedding, **kwargs)


def get_cached_embedding(text, uid, embedding_type="query"):
    user_id = ObjectId(uid) if not isinstance(uid, ObjectId) else uid
    key = {"user_id": user_id, "text": text, "type": embedding_type}

    cached = mongo.db.embeddings_cache.find_one(key)
    if cached:
        return cached["vector"]

    embedding_model = OllamaEmbeddings(model="nomic-embed-text")
    vector = embedding_model.embed_query(text)

    if not isinstance(vector, list):
        vector = vector.tolist() if hasattr(vector, "tolist") else list(vector)

    mongo.db.embeddings_cache.insert_one({
        **key,
        "vector": vector,
        "created_at": datetime.utcnow()
    })

    return vector

def delete_from_index(doc_path, user_id):
    """
    Supprime un document FAISS en fonction de son chemin (doc_path) et du user_id
    """
    index_path = os.path.join("documents", str(user_id), "faiss_index")

    if not os.path.exists(index_path):
        print(f"Index FAISS non trouvé pour l'utilisateur {user_id}")
        return False

    # Initialiser l'embedding
    embedding = OllamaEmbeddings(model="nomic-embed-text")

    try:
        # Charger l'index existant
        store = FAISS.load_local(index_path, embedding, allow_dangerous_deserialization=True)


        # Supprimer le document avec ce chemin
        initial_count = len(store.docstore._dict)

        # Garder seulement les documents dont le metadata source est différent
        keys_to_keep = [
            k for k, doc in store.docstore._dict.items()
            if doc.metadata.get("source") != doc_path
        ]

        # Mise à jour de l'index avec les documents à garder
        new_docstore = {k: store.docstore._dict[k] for k in keys_to_keep}
        store.docstore._dict = new_docstore
        store.index_to_docstore_id = {i: k for i, k in enumerate(new_docstore.keys())}

        # Sauvegarde
        store.save_local(index_path)
        print(f"Document supprimé de l'index FAISS: {doc_path}")

        return len(new_docstore) < initial_count  # True si suppression faite
    except Exception as e:
        print(f"Erreur suppression FAISS: {e}")
        return False
    
def update_faiss_index(index_folder, documents, user_id):
    try:
        embedding = CachedEmbeddingModel(user_id)

        # Create Langchain Document objects
        langchain_docs = []
        for doc in documents:
            langchain_docs.append(Document(
                page_content=doc["content"],
                metadata=doc["metadata"]
            ))

        index_file = os.path.join(index_folder, "index.faiss")
        if not os.path.exists(index_file):
            logger.info(f"Creating new FAISS index at: {index_file}")
            store = FAISS.from_documents(
                documents=langchain_docs,
                embedding=embedding
            )
            store.save_local(index_folder)
            return True

        logger.info(f"Loading existing FAISS index from: {index_file}")
        store = FAISS.load_local(index_folder, embedding, allow_dangerous_deserialization=True)

        logger.info(f"Adding {len(langchain_docs)} documents to index")
        store.add_documents(langchain_docs)

        store.save_local(index_folder)
        logger.info(f"Index updated successfully at: {index_file}")
        return True

    except Exception as e:
        logger.error(f"Index update error: {str(e)}")
        return False
