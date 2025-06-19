import os
from langchain_community.vectorstores import FAISS
from langchain_ollama import OllamaEmbeddings
from langchain.schema import Document
import logging

logger = logging.getLogger(__name__)

def get_user_index_folder(base_folder, user_id):
    return os.path.join(base_folder, "faiss_indexes", str(user_id))

def update_faiss_index(index_folder, documents):
    try:
        embedding = OllamaEmbeddings(model="nomic-embed-text")
        
        # Create Langchain Document objects
        langchain_docs = []
        for doc in documents:
            # Create Langchain Document from our document format
            langchain_docs.append(Document(
                page_content=doc["content"],
                metadata=doc["metadata"]
            ))
        
        # Create new index if none exists
        index_file = os.path.join(index_folder, "index.faiss")
        if not os.path.exists(index_file):
            logger.info(f"Creating new FAISS index at: {index_file}")
            # Create and save new index
            store = FAISS.from_documents(
                documents=langchain_docs,
                embedding=embedding
            )
            store.save_local(index_folder)
            return True
        
        # Load existing index
        logger.info(f"Loading existing FAISS index from: {index_file}")
        store = FAISS.load_local(index_folder, embedding, allow_dangerous_deserialization=True)
        
        # Add new documents
        logger.info(f"Adding {len(langchain_docs)} documents to index")
        store.add_documents(langchain_docs)
        
        # Save updated index
        store.save_local(index_folder)
        logger.info(f"Index updated successfully at: {index_file}")
        return True
    
    except Exception as e:
        logger.error(f"Index update error: {str(e)}")
        return False

def delete_from_index(doc_ids, user_id, base_folder):
    try:
        index_folder = get_user_index_folder(base_folder, str(user_id))
        index_file = os.path.join(index_folder, "index.faiss")
        if not os.path.exists(index_file):
            logger.warning(f"No index found at: {index_file}")
            return False

        embedding = OllamaEmbeddings(model="nomic-embed-text")
        store = FAISS.load_local(index_folder, embedding, allow_dangerous_deserialization=True)
        
        # Delete documents from index
        if doc_ids:
            logger.info(f"Deleting {len(doc_ids)} documents from index")
            store.delete(doc_ids)
            store.save_local(index_folder)
        
        return True
    except Exception as e:
        logger.error(f"Index deletion error: {str(e)}")
        return False