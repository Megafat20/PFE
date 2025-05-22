import os
from langchain.vectorstores import FAISS
from langchain_ollama import OllamaEmbeddings

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