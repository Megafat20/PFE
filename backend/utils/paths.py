import os

def get_user_folder(base_folder, user_id):
    # Convertir en string et nettoyer l'ID utilisateur
    user_id = str(user_id).replace("ObjectId('", "").replace("')", "")
    return os.path.join(base_folder, user_id)

def get_user_index_folder(base_folder, user_id):
    user_id = str(user_id).replace("ObjectId('", "").replace("')", "")
    return os.path.join(base_folder, "faiss_indexes", user_id)