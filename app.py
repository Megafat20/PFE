# app.py

import streamlit as st
import os
from pathlib import Path
from ingest import (
    fetch_arxiv_papers,
    fetch_pubmed_articles,
    load_documents,
    preprocess_documents,
    build_vector_index,
)

DATA_DIR = "./data"
USER_UPLOAD_DIR = "./data/user_uploads"
PERSIST_DIR = "./storage"

st.set_page_config(page_title="AI Research Assistant", layout="wide")
st.title("🔬 Assistant AI pour la recherche scientifique")

# Interface en deux onglets : 1. Upload utilisateur | 2. Recherche externe
tabs = st.tabs(["📁 Utiliser mes propres documents", "🔎 Recherche de sources externes"])

# ============ TAB 1: USER UPLOAD ============
with tabs[0]:
    st.subheader("📁 Téléchargez et indexez vos propres documents")
    uploaded_files = st.file_uploader(
        "Upload files (PDF, TXT, DOCX, PPTX, HTML, Markdown):",
        type=["pdf", "txt", "docx", "pptx", "md", "html"],
        accept_multiple_files=True
    )

    if uploaded_files:
        Path(USER_UPLOAD_DIR).mkdir(parents=True, exist_ok=True)
        for file in uploaded_files:
            file_path = os.path.join(USER_UPLOAD_DIR, file.name)
            with open(file_path, "wb") as f:
                f.write(file.getbuffer())
        st.success(f"{len(uploaded_files)} fichiers téléchargés avec succès!")

        if st.button("⚙️ Créer l'index de mon document"):
            with st.spinner("Traitement et indexation de vos documents..."):
                documents = load_documents(USER_UPLOAD_DIR)
                if documents:
                    nodes = preprocess_documents(documents)
                    index = build_vector_index(nodes, persist_dir=PERSIST_DIR + "/user")
                    st.session_state["user_query_engine"] = index.as_query_engine()
                    st.success("✅ L'index de votre document est prêt!")
                else:
                    st.warning("⚠️ Aucun document valide n'a été trouvé pour le traitement.")

    if "user_query_engine" in st.session_state:
        st.markdown("---")
        st.subheader("💬 Ask a Question About Your Documents")
        user_question = st.text_input("Your question:")
        if st.button("📖 Get Answer", key="user_answer"):
            with st.spinner("Generating answer..."):
                response = st.session_state["user_query_engine"].query(user_question)
                st.markdown("### 🧠 Response:")
                st.write(response.response)

# ============ TAB 2: EXTERNAL SOURCES ============
with tabs[1]:
    st.subheader("🔍 Recherche et indexation de documents de recherche externes")
    query = st.text_input("Enter your research topic (e.g. deep learning in medicine):", "large language models in healthcare")

    if st.button("🔄 Recherche et indexation de sources externes"):
        with st.spinner("Récupérer les articles d'ArXiv et de PubMed..."):
            fetch_arxiv_papers(query=query, max_results=5)
            fetch_pubmed_articles(query=query, max_results=5)

        with st.spinner("Traitement et indexation des documents..."):
            documents = load_documents(DATA_DIR)
            if documents:
                nodes = preprocess_documents(documents)
                index = build_vector_index(nodes, persist_dir=PERSIST_DIR + "/external")
                st.session_state["external_query_engine"] = index.as_query_engine()
                st.success("✅ External research index created!")
            else:
                st.warning("⚠️ Aucun document n'a été trouvé pour le sujet donné.")

    if "external_query_engine" in st.session_state:
        st.markdown("---")
        st.subheader("💬 Ask a Question About External Research")
        research_question = st.text_input("Your question:")
        if st.button("📖 Get Answer", key="external_answer"):
            with st.spinner("Generating answer..."):
                response = st.session_state["external_query_engine"].query(research_question)
                st.markdown("### 🧠 Response:")
                st.write(response.response)
