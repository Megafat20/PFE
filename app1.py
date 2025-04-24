# app.py

import streamlit as st
import os
from pathlib import Path
from ingest1 import (
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
st.title("🔬 AI Assistant for Scientific Research")

# Initialiser les mémoires de chat
if "chat_history_user" not in st.session_state:
    st.session_state["chat_history_user"] = []

if "chat_history_external" not in st.session_state:
    st.session_state["chat_history_external"] = []

# Onglets : 1. Documents utilisateur | 2. Sources externes
tabs = st.tabs(["📁 Use My Own Documents", "🔎 Search External Sources"])

# ============ TAB 1: USER UPLOAD ============
with tabs[0]:
    st.subheader("📁 Upload and Index Your Own Documents")
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
        st.success(f"{len(uploaded_files)} files uploaded successfully!")

        if st.button("⚙️ Build My Document Index"):
            with st.spinner("Processing and indexing your documents..."):
                documents = load_documents(USER_UPLOAD_DIR)
                if documents:
                    nodes = preprocess_documents(documents)
                    index = build_vector_index(nodes, persist_dir=PERSIST_DIR + "/user")
                    st.session_state["user_query_engine"] = index.as_query_engine()
                    st.success("✅ Your document index is ready!")
                else:
                    st.warning("⚠️ No valid documents found for processing.")

    if "user_query_engine" in st.session_state:
     st.markdown("---")
     st.subheader("💬 Ask Questions About Your Documents")

     if "user_questions" not in st.session_state:
        st.session_state["user_questions"] = [""]  # Start with one empty input

     for i, question in enumerate(st.session_state["user_questions"]):
        col1, col2 = st.columns([4, 1])
        with col1:
            new_q = st.text_input(f"Your question #{i + 1}:", value=question, key=f"user_q_{i}")
            st.session_state["user_questions"][i] = new_q
        with col2:
            if st.button("📖 Get Answer", key=f"user_answer_{i}"):
                with st.spinner("Generating answer..."):
                    response = st.session_state["user_query_engine"].query(new_q)
                    st.session_state["chat_history_user"].append(("You", new_q))
                    st.session_state["chat_history_user"].append(("Assistant", response.response))

    if st.button("➕ Add Another Question", key="add_user_q"):
        st.session_state["user_questions"].append("")

    # Affichage historique
    if st.session_state["chat_history_user"]:
        st.markdown("### 🗂️ Conversation History")
        for role, message in st.session_state["chat_history_user"]:
            st.markdown(f"**{role}:** {message}")

# ============ TAB 2: EXTERNAL SOURCES ============
with tabs[1]:
    st.subheader("🔍 Search and Index External Research Papers")
    query = st.text_input("Enter your research topic (e.g. deep learning in medicine):", "large language models in healthcare")

    if st.button("🔄 Fetch & Index External Sources"):
        with st.spinner("Fetching papers from ArXiv and PubMed..."):
            fetch_arxiv_papers(query=query, max_results=5)
            fetch_pubmed_articles(query=query, max_results=5)

        with st.spinner("Processing and indexing documents..."):
            documents = load_documents(DATA_DIR)
            if documents:
                nodes = preprocess_documents(documents)
                index = build_vector_index(nodes, persist_dir=PERSIST_DIR + "/external")
                st.session_state["external_query_engine"] = index.as_query_engine()
                st.success("✅ External research index created!")
            else:
                st.warning("⚠️ No documents were found for the given topic.")

    if "external_query_engine" in st.session_state:
        st.markdown("---")
        st.subheader("💬 Ask a Question About External Research")
        research_question = st.text_input("Your question:", key="external_input")

        if st.button("📖 Get Answer", key="external_answer"):
            with st.spinner("Generating answer..."):
                response = st.session_state["external_query_engine"].query(research_question)
                st.session_state["chat_history_external"].append(("You", research_question))
                st.session_state["chat_history_external"].append(("Assistant", response.response))

        # Affichage de l'historique
        if st.session_state["chat_history_external"]:
            st.markdown("### 🗂️ Conversation History")
            for role, message in st.session_state["chat_history_external"]:
                st.markdown(f"**{role}:** {message}")

# Optionnel : bouton pour réinitialiser toute la mémoire
st.sidebar.title("🧠 Memory Settings")
if st.sidebar.button("🧹 Clear All Chat History"):
    st.session_state["chat_history_user"] = []
    st.session_state["chat_history_external"] = []
    st.success("All chat history cleared!")
