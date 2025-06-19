import fitz  # PyMuPDF
import pdfplumber
from collections import defaultdict

def extract_pdf_content(filepath):
    doc = fitz.open(filepath)
    text_sections = []
    tables = []
    figures = []

    # 1. Extraction des sections (via les titres — heuristique simple)
    # Ici on récupère le texte page par page et détecte titres par taille police
    for page_num in range(len(doc)):
        page = doc.load_page(page_num)
        blocks = page.get_text("dict")["blocks"]
        for b in blocks:
            if b['type'] == 0:  # bloc de texte
                for line in b["lines"]:
                    for span in line["spans"]:
                        # Heuristique : si taille police > 14 -> titre / section
                        if span["size"] > 14:
                            text_sections.append(f"SECTION TITLE: {span['text']}")
                        else:
                            text_sections.append(span["text"])

    # 2. Extraction des tableaux avec pdfplumber (plus fiable que PyMuPDF)
    with pdfplumber.open(filepath) as pdf:
        for page in pdf.pages:
            page_tables = page.extract_tables()
            for table in page_tables:
                tables.append(table)

    # 3. Extraction des figures (images)
    for page_num in range(len(doc)):
        page = doc.load_page(page_num)
        images = page.get_images(full=True)
        for img in images:
            xref = img[0]
            base_image = doc.extract_image(xref)
            image_bytes = base_image["image"]
            image_ext = base_image["ext"]
            # Ici tu peux sauvegarder l’image si besoin, ou juste compter
            figures.append({"page": page_num, "ext": image_ext, "size": len(image_bytes)})

    # Préparer une synthèse (exemple simple)
    extracted = {
        "text_sections": "\n".join(text_sections),
        "tables_count": len(tables),
        "figures_count": len(figures),
        "tables": tables,  # optionnel à stocker ou analyser après
        "figures": figures
    }
    return extracted



def chunk_text(text, max_tokens=512):
    words = text.split()
    chunks = []
    current_chunk = []
    current_len = 0
    for word in words:
        current_len += len(tokenizer.tokenize(word))
        if current_len > max_tokens:
            chunks.append(" ".join(current_chunk))
            current_chunk = [word]
            current_len = len(tokenizer.tokenize(word))
        else:
            current_chunk.append(word)
    if current_chunk:
        chunks.append(" ".join(current_chunk))
    return chunks

from transformers import pipeline, AutoTokenizer

summarizer = pipeline("summarization", model="facebook/bart-large-cnn")
tokenizer = AutoTokenizer.from_pretrained("facebook/bart-large-cnn")

def summarize_text(text):
    chunks = chunk_text(text, max_tokens=512)
    summaries = []
    for chunk in chunks:
        output = summarizer(chunk, max_length=150, min_length=40, do_sample=False)
        summaries.append(output[0]["summary_text"])
    return " ".join(summaries)

def generate_table_of_contents(text):
    chunks = chunk_text(text, max_tokens=512)
    toc_parts = []
    for chunk in chunks:
        prompt = f"Génère une table des matières structurée à partir de ce texte (français) :\n\n{chunk}"
        output = summarizer(prompt, max_length=200, min_length=50, do_sample=False)
        toc_parts.append(output[0]["summary_text"])
    return "\n".join(toc_parts)

def extract_key_passages(text):
    prompt = f"Extrait les passages les plus pertinents (3 à 5) du texte suivant, en français :\n\n{text[:4000]}"
    output = summarizer(prompt, max_length=250, min_length=60, do_sample=False)
    return output[0]["summary_text"]    