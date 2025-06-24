import os
import requests
import smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from pymongo import MongoClient
from bson.objectid import ObjectId
from data_ingestion import search_arxiv,search_openalex,search_pubmed,merge_and_deduplicate,search_core
# --- CONFIGURATION ---

MONGO_URI = os.getenv("MONGO_URI", "mongodb://localhost:27017/PFE")
DB_NAME = "PFE"

SMTP_SERVER = "smtp.gmail.com"
SMTP_PORT = 587
SMTP_USER = "superfataou13@gmail.com"
SMTP_PASSWORD = "moor jmxf ivfd fhvv"
FROM_EMAIL = SMTP_USER

# Date depuis laquelle chercher les articles (ex: 1 jour avant maintenant)
SINCE_DAYS = 1

# --- Connexion MongoDB ---

client = MongoClient(MONGO_URI)
db = client[DB_NAME]

# --- Fonction pour envoyer un email ---

def send_email(to_email, subject, body):
    msg = MIMEText(body, "plain", "utf-8")
    msg["From"] = FROM_EMAIL
    msg["To"] = to_email
    msg["Subject"] = subject

    try:
        with smtplib.SMTP(SMTP_SERVER, SMTP_PORT) as server:
            server.starttls()
            server.login(SMTP_USER, SMTP_PASSWORD)
            print("Login SMTP OK")
            server.sendmail(FROM_EMAIL, [to_email], msg.as_string())
            print(f"Email envoyé à {to_email}")
    except Exception as e:
        print("Erreur login SMTP:", e)

# --- Fonction pour rechercher des articles Arxiv selon un mot-clé ---

def search_arxiv(keyword, max_results=5, since_days=1):
    base_url = "http://export.arxiv.org/api/query"
    date_from = (datetime.utcnow() - timedelta(days=since_days)).strftime("%Y%m%d%H%M%S")
    query = f"all:{keyword}"

    params = {
        "search_query": query,
        "start": 0,
        "max_results": max_results,
        "sortBy": "submittedDate",
        "sortOrder": "descending",
    }

    try:
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        return response.text  # Retour XML à parser
    except Exception as e:
        print(f"Erreur recherche Arxiv pour '{keyword}': {e}")
        return None

# --- Fonction simple pour parser le XML Arxiv (extraction titre + lien + date) ---

import xml.etree.ElementTree as ET

def parse_arxiv_xml(xml_text):
    ns = {"atom": "http://www.w3.org/2005/Atom"}
    root = ET.fromstring(xml_text)
    entries = []
    for entry in root.findall("atom:entry", ns):
        title = entry.find("atom:title", ns).text.strip()
        link = entry.find("atom:id", ns).text.strip()
        published = entry.find("atom:published", ns).text.strip()
        entries.append({"title": title, "link": link, "published": published})
    return entries

# --- Fonction principale d’automatisation ---

def run_notification_job():
    users = list(db.users.find({}, {"email": 1, "interests": 1}))
    print(f"{len(users)} utilisateurs trouvés")

    for user in users:
        email = user.get("email")
        interests = user.get("interests", [])
        if not email or not interests:
            continue

        new_notifications = []

        for keyword in interests:
            # Appeler toutes les sources
            xml_result = search_arxiv(keyword, max_results=3, since_days=SINCE_DAYS)
            if not xml_result:
                continue
            try:
                articles = parse_arxiv_xml(xml_result)
            except ET.ParseError as e:
                print(f"Erreur de parsing XML Arxiv pour '{keyword}': {e}")
                continue
            pubmed_results = search_pubmed(keyword, max_results=3)
            openalex_results = search_openalex(keyword, max_results=3)
            core_results = search_core(keyword, max_results=3)

            # Fusionner et dédupliquer (utilise ta fonction merge_and_deduplicate)
            all_results = merge_and_deduplicate(articles, pubmed_results, openalex_results, core_results)

            for art in all_results:
                # Vérifier si notification existe déjà pour ce lien et utilisateur
                exists = db.notifications.find_one({
                    "user_id": user["_id"],
                    "link": art.get("url") or art.get("link"),
                    "type": "nouvel_article",
                })
                if exists:
                    continue

                # Créer notification
                notif = {
                    "user_id": user["_id"],
                    "type": "nouvel_article",
                    "title": f"Nouveau article sur '{keyword}': {art['title']}",
                    "message": art.get("summary", "Pas de résumé disponible"),
                    "link": art.get("url") or art.get("link"),
                    "read": False,
                    "created_at": datetime.utcnow(),
                }
                notif_id = db.notifications.insert_one(notif).inserted_id
                new_notifications.append(notif)

        # Envoi d’un email résumé si nouvelles notifications
        if new_notifications:
            body_lines = [f"{n['title']}\nLien : {n['link']}" for n in new_notifications]
            body = "\n\n".join(body_lines)
            subject = f"[Assistant IA] {len(new_notifications)} Nouveaux articles pour vous"
            send_email(email, subject, body)


# --- Script principal ---

if __name__ == "__main__":
    print("Lancement du job de notification...")
    run_notification_job()
    print("Job terminé.")
