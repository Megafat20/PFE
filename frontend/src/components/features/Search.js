import React, { useState } from "react";
import { Search as SearchIcon, Globe, Loader2 } from "lucide-react";

// Fonction simple de distance cosinus pour similitude (pour regroupement)
function cosineSimilarity(vecA, vecB) {
  const dotProduct = vecA.reduce((sum, a, i) => sum + a * vecB[i], 0);
  const magA = Math.sqrt(vecA.reduce((sum, a) => sum + a * a, 0));
  const magB = Math.sqrt(vecB.reduce((sum, b) => sum + b * b, 0));
  return dotProduct / (magA * magB);
}

// Fonction simple de création d'embeddings dummy (à remplacer par ton API)
// Ici on fait juste une transformation en vecteur de longueur fixe avec char codes, pour l'exemple
function textToVector(text, length = 50) {
  const vec = new Array(length).fill(0);
  for (let i = 0; i < Math.min(text.length, length); i++) {
    vec[i] = text.charCodeAt(i) / 255;
  }
  return vec;
}

const Search = () => {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState([]);
  const [translatedResults, setTranslatedResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [translating, setTranslating] = useState(false);
  const [error, setError] = useState(null);
  const [language, setLanguage] = useState("fr");
  const [modalContent, setModalContent] = useState(null);
  const [summary, setSummary] = useState(""); // Synthèse globale
  const [citations, setCitations] = useState([]); // Liste citations simples
  const token = localStorage.getItem("authToken");

  // Recherche classique
  const handleSearch = async () => {
    if (!query.trim()) return;
    setLoading(true);
    setError(null);
    setResults([]);
    setTranslatedResults([]);
    setSummary("");
    setCitations([]);

    try {
      const response = await fetch("http://localhost:5000/multi/search", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          query: query.trim(),
          max_results: 10,
          llm: true,
        }),
      });

      const data = await response.json();

      if (!response.ok) {
        setError(data.error || data.search_error || "Erreur lors de la recherche");
      } else {
        setResults(data.documents || []);
        setTranslatedResults(data.documents || []);
        setSummary(data.summary || "");
        // Préparer citations simples
        const cits = (data.documents || []).map(
          (doc, i) => `${i + 1}. ${doc.title}${doc.authors ? " - " + doc.authors.join(", ") : ""}`
        );
        setCitations(cits);
      }
    } catch (err) {
      setError("Erreur de communication avec le serveur");
    } finally {
      setLoading(false);
    }
  };

  // Traduction des résultats
  const translateResults = async (lang) => {
    if (results.length === 0) return;
    setTranslating(true);

    const texts = results.map(
      (doc) => `${doc.title}\n\n${doc.summary || doc.abstract || ""}`
    );

    try {
      const res = await fetch("http://localhost:5000/multi/translate_bulk", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ texts, language: lang }),
      });

      const data = await res.json();
      if (data.translations) {
        const newResults = results.map((doc, i) => ({
          ...doc,
          translatedText: data.translations[i],
        }));
        setTranslatedResults(newResults);
      } else {
        alert("Erreur dans la traduction");
      }
    } catch (err) {
      alert("Erreur serveur: " + err.message);
    }

    setTranslating(false);
  };

  // Synthèse globale / Literature Review
  const handleLiteratureReview = async () => {
    if (results.length === 0) {
      alert("Aucun document à synthétiser");
      return;
    }
    setLoading(true);
    setError(null);

    try {
      const combinedText = results
        .map((doc) => doc.summary || doc.abstract || "")
        .join("\n\n");

      const response = await fetch("http://localhost:5000/multi/summarize", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ text: combinedText }),
      });

      const data = await response.json();

      if (!response.ok) {
        setError(data.error || "Erreur lors de la synthèse");
      } else {
        setSummary(data.summary || "Pas de synthèse disponible");
      }
    } catch (err) {
      setError("Erreur serveur lors de la synthèse");
    } finally {
      setLoading(false);
    }
  };

  // Regroupement simple par similarité (démonstration, fait sur résumé)
  const groupedResults = () => {
    if (translatedResults.length <= 1) return [translatedResults];

    // Calcul embeddings simples
    const embeddings = translatedResults.map((doc) =>
      textToVector(doc.translatedText || doc.summary || doc.abstract || "")
    );

    // Regroupement basique : on crée des clusters par similarité > 0.8
    const clusters = [];
    const assigned = new Array(translatedResults.length).fill(false);

    for (let i = 0; i < translatedResults.length; i++) {
      if (assigned[i]) continue;
      let cluster = [translatedResults[i]];
      assigned[i] = true;

      for (let j = i + 1; j < translatedResults.length; j++) {
        if (!assigned[j]) {
          const sim = cosineSimilarity(embeddings[i], embeddings[j]);
          if (sim > 0.8) {
            cluster.push(translatedResults[j]);
            assigned[j] = true;
          }
        }
      }
      clusters.push(cluster);
    }

    return clusters;
  };

  const openModal = (summaryText) => {
    setModalContent(summaryText);
  };

  const closeModal = () => {
    setModalContent(null);
  };

  return (
    <div className="max-w-6xl mx-auto px-4 py-8">
      <h1 className="text-4xl font-bold mb-8 text-blue-700 tracking-tight">
        🔎 Recherche scientifique assistée par IA - Literature Review
      </h1>

      {/* Langue + Barre de recherche */}
      <div className="flex flex-col md:flex-row items-center gap-4 mb-6">
        <div className="flex items-center gap-2 w-full md:w-1/3">
          <Globe className="text-blue-600" />
          <select
            value={language}
            onChange={async (e) => {
              const lang = e.target.value;
              setLanguage(lang);
              if (results.length > 0) {
                setTranslating(true);
                await translateResults(lang);
                setTranslating(false);
              }
            }}
            className="w-full border border-gray-300 rounded-lg px-3 py-2 text-gray-700 focus:ring-2 focus:ring-blue-500"
            disabled={loading || translating}
          >
            <option value="">Langue</option>
            <option value="fr">Français</option>
            <option value="en">Anglais</option>
            <option value="es">Espagnol</option>
            <option value="de">Allemand</option>
            <option value="it">Italien</option>
          </select>
        </div>

        <div className="flex-grow w-full md:w-2/3 flex items-center gap-2">
          <input
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && handleSearch()}
            placeholder="Ex: Large Language Models for Scientific Research"
            className="flex-grow border border-gray-300 rounded-lg px-4 py-3 focus:outline-none focus:ring-2 focus:ring-blue-500 text-lg"
            disabled={loading || translating}
          />
          <button
            onClick={handleSearch}
            disabled={loading || translating}
            className="flex items-center gap-2 px-5 py-3 rounded-lg bg-blue-600 hover:bg-blue-700 text-white font-semibold transition disabled:opacity-50"
          >
            {loading ? <Loader2 className="animate-spin w-5 h-5" /> : <SearchIcon />}
            {loading ? "Recherche..." : "Rechercher"}
          </button>
        </div>
      </div>

      {/* Bouton Literature Review / Synthèse globale */}
      <div className="mb-6">
        <button
          onClick={handleLiteratureReview}
          disabled={loading || translating || results.length === 0}
          className="px-6 py-3 rounded-lg bg-green-600 hover:bg-green-700 text-white font-semibold transition disabled:opacity-50"
        >
          Générer une synthèse globale (Literature Review)
        </button>
      </div>

      {/* Erreur */}
      {error && (
        <div className="text-red-600 font-medium mb-4 border-l-4 border-red-600 pl-4 py-2 bg-red-50 rounded">
          {error}
        </div>
      )}

      {/* Synthèse globale */}
      {summary && (
        <div className="mb-8 p-6 bg-gray-100 rounded-lg border border-gray-300 text-gray-800 whitespace-pre-line">
          <h2 className="text-2xl font-semibold mb-4 text-green-700">
            Synthèse globale (Literature Review)
          </h2>
          {summary}
        </div>
      )}

      {/* Citations */}
      {citations.length > 0 && (
        <div className="mb-8 p-6 bg-white rounded-lg border border-gray-300 text-gray-700">
          <h3 className="text-xl font-semibold mb-3 text-blue-700">Citations des documents</h3>
          <ul className="list-decimal list-inside max-h-48 overflow-y-auto">
            {citations.map((cit, i) => (
              <li key={i} className="mb-1">
                {cit}
              </li>
            ))}
          </ul>
        </div>
      )}

      {/* Résultats regroupés */}
      {translatedResults.length > 0 && (
        <div>
          <h2 className="text-3xl font-bold mb-6 text-blue-800">
            Résultats regroupés par similarité thématique
          </h2>
          {groupedResults().map((cluster, idx) => (
            <div key={idx} className="mb-8">
              <h3 className="text-2xl font-semibold mb-4">
                Groupe {idx + 1} ({cluster.length} documents)
              </h3>
              {cluster.map((doc, i) => (
                <div
                  key={i}
                  onClick={() => openModal(doc.translatedText || doc.summary || doc.abstract || "")}
                  className="mb-3 p-4 border rounded-lg cursor-pointer hover:bg-blue-50"
                  title="Cliquez pour voir le résumé complet"
                >
                  <h4 className="text-lg font-semibold text-blue-700">
                    {doc.title}
                  </h4>
                  <p className="text-gray-700 line-clamp-3">
                    {doc.translatedText || doc.summary || doc.abstract}
                  </p>
                </div>
              ))}
            </div>
          ))}
        </div>
      )}

      {/* Modal pour afficher le résumé complet */}
      {modalContent && (
        <div
          onClick={closeModal}
          className="fixed inset-0 bg-black bg-opacity-60 flex justify-center items-center z-50"
        >
          <div
            onClick={(e) => e.stopPropagation()}
            className="max-w-3xl bg-white p-6 rounded-lg shadow-lg max-h-[80vh] overflow-y-auto whitespace-pre-line"
          >
            <button
              onClick={closeModal}
              className="mb-4 font-bold text-red-600 hover:underline"
            >
              Fermer
            </button>
            {modalContent}
          </div>
        </div>
      )}
    </div>
  );
};

export default Search;
