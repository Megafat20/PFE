import React, { useState, useEffect, useRef } from "react";
import { Search as SearchIcon, Globe, Loader2 } from "lucide-react";
import ProgressBar from "../ProgressBar";
import ChatModal from "./ChatModal";
import DocumentInsightModal from "./DocumentInsightModal";
import LoadingOverlay from "../LoadingOverlay";
import { toast } from "react-toastify";

// Fonction de similarité cosinus
function cosineSimilarity(vecA, vecB) {
  const dotProduct = vecA.reduce((sum, a, i) => sum + a * vecB[i], 0);
  const magA = Math.sqrt(vecA.reduce((sum, a) => sum + a * a, 0));
  const magB = Math.sqrt(vecB.reduce((sum, b) => sum + b * b, 0));
  return dotProduct / (magA * magB);
}

// Dummy embedding
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
  const [translationProgress, setTranslationProgress] = useState(0);
  const [error, setError] = useState(null);
  const [language, setLanguage] = useState("fr");
  const [modalContent, setModalContent] = useState(null);
  const [summary, setSummary] = useState("");
  const [citations, setCitations] = useState([]);
  const token = localStorage.getItem("authToken");
  const [searchHistory, setSearchHistory] = useState(() => {
    return JSON.parse(localStorage.getItem("searchHistory")) || [];
  });
  const [showHistory, setShowHistory] = useState(false);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [suggestions, setSuggestions] = useState([]);
  const [chatDoc, setChatDoc] = useState(null);
  const [documents, setDocuments] = useState([]);
  const [loadingProgress, setLoadingProgress] = useState(0);
  const [isIndexing, setIsIndexing] = useState(false);
  const [loadingFavorites, setLoadingFavorites] = useState({});
  const [selectedFullText, setSelectedFullText] = useState("");
  const [selectedChunks, setSelectedChunks] = useState([]);
  const [showModal, setShowModal] = useState(false);
  const [loadingMessage, setLoadingMessage] = useState("");
  const [groupedResults, setGroupedResults] = useState([]);
  const allSuggestions = [
    "Large Language Models for Scientific Research",
    "Deep Learning for Medical Imaging",
    "AI in Agriculture",
    "Natural Language Processing in Education",
    "Graph Neural Networks",
    "Explainable AI",
    "Reinforcement Learning",
    "Climate Change Modeling with AI",
    "Biomedical Signal Processing",
    "Federated Learning in Healthcare",
  ];
  const inputRef = useRef(null);

  const updateSearchHistory = async (newQuery) => {
    if (!newQuery.trim()) return;

    // ⚡ LocalStorage
    const updated = [
      newQuery,
      ...searchHistory.filter((q) => q !== newQuery),
    ].slice(0, 5);
    setSearchHistory(updated);
    localStorage.setItem("searchHistory", JSON.stringify(updated));

    // ☁️ Backend (si authToken présent)
    const token = localStorage.getItem("authToken");
    if (token) {
      try {
        await fetch("http://localhost:5000/history", {
          method: "POST",
          headers: {
            Authorization: `Bearer ${token}`,
            "Content-Type": "application/json",
          },
          body: JSON.stringify({ query: newQuery }),
        });
      } catch (err) {
        console.error("Erreur lors de la sauvegarde de l'historique:", err);
      }
    }
  };

  const handleSearch = async (searchTerm = query) => {
    if (!searchTerm.trim()) return;
    setLoading(true);
    setError(null);
    setCitations([]);

    try {
      const response = await fetch("http://localhost:5000/multi/search", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          query: searchTerm.trim(),
          max_results: 10,
          llm: true,
          language,
        }),
      });

      const data = await response.json();

      if (!response.ok) {
        setError(
          data.error || data.search_error || "Erreur lors de la recherche"
        );
      } else {
        const unifiedDocs = (data.documents || []).map((doc) => ({
          ...doc,
          _id: doc._id || doc.url || doc.title,
          translatedTitle: doc.title,
          translatedSummary: doc.summary || doc.abstract || "",
          indexed: doc.indexed || false,
        }));

        setResults(unifiedDocs);
        setDocuments(unifiedDocs);
        setGroupedResults(data.grouped_documents || []);
        setSummary(data.summary || "");

        const cits = (data.documents || []).map(
          (doc, i) =>
            `${i + 1}. ${doc.title}${
              doc.authors ? " - " + doc.authors.join(", ") : ""
            }`
        );
        setCitations(cits);
        await new Promise((r) => setTimeout(r, 1000));

        updateSearchHistory(searchTerm);
        setQuery(searchTerm);
        setShowHistory(false);
        setShowSuggestions(false);
        setDocuments(data.documents || []);
      }
    } catch (err) {
      setError("Erreur de communication avec le serveur");
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (query.length > 1) {
      const filtered = allSuggestions.filter((s) =>
        s.toLowerCase().includes(query.toLowerCase())
      );
      setSuggestions(filtered);
    } else {
      setSuggestions([]);
    }
  }, [query]);

  const onFocus = () => {
    setShowHistory(true);
    setShowSuggestions(true);
  };

  const onBlur = () => {
    setTimeout(() => {
      setShowHistory(false);
      setShowSuggestions(false);
    }, 200);
  };


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

  // const groupedResults = () => {
  //   if (translatedResults.length <= 1) return [translatedResults];

  //   const embeddings = translatedResults.map((doc) =>
  //     textToVector(doc.translatedText || doc.summary || doc.abstract || "")
  //   );

  //   const clusters = [];
  //   const assigned = new Array(translatedResults.length).fill(false);

  //   for (let i = 0; i < translatedResults.length; i++) {
  //     if (assigned[i]) continue;
  //     let cluster = [translatedResults[i]];
  //     assigned[i] = true;

  //     for (let j = i + 1; j < translatedResults.length; j++) {
  //       if (!assigned[j]) {
  //         const sim = cosineSimilarity(embeddings[i], embeddings[j]);
  //         if (sim > 0.8) {
  //           cluster.push(translatedResults[j]);
  //           assigned[j] = true;
  //         }
  //       }
  //     }
  //     clusters.push(cluster);
  //   }

  //   return clusters;
  // };

  const openModal = (summaryText) => {
    setModalContent(summaryText);
  };

  const closeModal = () => {
    setModalContent(null);
  };
  const handleDownloadPDF = async (doc) => {
    let pdfUrl = getPdfUrl(doc);
    if (!pdfUrl) {
      alert("URL PDF non disponible pour ce document.");
      return;
    }
  
    // Encode l'URL pour passer en paramètre
    const proxyUrl = `http://localhost:5000/proxy_pdf?url=${encodeURIComponent(pdfUrl)}`;
  
    try {
      const response = await fetch(proxyUrl);
  
      const contentType = response.headers.get("content-type") || "";
      if (!response.ok || !contentType.includes("pdf")) {
        alert("Le fichier PDF n'est pas disponible ou l'URL est invalide.");
        return;
      }
  
      const blob = await response.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
  
      const fileName = (doc.title || "document").replace(/[^\w\d]+/g, "_") + ".pdf";
      a.download = fileName;
      document.body.appendChild(a);
      a.click();
  
      a.remove();
      window.URL.revokeObjectURL(url);
  
    } catch (error) {
      alert("Erreur lors du téléchargement du PDF.");
      console.error(error);
    }
  };

const getPdfUrl = (doc) => {
  // 1. Priorité au pdf_url s’il existe
  if (doc.pdf_url) return doc.pdf_url;

  // 2. ArXiv : convertir abs en pdf
  if (doc.url && doc.url.includes("arxiv.org/abs/")) {
    return doc.url.replace("arxiv.org/abs/", "arxiv.org/pdf/") + ".pdf";
  }

  // 3. CORE / OpenAlex : ils ont souvent pdf_url déjà géré en 1.

  // 4. PubMed et autres : généralement pas de PDF direct, retourner page
  if (doc.url && (doc.url.includes("pubmed.ncbi.nlm.nih.gov") || doc.url.includes("ncbi.nlm.nih.gov"))) {
    // Pas de PDF direct dispo, retourne la page
    return doc.url;
  }

  // 5. Par défaut, retourne url (page d'article)
  return doc.url || "";
};

  const openInsightModal = (doc, chunks, fullText) => {
    setChatDoc(doc);
    setSelectedChunks(chunks);
    setSelectedFullText(fullText);
    setShowModal(true);
  };

  const closeInsightModal = () => {
    setShowModal(false);
    setChatDoc(null);
    setSelectedChunks([]);
    setSelectedFullText("");
  };

  // 🔄 Fusionne ta logique actuelle dans cette méthode
  const handleQuestionDocument = async (doc) => {
    setIsIndexing(true);
    setLoadingProgress(0);

    try {
      const token = localStorage.getItem("authToken");

      let fileUrl = doc.url;
      let filename =
        (doc.filename || doc.title || "document")
          .replace(/\s+/g, "_")
          .replace(/\.pdf$/, "") + ".pdf";

      if (fileUrl.includes("arxiv.org/abs/")) {
        fileUrl = fileUrl.replace("/abs/", "/pdf/") + ".pdf";
      }

      const simulateProgress = () => {
        setLoadingProgress((prev) => {
          if (prev >= 100) return 100;
          const next = prev + Math.floor(Math.random() * 15) + 5;
          return next > 100 ? 100 : next;
        });
      };
      const progressInterval = setInterval(simulateProgress, 400);

      const res = await fetch(
        "http://localhost:5000/download_external_document",
        {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            Authorization: `Bearer ${token}`,
          },
          body: JSON.stringify({
            file_url: fileUrl,
            filename,
            source: doc.source || "unknown",
            original_url: doc.url,
          }),
        }
      );

      clearInterval(progressInterval);
      setLoadingProgress(100);

      const data = await res.json();
      if (!res.ok || !data.document) {
        alert(data.error || "Erreur lors de l’indexation");
      } else {
        // Ensuite, appelle automatiquement handleShowChunks
        await handleShowChunks(data.document);
      }
    } catch (err) {
      alert("❌ Erreur : " + err.message);
    }

    setIsIndexing(false);
    setLoadingProgress(0);
  };

  const handleShowChunks = async (doc) => {
    const token = localStorage.getItem("authToken");
    if (!token) return alert("Token non trouvé");

    const documentId = doc._id || doc.url;

    try {
      setIsIndexing(true);
      setLoadingProgress(0);
      setLoadingMessage("Préparation du document en cours...");

      const chunksRes = await fetch("http://localhost:5000/document_chunks", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          document_id: documentId,
          filename: doc.title || "document.pdf",
          source: doc.source || "unknown",
          language ,
        }),
      });

      if (!chunksRes.ok) {
        const errData = await chunksRes.json();
        throw new Error(
          errData.error || "Erreur lors de la récupération des extraits"
        );
      }

      const chunksData = await chunksRes.json();

      // Ouvre le modal fusionné avec tous les éléments nécessaires
      openInsightModal(
        doc,
        chunksData.chunks || [],
        chunksData.full_text || ""
      );
    } catch (err) {
      alert(err.message);
    }

    setIsIndexing(false);
    setLoadingProgress(0);
    setLoadingMessage("");
  };

  const handleAddToLibrary = async (doc) => {
    setLoadingFavorites((prev) => ({ ...prev, [doc._id]: true }));
    try {
      const token = localStorage.getItem("authToken");
      const response = await fetch("http://localhost:5000/api/favorites", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          documentId: doc._id || doc.url || doc.title, // 👈 À adapter si _id manquant
          documentData: doc, // 👈 Requis si le doc n'est pas dans MongoDB
        }),
      });

      if (!response.ok) throw new Error("Échec de l'ajout aux favoris.");
      toast.success("Ajouté aux favoris !");
    } catch (error) {
      toast.error("Impossible d'ajouter aux favoris.");
    } finally {
      setLoadingFavorites((prev) => ({ ...prev, [doc._id]: false }));
    }
  };
  return (
    <>
      {(loading || translating) && (
        <ProgressBar
          progress={translating ? translationProgress : loading ? 100 : 0}
        />
      )}
      {isIndexing && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex flex-col justify-center items-center z-50">
          <div className="loader ease-linear rounded-full border-8 border-t-8 border-gray-200 h-20 w-20 mb-4"></div>
          <p className="text-white text-lg">{loadingMessage}</p>
        </div>
      )}
      <div className="w-full px-4 py-6 sm:px-6 lg:px-8 max-w-7xl mx-auto">
        {/* Titre principal */}
        <h1 className="text-2xl sm:text-3xl md:text-4xl font-bold mb-6 sm:mb-8 text-blue-700 tracking-tight text-center sm:text-left">
          🔎 Recherche scientifique assistée par IA - Literature Review
        </h1>

        <div className="w-full font-sans">
          <div className="flex flex-col md:flex-row items-stretch md:items-center gap-4 mb-6">
            {/* Sélecteur de langue */}
            <div className="flex items-center gap-2 w-full md:w-1/4">
              <label htmlFor="language-select" className="sr-only">
                Sélectionner la langue
              </label>
              <button
                type="button"
                className="flex items-center gap-2 px-4 py-2 rounded-md border border-gray-400 text-gray-700 font-medium hover:bg-gray-100 focus:outline-none focus:ring-2 focus:ring-gray-300"
                onClick={() =>
                  document.getElementById("language-select").focus()
                }
              >
                <Globe className="w-5 h-5 text-gray-600" />
              </button>

              <select
                id="language-select"
                value={language}
                onChange={async (e) => {
                  const lang = e.target.value;
                  setLanguage(lang);
                }}
                className="w-full rounded-md border border-gray-400 text-gray-800 px-4 py-2 focus:outline-none focus:ring-2 focus:ring-gray-300"
                disabled={loading || translating}
              >
                <option value="">Langue</option>
                <option value="en">Anglais</option>
                <option value="fr">Français</option>
                <option value="ar">Arabe</option>
                <option value="de">Allemand</option>
                <option value="it">Italien</option>
              </select>
            </div>

            {/* Barre de recherche */}
            <div className="flex-grow w-full md:w-3/4 relative">
              <input
                ref={inputRef}
                type="text"
                value={query}
                onChange={(e) => setQuery(e.target.value)}
                onKeyDown={(e) => e.key === "Enter" && handleSearch()}
                onFocus={onFocus}
                onBlur={onBlur}
                placeholder="Ex: Large Language Models for Scientific Research"
                className="w-full rounded-full border border-gray-300 bg-white px-6 py-4 text-base sm:text-lg text-gray-900 placeholder-gray-400 shadow-md focus:outline-none focus:ring-4 focus:ring-blue-400 focus:border-blue-400 transition"
                aria-label="Champ de recherche"
                disabled={loading}
                autoComplete="off"
              />

              {/* Icone */}
              <div className="absolute right-4 top-1/2 -translate-y-1/2 flex items-center">
                {loading ? (
                  <span className="animate-spin h-6 w-6 border-4 border-blue-500 border-t-transparent rounded-full"></span>
                ) : (
                  <button
                    onClick={() => handleSearch()}
                    disabled={loading || query.trim() === ""}
                    className="text-gray-700 hover:text-gray-900 focus:outline-none"
                    aria-label="Lancer la recherche"
                    title="Lancer la recherche"
                  >
                    <SearchIcon className="w-6 h-6" />
                  </button>
                )}
              </div>

              {/* Suggestions glassmorphism */}
              {showSuggestions && suggestions.length > 0 && (
                <ul
                  className="absolute top-full left-0 right-0 mt-3 bg-white/60 backdrop-blur-sm border border-blue-300 rounded-2xl shadow-lg max-h-72 overflow-y-auto scrollbar-thin scrollbar-thumb-blue-400 scrollbar-track-blue-100 z-50"
                  role="listbox"
                  aria-label="Suggestions de recherche"
                >
                  {suggestions.map((s, i) => (
                    <li
                      key={i}
                      role="option"
                      aria-selected="false"
                      className="px-6 py-3 hover:bg-blue-100 cursor-pointer text-blue-900 font-medium select-none transition"
                      onMouseDown={() => handleSearch(s)}
                    >
                      🔍 {s}
                    </li>
                  ))}
                </ul>
              )}

              {/* Historique glassmorphism */}
              {showHistory && searchHistory.length > 0 && (
                <div className="absolute top-full left-0 right-0 mt-3 bg-white/50 backdrop-blur-md border border-blue-200 rounded-2xl shadow-xl max-h-72 overflow-y-auto scrollbar-thin scrollbar-thumb-blue-300 scrollbar-track-blue-50 z-40">
                  <div className="p-3 text-sm text-blue-600 font-semibold border-b border-blue-300 select-none">
                    🕘 Historique de recherche
                  </div>
                  <ul>
                    {searchHistory.map((item, index) => (
                      <li
                        key={index}
                        onMouseDown={() => {
                          setQuery(item);
                          handleSearch(item);
                        }}
                        className="cursor-pointer px-6 py-3 hover:bg-blue-200 text-blue-800 font-semibold select-none transition"
                      >
                        {item}
                      </li>
                    ))}
                  </ul>
                </div>
              )}
            </div>
          </div>
        </div>

        {/* Bouton Literature Review */}
        <div className="mb-6">
          <button
            onClick={handleLiteratureReview}
            disabled={loading || translating || results.length === 0}
            className="px-6 py-3 rounded-md bg-neutral-800 text-white hover:bg-neutral-700 font-medium transition disabled:opacity-50"
            title="Analyser les résultats et générer une synthèse"
          >
            Générer une synthèse globale
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
            <h3 className="text-xl font-semibold mb-3 text-blue-700">
              Citations des documents
            </h3>
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
        {groupedResults.length > 0 && (
          <div className="mt-10">
            <h2 className="text-3xl font-bold mb-6 text-blue-800">
              Résultats regroupés par thème
            </h2>

            {groupedResults.map((group, idx) => (
              <div key={idx} className="mb-10">
                <h3 className="text-2xl font-semibold mb-4 text-indigo-700">
                  🧠 {group.theme} ({group.documents.length} articles)
                </h3>

                {group.documents.map((doc, i) => (
                  <div
                    key={doc._id || doc.url || doc.title || i}
                    onClick={() =>
                      openModal(
                        doc.translatedText || doc.summary || doc.abstract || ""
                      )
                    }
                    className="mb-4 p-5 border rounded-xl hover:bg-blue-50 hover:shadow-md transition cursor-pointer flex flex-col md:flex-row justify-between items-start gap-4"
                    title="Cliquez pour voir le résumé complet"
                  >
                    {/* Partie gauche : texte */}
                    <div className="flex-1 min-w-0">
                      <h4 className="text-xl font-bold text-gray-900 mb-2 break-words">
                        {doc.translatedTitle || doc.title}
                      </h4>
                      <p className="text-gray-700 text-sm leading-relaxed line-clamp-4 break-words">
                        {doc.translatedSummary || doc.summary || doc.abstract}
                      </p>

                      <div className="flex flex-wrap gap-2 mt-4">
                        <button
                          onClick={(e) => {
                            e.stopPropagation();
                            openModal(doc.summary);
                          }}
                          title="Voir les détails"
                          className="flex items-center gap-1 px-3 py-1.5 rounded-full bg-blue-100 text-blue-700 hover:bg-blue-200 transition text-sm"
                        >
                          🔍 <span className="hidden sm:inline">Détails</span>
                        </button>

                        <button
                          onClick={(e) => {
                            e.stopPropagation();
                            handleDownloadPDF(doc);
                          }}
                          title="Télécharger le PDF"
                          className="flex items-center gap-1 px-3 py-1.5 rounded-full bg-purple-100 text-purple-700 hover:bg-purple-200 transition text-sm"
                        >
                          ⬇️{" "}
                          <span className="hidden sm:inline">Télécharger</span>
                        </button>
                        <button
                          onClick={(e) => {
                            e.stopPropagation();
                            handleAddToLibrary(doc);
                          }}
                          disabled={loadingFavorites[doc._id]}
                          className={`px-3 py-1.5 rounded text-sm ${
                            loadingFavorites[doc._id]
                              ? "bg-gray-400 cursor-not-allowed"
                              : "bg-yellow-400 text-yellow-900 hover:bg-yellow-500"
                          } transition`}
                        >
                          {loadingFavorites[doc._id]
                            ? "Ajout en cours..."
                            : "⭐ Ajouter aux favoris"}
                        </button>
                      </div>
                    </div>

                    {/* Partie droite : actions secondaires */}
                    <div className="flex flex-col gap-2 shrink-0 md:items-end w-full md:w-auto">
                      {doc.url && (
                        <a
                          href={doc.url}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="px-3 py-2 bg-gray-700 text-white rounded hover:bg-gray-800 transition flex items-center gap-2 text-sm w-full md:w-auto justify-center"
                          onClick={(e) => e.stopPropagation()}
                        >
                          📄 Lire PDF
                        </a>
                      )}

                      <button
                        onClick={(e) => {
                          e.stopPropagation();
                          handleQuestionDocument(doc);
                        }}
                        className="px-3 py-2 bg-yellow-100 text-yellow-800 rounded hover:bg-yellow-200 transition flex items-center gap-2 text-sm w-full md:w-auto justify-center"
                        title="Afficher les passages utilisés par l'IA et démarrer le chat"
                      >
                        💬 Discussion & Extraits
                      </button>

                  
                    </div>
                  </div>
                ))}
              </div>
            ))}

            {showModal && chatDoc && (
              <DocumentInsightModal
                isOpen={showModal}
                onClose={closeInsightModal}
                document={chatDoc}
                chunks={selectedChunks}
                fullText={selectedFullText}
              />
            )}
          </div>
        )}

        {/* Modal */}
        {modalContent && (
          <div
            onClick={closeModal}
            className="fixed inset-0 bg-black bg-opacity-60 flex justify-center items-center z-50"
          >
            <div
              role="dialog"
              tabIndex={-1}
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
    </>
  );
};

export default Search;
