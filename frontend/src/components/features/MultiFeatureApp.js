import React, { useState } from "react";

const MultiFeatureApp = () => {
  const [selectedFeature, setSelectedFeature] = useState("");
  const [text, setText] = useState("");
  const [file, setFile] = useState(null);
  const [docTitle, setDocTitle] = useState("");
  const [translateLang, setTranslateLang] = useState("swahili");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const token = localStorage.getItem("authToken");

  const handleFileChange = (e) => {
    setFile(e.target.files[0]);
    setText("");
    setResult(null);
  };

  const callSummarize = async () => {
    setLoading(true);
    setError(null);
    setResult(null);

    try {
      if (file) {
        const formData = new FormData();
        formData.append("file", file);

        const res = await fetch("/multi/summarize_pdf", {
          method: "POST",
          headers: {
            Authorization: `Bearer ${token}`,
          },
          body: formData,
        });

        const data = await res.json();
        if (res.ok) setResult(data.summary);
        else setError(data.error || "Erreur lors du résumé PDF");
      } else if (text.trim()) {
        const res = await fetch("http://localhost:5000/multi/summarize", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            Authorization: `Bearer ${token}`,
          },
          body: JSON.stringify({ documents: text }),
        });

        const data = await res.json();
        if (res.ok) setResult(data.summary);
        else setError(data.error || "Erreur lors du résumé texte");
      } else {
        alert("Veuillez entrer un texte ou choisir un fichier PDF");
      }
    } catch (e) {
      setError(e.message);
    }

    setLoading(false);
  };

  const handleTranslate = async () => {
    if (!text.trim()) return alert("Veuillez saisir du texte à traduire.");
    setLoading(true);
    setError(null);
    setResult(null);

    try {
      const res = await fetch("http://localhost:5000/multi/translate", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ text, language: translateLang }),
      });

      const data = await res.json();
      if (res.ok) setResult(data.translation || "Aucune traduction reçue.");
      else setError(data.error || "Erreur lors de la traduction");
    } catch (err) {
      setError(err.message);
    }

    setLoading(false);
  };

  const handleParaphrase = async () => {
    if (!text.trim()) return alert("Veuillez saisir un texte.");
    setLoading(true); setError(null); setResult(null);
    try {
      const res = await fetch("/multi/paraphrase", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ texte: text }),
      });
      const data = await res.json();
      if (res.ok) setResult(data.paraphrase);
      else setError(data.error || "Erreur lors de la paraphrase");
    } catch (err) {
      setError(err.message);
    }
    setLoading(false);
  };
  


  return (
    <div className="min-h-screen pt-16 bg-gradient-to-r from-slate-100 via-white to-slate-200 py-10 px-4">
      <div className="max-w-3xl mx-auto bg-white p-8 rounded-lg shadow-md border border-gray-200">
        <h1 className="text-3xl font-bold text-center text-blue-700 mb-8">
          🤖 Assistant IA Multi-Fonctions
        </h1>

          <select
            className="border rounded p-2 mb-4 w-full"
            value={selectedFeature}
            onChange={(e) => {
              setSelectedFeature(e.target.value);
              setResult(null);
              setText("");
              setFile(null);
              setDocTitle("");
              setError(null);
            }}
          >
            <option value="">-- Sélectionnez une fonctionnalité --</option>
            <option value="summary">📄 Résumé de texte ou PDF</option>
            <option value="translate">🌍 Traduction en langues africaines</option>
            <option value="generate">📝 Générateur de documents (CV, lettre...)</option>
            <option value="paraphrase">🔁 Paraphraser un texte</option>
            <option value="rapport">📑 Générer un rapport</option>
            <option value="intro">📘 Introduction académique</option>
            <option value="extract">🔍 Extraction d'infos clés</option>
            <option value="grammar">✍️ Correction grammaticale</option>
            <option value="trends">📊 Analyse de tendances scientifiques</option>
          </select>


        {selectedFeature === "summary" && (
          <div>
            <textarea
              rows={6}
              placeholder="Tapez ou collez le texte à résumer ici..."
              value={text}
              onChange={(e) => setText(e.target.value)}
              className="w-full border border-gray-300 rounded-md p-3 mb-4 focus:ring focus:ring-blue-200"
              disabled={!!file}
            />
            <input
              type="file"
              accept="application/pdf"
              onChange={handleFileChange}
              className="block mb-4"
              disabled={!!text.trim()}
            />
            <button
              onClick={callSummarize}
              disabled={loading || (!text.trim() && !file)}
              className="bg-blue-600 hover:bg-blue-700 text-white px-5 py-2 rounded-md transition"
            >
              {loading ? "Génération..." : "📄 Générer le résumé"}
            </button>
          </div>
        )}

        {selectedFeature === "translate" && (
          <div>
            <textarea
              rows={6}
              placeholder="Tapez ou collez le texte à traduire ici..."
              value={text}
              onChange={(e) => setText(e.target.value)}
              className="w-full border border-gray-300 rounded-md p-3 mb-4 focus:ring focus:ring-blue-200"
            />
            <select
              className="w-full border border-gray-300 rounded-md p-3 mb-4"
              value={translateLang}
              onChange={(e) => setTranslateLang(e.target.value)}
            >
                <option value="haoussa">Haoussa</option>
                <option value="swahili">Swahili</option>
                <option value="zoulou">Zoulou</option>
                <option value="peul">Peul</option>
                <option value="malgache">Malgache</option>
                <option value="yoruba">Yoruba</option>
                <option value="amazigh">Amazigh</option>

                {/* Langues européennes */}
                <option value="anglais">Anglais</option>
                <option value="français">Français</option>
                <option value="espagnol">Espagnol</option>
                <option value="allemand">Allemand</option>
                <option value="italien">Italien</option>
                <option value="portugais">Portugais</option>
                <option value="néerlandais">Néerlandais</option>

                {/* Langues asiatiques */}
                <option value="chinois">Chinois</option>
                <option value="japonais">Japonais</option>
                <option value="coréen">Coréen</option>
                <option value="hindi">Hindi</option>
                <option value="arabe">Arabe</option>
                <option value="turc">Turc</option>

                {/* Autres */}
                <option value="russe">Russe</option>
                <option value="grec">Grec</option>
                <option value="latin">Latin</option>
            </select>
            <button
              onClick={handleTranslate}
              disabled={loading || !text.trim()}
              className="bg-green-600 hover:bg-green-700 text-white px-5 py-2 rounded-md transition"
            >
              {loading ? "Traduction..." : "🌍 Traduire"}
            </button>
          </div>
        )}

        {selectedFeature === "paraphrase" && (
          <>
            <textarea
              rows={6}
              placeholder="Texte à reformuler"
              value={text}
              onChange={(e) => setText(e.target.value)}
              className="border rounded p-2 w-full mb-2"
            />
            <button className="btn-primary" onClick={handleParaphrase} disabled={loading}>
              {loading ? "Paraphrase en cours..." : "Paraphraser"}
            </button>
          </>
        )}


        {result && (
          <div className="mt-6 p-4 bg-green-50 border border-green-300 text-gray-800 rounded-md whitespace-pre-wrap">
            <strong>Résultat :</strong>
            <p>{result}</p>
          </div>
        )}
      </div>
    </div>
  );
};

export default MultiFeatureApp;
