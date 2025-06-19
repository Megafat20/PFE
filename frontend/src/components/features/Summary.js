import React, { useState } from "react";

const MAX_CHAR = 5000;

const AutoSummary = () => {
  const [text, setText] = useState("");
  const [file, setFile] = useState(null);
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const token = localStorage.getItem("authToken");

  const handleFileChange = async (e) => {
    const selectedFile = e.target.files[0];
    if (!selectedFile || selectedFile.type !== "application/pdf") {
      setError("Veuillez sélectionner un fichier PDF valide.");
      return;
    }

    setFile(selectedFile);
    setText("");
    setResult(null);
    setError(null);
    await sendFileForAnalysis(selectedFile);
  };

  const handleTextChange = (e) => {
    if (e.target.value.length <= MAX_CHAR) {
      setText(e.target.value);
      setFile(null);
      setResult(null);
      setError(null);
    }
  };

  const sendFileForAnalysis = async (file) => {
    setLoading(true);
    try {
      const formData = new FormData();
      formData.append("file", file);

      const response = await fetch("http://localhost:5000/multi/upload_summarize", {
        method: "POST",
        headers: {
          Authorization: `Bearer ${token}`,
        },
        body: formData,
      });

      const data = await response.json();
      if (!response.ok) throw new Error(data.error || "Erreur serveur");

      setResult(data);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  const sendTextForSummary = async () => {
    if (!text.trim()) {
      setError("Veuillez saisir un texte à résumer.");
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);

    try {
      const response = await fetch("http://localhost:5000/multi/summary", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ documents: text }),
      });

      const data = await response.json();
      if (!response.ok) throw new Error(data.error || "Erreur serveur");

      setResult(data);
    } catch (e) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-gray-50 dark:bg-gray-900 py-10 px-4 flex items-center justify-center">
      <div className="max-w-4xl w-full bg-white dark:bg-gray-800 rounded-2xl shadow-lg p-8">
        <h1 className="text-3xl font-bold mb-6 text-center text-blue-700 dark:text-blue-300">📄 Résumé Automatique</h1>

        <textarea
          value={text}
          onChange={handleTextChange}
          placeholder="Saisissez un texte ici..."
          rows={8}
          disabled={!!file}
          className="w-full resize-none mb-3 p-4 rounded-xl border border-gray-300 dark:border-gray-700 focus:ring-2 focus:ring-blue-400 disabled:bg-gray-200 dark:disabled:bg-gray-800"
        />
        <div className="text-right mb-4 text-sm text-gray-500 dark:text-gray-400">
          {text.length} / {MAX_CHAR} caractères
        </div>

        <input
          type="file"
          accept="application/pdf"
          onChange={handleFileChange}
          className="mb-6 block w-full text-sm file:mr-4 file:py-2 file:px-4 file:rounded-xl file:border-0 file:font-semibold file:bg-blue-600 file:text-white hover:file:bg-blue-700"
          disabled={!!text.trim()}
        />

        {text.trim() && (
          <button
            onClick={sendTextForSummary}
            disabled={loading}
            className="w-full mb-6 py-3 rounded-xl text-white font-bold bg-blue-600 hover:bg-blue-700 disabled:bg-blue-300"
          >
            📄 Résumer le texte saisi
          </button>
        )}

        {loading && <p className="text-center text-blue-500 animate-pulse">⏳ Traitement en cours...</p>}

        {error && (
          <div className="mt-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded-lg">
            ❌ {error}
          </div>
        )}

        {result && (
          <div className="mt-6 space-y-6">
            <div>
              <h2 className="text-xl font-semibold text-green-700 dark:text-green-400 mb-2">Résumé</h2>
              <p className="whitespace-pre-wrap">{result.summary}</p>
            </div>
            <div>
              <h2 className="text-xl font-semibold text-blue-700 dark:text-blue-400 mb-2">🧭 Table des matières</h2>
              <p className="whitespace-pre-wrap">{result.toc}</p>
            </div>
            <div>
              <h2 className="text-xl font-semibold text-indigo-700 dark:text-indigo-300 mb-2">🌟 Extraits pertinents</h2>
              <p className="whitespace-pre-wrap">{result.highlights}</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};

export default AutoSummary;
