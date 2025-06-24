import React, { useState } from "react";

const MAX_CHAR = 5000;

const Spinner = () => (
  <div className="flex justify-center my-8">
    <svg
      className="animate-spin h-10 w-10 text-blue-600 dark:text-indigo-400"
      xmlns="http://www.w3.org/2000/svg"
      fill="none"
      viewBox="0 0 24 24"
      aria-label="Chargement en cours"
    >
      <circle
        className="opacity-25"
        cx="12"
        cy="12"
        r="10"
        stroke="currentColor"
        strokeWidth="4"
      ></circle>
      <path
        className="opacity-75"
        fill="currentColor"
        d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z"
      ></path>
    </svg>
  </div>
);

const Alert = ({ type = "error", children }) => {
  const colors = {
    error: "bg-red-100 border-red-400 text-red-700",
    success: "bg-green-100 border-green-400 text-green-700",
    info: "bg-blue-100 border-blue-400 text-blue-700",
  };
  return (
    <div
      role="alert"
      className={`mb-6 p-4 border rounded-lg font-semibold text-center ${colors[type]}`}
    >
      {children}
    </div>
  );
};

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
      setFile(null);
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

  const textToList = (text) => {
    if (!text) return [];
    // On split par ligne, on nettoie les lignes vides
    return text
      .split('\n')
      .map(line => line.trim())
      .filter(line => line.length > 0);
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
    <div className="min-h-screen bg-gradient-to-r from-indigo-100 to-blue-100 dark:from-gray-900 dark:to-gray-800 flex items-center justify-center p-6">
      <div className="w-full max-w-5xl bg-white dark:bg-gray-900 rounded-3xl shadow-2xl p-10 space-y-8">
        <h1 className="text-5xl font-extrabold text-center bg-gradient-to-r from-indigo-600 to-blue-600 bg-clip-text text-transparent select-none mb-10">
          📄 Résumé Automatique
        </h1>

        {/* Textarea Input */}
        <div className="relative">
          <label htmlFor="text-input" className="block mb-2 font-semibold text-gray-700 dark:text-gray-300">
            Saisissez un texte (max {MAX_CHAR} caractères)
          </label>
          <textarea
            id="text-input"
            value={text}
            onChange={handleTextChange}
            rows={12}
            disabled={!!file}
            placeholder="Entrez votre texte ici..."
            className="w-full rounded-2xl border border-gray-300 dark:border-gray-700 p-6 resize-none shadow-sm focus:outline-none focus:ring-4 focus:ring-indigo-400 dark:focus:ring-blue-600 bg-gray-50 dark:bg-gray-800 text-gray-900 dark:text-gray-100 transition"
            aria-describedby="text-char-count"
          />
          <div
            id="text-char-count"
            className="absolute bottom-3 right-6 text-sm text-gray-500 dark:text-gray-400 select-none"
          >
            {text.length} / {MAX_CHAR}
          </div>
        </div>

        {/* File Upload */}
        <div>
          <label
            htmlFor="file-upload"
            className={`block mb-2 font-semibold cursor-pointer select-none ${
              text.trim() ? "opacity-50 cursor-not-allowed" : "hover:text-indigo-700 dark:hover:text-blue-400"
            } text-gray-700 dark:text-gray-300`}
          >
            📁 Télécharger un fichier PDF
          </label>
          <input
            id="file-upload"
            type="file"
            accept="application/pdf"
            onChange={handleFileChange}
            disabled={!!text.trim()}
            className="hidden"
          />
          <div
            className={`mt-2 p-4 rounded-xl border-2 border-dashed ${
              file ? "border-indigo-500 bg-indigo-50 dark:bg-indigo-900" : "border-gray-300 dark:border-gray-700"
            } text-center text-gray-600 dark:text-gray-400 select-none transition-colors duration-300`}
          >
            {file ? file.name : "Aucun fichier sélectionné"}
          </div>
        </div>

        {/* Action Buttons */}
        <div className="flex flex-col sm:flex-row gap-4 justify-center">
          {text.trim() && (
            <button
              onClick={sendTextForSummary}
              disabled={loading}
              className="flex-1 py-4 rounded-2xl bg-gradient-to-r from-indigo-600 to-blue-600 hover:from-indigo-700 hover:to-blue-700 text-white font-bold transition disabled:opacity-50 disabled:cursor-not-allowed shadow-lg"
            >
              📄 Résumer le texte saisi
            </button>
          )}
        </div>

        {/* Loading */}
        {loading && <Spinner />}

        {/* Error */}
        {error && <Alert type="error">❌ {error}</Alert>}

        {/* Result */}
        {result && (
  <div className="space-y-12">
    {/* Résumé simple */}
    <section
      className="p-6 rounded-2xl bg-gradient-to-tr from-indigo-50 to-blue-50 dark:from-gray-800 dark:to-gray-900 border border-indigo-200 dark:border-gray-700 shadow-lg"
    >
      <h2 className="text-3xl font-extrabold mb-4 text-indigo-700 dark:text-indigo-300 select-none">
        Résumé
      </h2>
      <p className="whitespace-pre-wrap text-gray-900 dark:text-gray-100 leading-relaxed tracking-wide">
        {result.summary}
      </p>
    </section>

    {/* Table des matières en liste */}
    <section
      className="p-6 rounded-2xl bg-gradient-to-tr from-indigo-50 to-blue-50 dark:from-gray-800 dark:to-gray-900 border border-indigo-200 dark:border-gray-700 shadow-xl transition-all"
    >
      <h2 className="text-3xl font-extrabold mb-5 text-blue-700 dark:text-blue-400 select-none">
        🧭 Table des matières
      </h2>
      <ul className="list-decimal list-inside space-y-2 text-gray-800 dark:text-gray-100 leading-relaxed tracking-wide max-h-96 overflow-y-auto scrollbar-thin scrollbar-thumb-blue-400 dark:scrollbar-thumb-blue-600 pr-2">
        {textToList(result.toc).map((item, idx) => (
          <li key={idx} className="hover:pl-2 transition-all duration-200">
            {item}
          </li>
        ))}
      </ul>
    </section>

    {/* Extraits pertinents en liste numérotée */}
    <section
      className="p-6 rounded-2xl bg-gradient-to-tr from-indigo-50 to-blue-50 dark:from-gray-800 dark:to-gray-900 border border-indigo-200 dark:border-gray-700 shadow-lg"
    >
      <h2 className="text-3xl font-extrabold mb-4 text-indigo-700 dark:text-indigo-300 select-none">
        🌟 Extraits pertinents
      </h2>
      <ol className="list-decimal list-inside space-y-2 text-gray-900 dark:text-gray-100 leading-relaxed tracking-wide max-h-96 overflow-y-auto">
        {textToList(result.highlights).map((item, idx) => (
          <li key={idx}>{item}</li>
        ))}
      </ol>
    </section>
  </div>
)}

      </div>
    </div>
  );
};

export default AutoSummary;
