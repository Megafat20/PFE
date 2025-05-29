import React, { useState } from "react";

const Summary = () => {
  const [text, setText] = useState("");
  const [file, setFile] = useState(null);
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
        const res = await fetch("/multi/summarize", {
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

  return (
    <div className="max-w-3xl mx-auto p-8 bg-white rounded shadow-md border border-gray-200">
      <h2 className="text-2xl font-bold mb-6 text-blue-700">📄 Résumé de texte ou PDF</h2>
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

      {error && <p className="mt-4 text-red-600">{error}</p>}

      {result && (
        <div className="mt-6 p-4 bg-green-50 border border-green-300 text-gray-800 rounded-md whitespace-pre-wrap">
          <strong>Résumé :</strong>
          <p>{result}</p>
        </div>
      )}
    </div>
  );
};

export default Summary;
