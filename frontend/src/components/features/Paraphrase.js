import React, { useState } from "react";

const Paraphrase = () => {
  const [text, setText] = useState("");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const token = localStorage.getItem("authToken");

  const handleParaphrase = async () => {
    if (!text.trim()) return alert("Veuillez saisir un texte.");
    setLoading(true);
    setError(null);
    setResult(null);
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
    <div className="max-w-3xl mx-auto p-8 bg-white rounded shadow-md border border-gray-200">
      <h2 className="text-2xl font-bold mb-6 text-purple-700">🔁 Paraphraser un texte</h2>
      <textarea
        rows={6}
        placeholder="Texte à reformuler"
        value={text}
        onChange={(e) => setText(e.target.value)}
        className="border rounded p-3 w-full mb-4"
      />
      <button
        onClick={handleParaphrase}
        disabled={loading}
        className="bg-purple-600 hover:bg-purple-700 text-white px-5 py-2 rounded-md transition"
      >
        {loading ? "Paraphrase en cours..." : "Paraphraser"}
      </button>

      {error && <p className="mt-4 text-red-600">{error}</p>}

      {result && (
        <div className="mt-6 p-4 bg-purple-50 border border-purple-300 text-gray-800 rounded-md whitespace-pre-wrap">
          <strong>Paraphrase :</strong>
          <p>{result}</p>
        </div>
      )}
    </div>
  );
};

export default Paraphrase;
