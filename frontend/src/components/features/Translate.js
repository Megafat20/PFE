import React, { useState } from "react";

const Translate = () => {
  const [text, setText] = useState("");
  const [translateLang, setTranslateLang] = useState("swahili");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const token = localStorage.getItem("authToken");

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

  return (
    <div className="max-w-3xl mx-auto p-8 bg-white rounded shadow-md border border-gray-200">
      <h2 className="text-2xl font-bold mb-6 text-green-700">🌍 Traduction en langues africaines</h2>
      <textarea
        rows={6}
        placeholder="Tapez ou collez le texte à traduire ici..."
        value={text}
        onChange={(e) => setText(e.target.value)}
        className="w-full border border-gray-300 rounded-md p-3 mb-4 focus:ring focus:ring-green-200"
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

      {error && <p className="mt-4 text-red-600">{error}</p>}

      {result && (
        <div className="mt-6 p-4 bg-green-50 border border-green-300 text-gray-800 rounded-md whitespace-pre-wrap">
          <strong>Traduction :</strong>
          <p>{result}</p>
        </div>
      )}
    </div>
  );
};

export default Translate;
