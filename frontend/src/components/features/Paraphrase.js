import React, { useState, useEffect, useRef } from "react";
import ProgressBar from "../ProgressBar";
const PARAPHRASE_STYLES = [
  { value: "fluent", label: "Fluent" },
  { value: "formal", label: "Formal" },
  { value: "creative", label: "Creative" },
  { value: "academic", label: "Academic" },
  { value: "professional", label: "Professional" },
];

const PRIMARY_COLOR = "#1E40AF"; // bleu foncé
const SECONDARY_COLOR = "#3B82F6"; // bleu moyen
const LIGHT_BG_COLOR = "#DBEAFE"; // bleu très clair
const BORDER_COLOR = "#93C5FD"; // bleu clair pour bordures

const Paraphrase = () => {
  const [text, setText] = useState("");
  const [style, setStyle] = useState("fluent");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const inputRef = useRef(null);

  const token = localStorage.getItem("authToken");

  useEffect(() => {
    if (inputRef.current) {
      inputRef.current.style.height = "auto";
      inputRef.current.style.height = inputRef.current.scrollHeight + "px";
    }
  }, [text]);

  const handleParaphrase = async () => {
    if (!text.trim()) return alert("Veuillez saisir un texte à paraphraser.");
    setLoading(true);
    setError(null);
    setResult(null);

    try {
      const res = await fetch("http://localhost:5000//multi/paraphrase", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ texte: text, style }),
      });

      const data = await res.json();

      if (res.ok) setResult(data.paraphrase || "Aucune paraphrase reçue.");
      else setError(data.error || "Erreur lors de la paraphrase.");
    } catch (err) {
      setError(err.message);
    }

    setLoading(false);
  };

  return (
      <>
        <ProgressBar loading={loading} />
        <div className="min-h-screen bg-gradient-to-br from-blue-50 to-white flex flex-col items-center justify-center p-6">

      {/* Titre */}
      <div className="max-w-7xl w-full mb-10 px-6 text-center">
        <h1 className="text-4xl font-extrabold text-blue-800 mb-2 tracking-tight select-none">
          AI Paraphrases Features
        </h1>
        <p className="text-lg text-blue-600 font-medium select-none">
          || Free Make your academic writing clear and original.
        </p>
        <div className="mt-3 h-1 w-24 bg-blue-400 mx-auto rounded-full"></div>
      </div>

      {/* Carte principale */}
      <div className="max-w-7xl w-full bg-white rounded-xl shadow-xl border border-blue-300 p-8 grid grid-cols-1 md:grid-cols-2 gap-8">
        {/* Input */}
        <div className="flex flex-col">
          <h2 className="text-3xl font-extrabold text-blue-800 mb-6 select-none">
            🔁 Paraphraser un texte
          </h2>
          <textarea
            ref={inputRef}
            rows={6}
            placeholder="Tapez ou collez votre texte ici..."
            value={text}
            onChange={(e) => setText(e.target.value)}
            className="resize-none w-full rounded-lg border border-blue-300 p-4
                      text-gray-900 placeholder-blue-400
                      focus:outline-none focus:ring-4 focus:ring-blue-400
                      transition duration-300 shadow-sm"
            spellCheck="false"
          />

          <label htmlFor="style" className="mt-6 mb-2 font-semibold text-blue-700">
            Choisissez le style de paraphrase :
          </label>
          <select
            id="style"
            value={style}
            onChange={(e) => setStyle(e.target.value)}
            className="w-full rounded-lg border border-blue-300 p-3
                      focus:outline-none focus:ring-4 focus:ring-blue-400
                      transition duration-300 shadow-sm"
          >
            {PARAPHRASE_STYLES.map(({ value, label }) => (
              <option key={value} value={value}>
                {label}
              </option>
            ))}
          </select>

          <button
            onClick={handleParaphrase}
            disabled={loading || !text.trim()}
            className="mt-6 bg-gradient-to-r from-blue-600 to-blue-700
                      hover:from-blue-700 hover:to-blue-800
                      text-white font-bold py-3 rounded-lg
                      shadow-lg disabled:opacity-50 disabled:cursor-not-allowed
                      transition-transform duration-300 hover:scale-[1.05]"
          >
            {loading ? "Paraphrase en cours..." : "Paraphraser"}
          </button>

          {error && <p className="mt-4 text-red-600 font-semibold">{error}</p>}
        </div>

        {/* Output */}
        <div className="flex flex-col">
          <h3 className="text-2xl font-semibold text-blue-800 mb-4 select-none">
            Résultat
          </h3>
          <div
            className="flex-grow rounded-lg border border-blue-300 bg-blue-50
                      p-6 text-gray-900 whitespace-pre-wrap
                      overflow-auto min-h-[200px] shadow-inner"
          >
            {loading ? (
              <p className="italic text-blue-600">Génération en cours...</p>
            ) : result ? (
              <>{result}</>
            ) : (
              <p className="text-blue-400 italic">Le texte reformulé apparaîtra ici.</p>
            )}
          </div>
        </div>
      </div>
    </div>
    </>
  );
  
};

export default Paraphrase;
