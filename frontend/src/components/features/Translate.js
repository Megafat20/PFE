import React, { useState, useEffect, useRef } from "react";
import ProgressBar from "../ProgressBar";
const LANGUAGES = [
  { value: "haoussa", label: "Haoussa" },
  { value: "swahili", label: "Swahili" },
  { value: "zoulou", label: "Zoulou" },
  { value: "peul", label: "Peul" },
  { value: "malgache", label: "Malgache" },
  { value: "yoruba", label: "Yoruba" },
  { value: "amazigh", label: "Amazigh" },
  { value: "anglais", label: "Anglais" },
  { value: "français", label: "Français" },
  { value: "espagnol", label: "Espagnol" },
  { value: "allemand", label: "Allemand" },
  { value: "italien", label: "Italien" },
  { value: "portugais", label: "Portugais" },
  { value: "néerlandais", label: "Néerlandais" },
  { value: "chinois", label: "Chinois" },
  { value: "japonais", label: "Japonais" },
  { value: "coréen", label: "Coréen" },
  { value: "hindi", label: "Hindi" },
  { value: "arabe", label: "Arabe" },
  { value: "turc", label: "Turc" },
  { value: "russe", label: "Russe" },
  { value: "grec", label: "Grec" },
  { value: "latin", label: "Latin" },
];
const mainColor = "#1E40AF";
const Translate = () => {
  const [text, setText] = useState("");
  const [translateLang, setTranslateLang] = useState("swahili");
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [progress, setProgress] = useState(0);
  const [error, setError] = useState(null);
  const textareaRef = useRef(null);
  
  const token = localStorage.getItem("authToken");

  useEffect(() => {
    let interval;
    if (loading) {
      setProgress(0);
      interval = setInterval(() => {
        setProgress((oldProgress) => {
          if (oldProgress >= 90) return oldProgress; // stoppe à 90% (on attend la fin)
          return oldProgress + 10; // incrémente de 10% toutes les 200ms
        });
      }, 200);
    } else {
      setProgress(100); // passe à 100% quand c'est fini
      const timeout = setTimeout(() => setProgress(0), 300); // reset après un petit délai
      return () => clearTimeout(timeout);
    }
    return () => clearInterval(interval);
  }, [loading]);

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
    <>
    <ProgressBar progress={progress} />
    <div className="min-h-screen bg-[#F0F4FF] flex flex-col items-center justify-center p-6">
    {/* Feature Title Section */}
    <div className="max-w-7xl w-full mb-10 px-6 text-center">
      <h1
        className="text-4xl font-extrabold mb-2 tracking-tight select-none"
        style={{ color: mainColor }}
      >
        AI Translator Features
      </h1>
      <p
        className="text-lg font-medium select-none"
        style={{ color: mainColor }}
      >
        || Translate your Texte in any language.
      </p>
      <div
        className="mt-3 h-1 w-24 mx-auto rounded-full"
        style={{ backgroundColor: mainColor }}
      />
    </div>
  
    <div
      className="max-w-7xl w-full grid grid-cols-1 md:grid-cols-2 gap-10 bg-white rounded-xl shadow-lg border p-8"
      style={{ borderColor: "#93A6F9" }} // une teinte plus claire de bleu
    >
      {/* Colonne gauche */}
      <section className="flex flex-col">
      <h2 className="text-3xl font-extrabold text-blue-800 mb-6 select-none">
        🔁 Texte à traduire 
      </h2>
        <textarea
          id="input-text"
          ref={textareaRef}
          rows={8}
          placeholder="Tapez ou collez votre texte ici..."
          value={text}
          onChange={(e) => setText(e.target.value)}
          className="resize-none w-full rounded-lg border p-4
            text-gray-900 placeholder-[#1E40AF]/70 bg-[#E5E9FF]
            focus:outline-none transition duration-300 shadow-sm"
          style={{
            borderColor: mainColor,
            boxShadow: `0 0 8px 2px ${mainColor}44`,
          }}
          spellCheck="false"
        />
  
        <label
          htmlFor="language"
          className="mt-6 mb-2 font-semibold"
          style={{ color: "#3B59D1" }} // une nuance un peu plus claire
        >
          Choisissez la langue de traduction :
        </label>
        <select
          id="language"
          value={translateLang}
          onChange={(e) => setTranslateLang(e.target.value)}
          className="w-full rounded-lg border p-3
            bg-[#E5E9FF] text-[#1E40AF]
            focus:outline-none transition duration-300 shadow-sm"
          style={{ borderColor: "#93A6F9" }}
        >
          {LANGUAGES.map(({ value, label }) => (
            <option key={value} value={value}>
              {label}
            </option>
          ))}
        </select>
  
        <button
          onClick={handleTranslate}
          disabled={loading || !text.trim()}
          className="mt-8 text-white font-bold py-3 rounded-lg shadow-lg
            disabled:opacity-50 disabled:cursor-not-allowed transition-transform duration-300 hover:scale-[1.05]"
          style={{
            background: `linear-gradient(to right, ${mainColor}, #3B59D1)`,
          }}
        >
          {loading ? "Traduction en cours..." : "🌍 Traduire"}
        </button>
  
        {error && (
          <p className="mt-6 font-semibold select-text" style={{ color: "#D14343" }}>
            {error}
          </p>
        )}
      </section>
  
      {/* Colonne droite */}
      <section className="flex flex-col">
        <h2 className="text-2xl font-semibold mb-4 select-none" style={{ color: "#27408B" }}>
          Traduction
        </h2>
        <div
          className="flex-grow rounded-lg border p-6 text-gray-900 whitespace-pre-wrap overflow-auto min-h-[250px] shadow-inner"
          style={{ borderColor: "#93A6F9", backgroundColor: "#E5E9FF" }}
        >
          {loading ? (
            <p className="italic" style={{ color: "#27408B" }}>
              Génération en cours...
            </p>
          ) : result ? (
            <>{result}</>
          ) : (
            <p className="italic select-none" style={{ color: "#aac2ff" }}>
              La traduction apparaîtra ici.
            </p>
          )}
        </div>
      </section>
    </div>
  </div>
  </>
  );
};

export default Translate;
