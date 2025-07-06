import React, { useState, useRef, useEffect } from "react";
import { motion } from "framer-motion";

const ReportGenerator = () => {
  const [objective, setObjective] = useState("");
  const [language, setLanguage] = useState("fr");
  const [style, setStyle] = useState("scientifique");
  const [report, setReport] = useState("");
  const [loading, setLoading] = useState(false);
  const [copied, setCopied] = useState(false);
  const reportRef = useRef(null);

  const handleGenerate = async () => {
    if (!objective.trim()) {
      alert("Merci de renseigner l'objectif.");
      return;
    }

    setReport("");
    setLoading(true);

    try {
      const token = localStorage.getItem("authToken");
      const res = await fetch("http://localhost:5000/multi/generate_report", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          objective,
          language,
          style,
        }),
      });

      if (!res.body) throw new Error("Pas de flux reçu.");

      const reader = res.body.getReader();
      const decoder = new TextDecoder("utf-8");

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        const chunk = decoder.decode(value);
        setReport((prev) => prev + chunk);
      }
    } catch (err) {
      alert("Erreur lors de la génération du rapport : " + err.message);
    }

    setLoading(false);
  };

  useEffect(() => {
    if (reportRef.current) {
      reportRef.current.scrollTop = reportRef.current.scrollHeight;
    }
  }, [report]);

  const handleCopy = () => {
    navigator.clipboard.writeText(report);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handleDownloadPDF = async () => {
    if (!report.trim()) return;
    const token = localStorage.getItem("authToken");
    const res = await fetch("http://localhost:5000/multi/generate_report_pdf", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${token}`,
      },
      body: JSON.stringify({
        content: report,
        filename: "rapport_final",
      }),
    });

    if (!res.ok) {
      alert("Erreur lors du téléchargement du PDF");
      return;
    }

    const blob = await res.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "rapport.pdf";
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(url);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-indigo-100 to-blue-100 dark:from-gray-900 dark:to-gray-950 py-10 px-4 md:px-12">
      <div className="max-w-7xl mx-auto space-y-8">
        {/* Header */}
        <h1 className="text-4xl font-extrabold text-indigo-700 dark:text-white mb-4 text-center">
          📘 Générateur de Rapport Intelligent
        </h1>

        {/* Main Layout */}
        <div className="flex flex-col md:flex-row gap-8">
          {/* Form */}
          <div className="w-full md:w-1/2 bg-white dark:bg-gray-800 border dark:border-gray-700 rounded-2xl p-6 shadow space-y-5">
            <h2 className="text-2xl font-semibold text-indigo-600 dark:text-indigo-300">
              ⚙️ Configuration
            </h2>

            <textarea
              placeholder="Décris l'objectif du rapport..."
              value={objective}
              onChange={(e) => setObjective(e.target.value)}
              rows={4}
              disabled={loading}
              className="w-full p-3 rounded-xl border dark:border-gray-600 dark:bg-gray-900 dark:text-white"
            />

            <div className="flex flex-col gap-3">
              <select
                value={language}
                onChange={(e) => setLanguage(e.target.value)}
                disabled={loading}
                className="p-3 rounded-xl border dark:border-gray-600 dark:bg-gray-900 dark:text-white"
              >
                <option value="fr">🇫🇷 Français</option>
                <option value="en">🇬🇧 Anglais</option>
              </select>

              <select
                value={style}
                onChange={(e) => setStyle(e.target.value)}
                disabled={loading}
                className="p-3 rounded-xl border dark:border-gray-600 dark:bg-gray-900 dark:text-white"
              >
                <option value="scientifique">🔬 Scientifique</option>
                <option value="technique">⚙️ Technique</option>
                <option value="commercial">📈 Commercial</option>
              </select>
            </div>

            <button
              onClick={handleGenerate}
              disabled={loading}
              className={`w-full py-3 rounded-xl font-semibold transition ${
                loading
                  ? "bg-gray-400 cursor-not-allowed text-white"
                  : "bg-indigo-600 hover:bg-indigo-700 text-white"
              }`}
            >
              {loading ? "Génération en cours..." : "✨ Générer le rapport"}
            </button>
          </div>

          {/* Report Preview */}
          <div className="w-full md:w-1/2 bg-white dark:bg-gray-800 border dark:border-gray-700 rounded-2xl p-6 shadow relative flex flex-col max-h-[600px]">
            <div className="flex items-center justify-between mb-3">
              <h3 className="text-xl font-semibold text-gray-800 dark:text-white">📄 Rapport généré</h3>
              <div className="flex gap-2">
                <button
                  onClick={handleCopy}
                  disabled={!report}
                  className="text-sm px-3 py-1 rounded-full bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 transition"
                >
                  {copied ? "✅ Copié !" : "📋 Copier"}
                </button>
                <button
                  onClick={handleDownloadPDF}
                  disabled={!report}
                  className="text-sm px-3 py-1 rounded-full bg-green-100 text-green-700 hover:bg-green-200 transition"
                >
                  📥 PDF
                </button>
              </div>
            </div>

            <motion.pre
              ref={reportRef}
              className="flex-grow overflow-auto whitespace-pre-wrap text-sm text-gray-700 dark:text-gray-100 scrollbar-thin scrollbar-thumb-indigo-400 scrollbar-track-gray-200 hover:scrollbar-thumb-indigo-600"
              initial={{ opacity: 0.4 }}
              animate={{ opacity: 1 }}
              transition={{ duration: 0.4 }}
            >
              {report || "⏳ En attente d'une génération de rapport..."}
            </motion.pre>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ReportGenerator;
