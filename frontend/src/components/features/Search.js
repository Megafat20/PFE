import React, { useEffect, useState } from "react";
import { io } from "socket.io-client";


const token = localStorage.getItem("authToken");
const socket = io("http://localhost:5000", {
  autoConnect: false,
  transports: ["websocket"],
  auth: {
    token: token
  }
});

const Search = () => {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);



  useEffect(() => {
    // Connect socket only once
    socket.connect();

    socket.on("connect", () => {
      console.log("Socket connecté avec ID :", socket.id);
    });
  
    socket.on("disconnect", () => {
      console.log("Socket déconnecté");
    });
    socket.on("document", (doc) => {
      console.log("Document reçu :", doc);
      setResults((prev) => [...prev, doc]);
    });

    socket.on("search_complete", () => {
      setLoading(false);
    });

    socket.on("search_error", (msg) => {
      setError(msg);
      setLoading(false);
    });

    socket.on("connect_error", (err) => {
      setError("Erreur de connexion au serveur");
      setLoading(false);
    });

    return () => {
      socket.off("document");
      socket.off("search_complete");
      socket.off("search_error");
      socket.off("connect_error");
      // socket.disconnect();
    };
  }, []);

  const handleSearch = () => {
    if (!query.trim()) return;

    setLoading(true);
    setError(null);
    setResults([]);

    socket.emit("start_document_search", {
      query: query.trim(),
      max_results: 10,
      token: token,
    });
  };

  return (
    <div className="max-w-4xl mx-auto p-4">
      <h1 className="text-3xl font-extrabold mb-6 text-blue-700 select-none">
        🔍 Recherche de documents scientifiques
      </h1>

      <div className="flex gap-3 mb-6">
        <input
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && handleSearch()}
          className="flex-grow border border-gray-300 rounded-lg px-4 py-3 text-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
          placeholder="Ex : Large Language Models for Scientific Research"
          disabled={loading}
        />
        <button
          onClick={handleSearch}
          disabled={loading || !query.trim()}
          className="bg-blue-600 disabled:bg-blue-300 text-white font-semibold px-6 py-3 rounded-lg hover:bg-blue-700 transition"
        >
          {loading ? "Recherche..." : "Rechercher"}
        </button>
      </div>

      {error && (
        <p className="text-red-600 mb-6 font-medium select-none">{error}</p>
      )}

      {!loading && results.length === 0 && query.trim() !== "" && !error && (
        <p className="text-gray-600 mb-6 select-none">Aucun document trouvé.</p>
      )}

      <ul className="space-y-6">
        {results.map((doc, index) => (
          <li
            key={index}
            className="border border-gray-200 p-6 rounded-lg shadow hover:shadow-lg transition"
          >
            <h2 className="font-bold text-xl mb-1">{doc.title}</h2>
            <p className="text-sm text-gray-500 italic mb-2">
              {doc.authors?.join(", ")}
            </p>
            <p className="text-gray-700 mb-3">{doc.summary}</p>
            <a
              href={doc.pdf_url}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-block text-blue-600 hover:text-blue-800 font-semibold"
            >
              Voir le PDF 📄
            </a>
          </li>
        ))}
      </ul>
    </div>
  );
};

export default Search;
