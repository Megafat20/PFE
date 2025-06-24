import React, { useEffect, useState } from "react";
import { BookOpen, Globe, Loader2 } from "lucide-react";

export default function FavoritesList() {
  const [favorites, setFavorites] = useState([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    const fetchFavorites = async () => {
      try {
        const token = localStorage.getItem("authToken");
        const res = await fetch("http://localhost:5000/favorites", {
          method: "GET",
          headers: {
            Authorization: `Bearer ${token}`,
          },
        });
        if (!res.ok) {
          const err = await res.json();
          throw new Error(err.error || "Erreur lors du chargement des favoris");
        }
        const data = await res.json();
        setFavorites(data.favorites || []);
      } catch (err) {
        setError(err.message);
      } finally {
        setLoading(false);
      }
    };

    fetchFavorites();
  }, []);

  if (loading) {
    return (
      <div className="flex items-center justify-center mt-10 text-blue-600">
        <Loader2 className="animate-spin w-6 h-6 mr-2" />
        Chargement des favoris...
      </div>
    );
  }

  if (error) {
    return (
      <div className="text-red-600 mt-10 text-center">
        ❌ Erreur : {error}
      </div>
    );
  }

  if (favorites.length === 0) {
    return (
      <div className="text-gray-500 mt-10 text-center">
        📭 Aucun document favori pour le moment.
      </div>
    );
  }

  return (
    <div className="p-6 max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-blue-800 mb-6">⭐ Mes Favoris</h1>

      <div className="grid grid-cols-1 sm:grid-cols-2 gap-6">
        {favorites.map((doc) => (
          <div
            key={doc._id}
            onClick={() => window.open(doc.url || "#", "_blank")}
            className="p-5 bg-white rounded-2xl border shadow hover:shadow-md transition cursor-pointer"
            title="Cliquez pour ouvrir"
          >
            <h2 className="text-xl font-semibold text-gray-900 mb-2 line-clamp-2">
              <BookOpen className="inline mr-2 text-blue-500" />
              {doc.title || doc.filename || "Document sans titre"}
            </h2>
            <p className="text-sm text-gray-600 line-clamp-3 mb-3">
              {doc.summary || doc.abstract || "Aucun résumé disponible."}
            </p>
            <span className="inline-flex items-center text-xs bg-blue-100 text-blue-800 px-3 py-1 rounded-full">
              <Globe className="mr-1 h-4 w-4" />
              {doc.source || "Source inconnue"}
            </span>
          </div>
        ))}
      </div>
    </div>
  );
}
