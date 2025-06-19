import React, { useEffect, useState } from "react";

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

  if (loading) return <p>Chargement des favoris...</p>;
  if (error) return <p className="text-red-600">Erreur : {error}</p>;
  if (favorites.length === 0) return <p>Aucun favori pour le moment.</p>;

  return (
    <div className="p-4">
      <h1 className="text-2xl font-bold mb-4">Mes Favoris</h1>
      <ul className="space-y-3">
        {favorites.map((doc) => (
          <li
            key={doc._id}
            className="p-3 border rounded cursor-pointer hover:bg-gray-100 transition"
            onClick={() => window.open(doc.url || "#", "_blank")}
            title="Ouvrir le document"
          >
            <h2 className="font-semibold text-lg">{doc.title || doc.filename || "Document sans titre"}</h2>
            <p className="text-sm text-gray-600">{doc.source || "Source inconnue"}</p>
          </li>
        ))}
      </ul>
    </div>
  );
}
