import React, { useState, useEffect } from "react";
import { FiSearch, FiFolder, FiPlus, FiX, FiTrash2 } from "react-icons/fi";

function Sidebar({
  onSelectConversation,
  selectedConversationId,
  onSelectDocument,
  chatSearchQuery,
  setChatSearchQuery,
  handleChatSearch,
}) {
  const [history, setHistory] = useState([]);
  const [documents, setDocuments] = useState([]);
  const [showSearchModal, setShowSearchModal] = useState(false);
  const [showLibraryModal, setShowLibraryModal] = useState(false);
  const [conversationSearchQuery, setConversationSearchQuery] = useState("");

  useEffect(() => {
    const fetchHistory = async () => {
      const token = localStorage.getItem("authToken");
      try {
        const res = await fetch("http://localhost:5000/conversations", {
          headers: { Authorization: `Bearer ${token}` },
        });
        if (res.ok) setHistory(await res.json());
      } catch (err) {
        console.error("Erreur réseau :", err);
      }
    };
    fetchHistory();
  }, []);

  const fetchDocuments = async () => {
    const token = localStorage.getItem("authToken");
    try {
      const res = await fetch("http://localhost:5000/documents", {
        headers: { Authorization: `Bearer ${token}` },
      });
      if (res.ok) setDocuments(await res.json());
    } catch (err) {
      console.error("Erreur documents :", err);
    }
  };

  const handleDeleteConversation = async (convId) => {
    if (!window.confirm("Voulez-vous vraiment supprimer cette conversation ?"))
      return;

    const token = localStorage.getItem("authToken");
    try {
      const res = await fetch(`http://localhost:5000/conversations/${convId}`, {
        method: "DELETE",
        headers: { Authorization: `Bearer ${token}` },
      });

      if (res.ok) {
        setHistory((prev) => prev.filter((c) => c.id !== convId));
        if (convId === selectedConversationId) {
          onSelectConversation(null);
        }
      } else {
        const errData = await res.json();
        alert(
          `Erreur : ${
            errData.error || "Impossible de supprimer la conversation."
          }`
        );
      }
    } catch (err) {
      alert("Erreur réseau lors de la suppression");
      console.error(err);
    }
  };

  const handleNewConversation = async () => {
    const token = localStorage.getItem("authToken");
    try {
      const res = await fetch("http://localhost:5000/conversations", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({ title: "Nouvelle conversation" }),
      });
      if (res.ok) {
        const data = await res.json();
        onSelectConversation(data.conversation_id);
        const updated = await fetch("http://localhost:5000/conversations", {
          headers: { Authorization: `Bearer ${token}` },
        });
        setHistory(await updated.json());
      }
    } catch (err) {
      console.error("Erreur création conversation :", err);
    }
  };

  const openLibraryModal = async () => {
    await fetchDocuments();
    setShowLibraryModal(true);
  };

  return (
    <aside className="w-80 bg-gray-900 text-white flex-shrink-0 border-r border-gray-800 p-5 space-y-6 overflow-y-auto">
      {/* Actions */}
      <div className="space-y-3">
        <button
          onClick={() => setShowSearchModal(true)}
          className="w-full flex items-center gap-3 bg-blue-600 hover:bg-blue-700 py-3 px-4 rounded-lg text-sm font-semibold shadow-md transition duration-300 ease-in-out focus:outline-none focus:ring-2 focus:ring-blue-400"
          aria-label="Rechercher dans les conversations"
        >
          <FiSearch className="text-lg" /> Rechercher dans les conversations
        </button>
        <button
          onClick={handleNewConversation}
          className="w-full flex items-center gap-3 bg-green-600 hover:bg-green-700 py-3 px-4 rounded-lg text-sm font-semibold shadow-md transition duration-300 ease-in-out focus:outline-none focus:ring-2 focus:ring-green-400"
          aria-label="Nouvelle conversation"
        >
          <FiPlus className="text-lg" /> Nouvelle conversation
        </button>
        <button
          onClick={openLibraryModal}
          className="w-full flex items-center gap-3 bg-indigo-600 hover:bg-indigo-700 py-3 px-4 rounded-lg text-sm font-semibold shadow-md transition duration-300 ease-in-out focus:outline-none focus:ring-2 focus:ring-indigo-400"
          aria-label="Bibliothèque"
        >
          <FiFolder className="text-lg" /> Bibliothèque
        </button>
      </div>

      <hr className="border-gray-700" />

      {/* Conversations */}
      <div>
        <h2 className="text-lg font-semibold mb-3">💬 Conversations</h2>
        <ul className="space-y-2 max-h-[55vh] overflow-y-auto pr-1">
          {history.map((conv) => (
            <li
              key={conv.id}
              className={`p-3 rounded-lg text-sm flex justify-between items-center cursor-pointer transition-all duration-300 ease-in-out ${
                conv.id === selectedConversationId
                  ? "bg-blue-600 text-white shadow-lg"
                  : "bg-gray-800 hover:bg-blue-700 hover:text-white text-gray-300"
              }`}
            >
              <div
                onClick={() => onSelectConversation(conv.id)}
                className="truncate font-medium flex-1"
                title={conv.title}
              >
                {conv.title}
                <div className="text-xs text-gray-400 mt-1">
                  {conv.created_at
                    ? new Date(conv.created_at).toLocaleString()
                    : "Date inconnue"}
                </div>
              </div>
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  handleDeleteConversation(conv.id);
                }}
                className="ml-3 text-red-500 hover:text-red-700 font-bold text-lg p-1 rounded-lg transition hover:bg-red-100 hover:bg-opacity-30"
                title="Supprimer la conversation"
                aria-label="Supprimer la conversation"
              >
                <FiTrash2 />
              </button>
            </li>
          ))}
        </ul>
      </div>

      {/* MODALS */}

      {/* 🔍 Modal de recherche */}
      {showSearchModal && (
        <div className="fixed inset-0 bg-black bg-opacity-60 z-50 flex justify-center items-start pt-20 px-4">
          <div className="bg-white text-gray-800 w-full max-w-xl rounded-2xl shadow-2xl p-8 relative animate-fadeIn">
            <button
              onClick={() => setShowSearchModal(false)}
              className="absolute top-4 right-4 text-gray-600 hover:text-red-600 text-3xl transition"
              aria-label="Fermer la recherche"
            >
              <FiX />
            </button>
            <h3 className="text-2xl font-bold mb-6">
              🔍 Rechercher une conversation
            </h3>

            <input
              type="text"
              value={conversationSearchQuery}
              onChange={(e) => setConversationSearchQuery(e.target.value)}
              placeholder="Titre de la conversation..."
              className="w-full mb-6 px-5 py-3 border border-gray-300 rounded-lg shadow-sm focus:ring-2 focus:ring-blue-600 transition"
              autoFocus
            />

            <ul className="max-h-[300px] overflow-y-auto space-y-3 scrollbar-thin scrollbar-thumb-blue-600 scrollbar-track-gray-200">
              {history
                .filter((c) =>
                  c.title
                    .toLowerCase()
                    .includes(conversationSearchQuery.toLowerCase())
                )
                .map((conv) => (
                  <li
                    key={conv.id}
                    onClick={() => {
                      onSelectConversation(conv.id);
                      setShowSearchModal(false);
                    }}
                    className="cursor-pointer p-4 rounded-lg hover:bg-blue-100 transition-shadow border border-transparent hover:border-blue-300"
                  >
                    <div className="font-semibold">{conv.title}</div>
                    <div className="text-xs text-gray-500 mt-1">
                      {new Date(conv.created_at).toLocaleString()}
                    </div>
                  </li>
                ))}
              {history.filter((c) =>
                c.title
                  .toLowerCase()
                  .includes(conversationSearchQuery.toLowerCase())
              ).length === 0 && (
                <li className="text-center text-gray-500 py-6">
                  Aucune conversation trouvée.
                </li>
              )}
            </ul>
          </div>
        </div>
      )}

      {/* 📁 Modal documents */}
      {showLibraryModal && (
        <div className="fixed inset-0 bg-black bg-opacity-60 z-50 flex justify-center items-start pt-20 px-4">
          <div className="bg-white text-gray-900 w-full max-w-3xl max-h-[80vh] overflow-y-auto p-8 rounded-2xl shadow-2xl relative">
            <button
              onClick={() => setShowLibraryModal(false)}
              className="absolute top-4 right-4 text-gray-600 hover:text-red-600 text-3xl transition"
              aria-label="Fermer la bibliothèque"
            >
              <FiX />
            </button>
            <h3 className="text-2xl font-bold mb-6">📁 Vos documents</h3>
            {documents.length > 0 ? (
              <ul className="space-y-4">
                {documents.map((doc) => (
                  <li
                    key={doc.id}
                    onClick={() => {
                      setShowLibraryModal(false);
                      onSelectDocument?.(doc);
                    }}
                    className="border-b py-3 px-3 cursor-pointer hover:bg-gray-100 rounded-lg transition"
                  >
                    <div className="font-semibold flex items-center gap-2">
                      <FiFolder className="text-indigo-600" />{" "}
                      {doc.title || doc.filename}
                    </div>
                    <div className="text-xs text-gray-500 mt-1">
                      {doc.uploaded_at
                        ? `Indexé le ${new Date(
                            doc.uploaded_at
                          ).toLocaleDateString()}`
                        : "Date inconnue"}
                    </div>
                  </li>
                ))}
              </ul>
            ) : (
              <p className="text-center text-gray-600 text-sm mt-8">
                Aucun document n’a encore été indexé.
              </p>
            )}
          </div>
        </div>
      )}
    </aside>
  );
}

export default Sidebar;
