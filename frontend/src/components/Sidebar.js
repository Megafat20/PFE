import React, { useState, useEffect } from 'react';

function Sidebar({
    searchQuery,
    setSearchQuery,
    handleSearchSubmit,
    showResultsBubble,
    setShowResultsBubble,
    searchResults,
    onSelectConversation,
    
    selectedConversationId  // <-- on récupère l'ID sélectionné en prop
}) {
    const [history, setHistory] = useState([]);
    const [selectedSource, setSelectedSource] = useState('');
    const sourceIcons = {
        arxiv: '🅰️',
        pubmed: '🅿️',
        semantic_scholar: '🔍',
        OpenAlex: '📚',
        unknown: '❓'
      };
    const fetchHistory = async () => {
        const token = localStorage.getItem('authToken');
        try {
            const response = await fetch('http://localhost:5000/conversations', {
                headers: { Authorization: `Bearer ${token}` },
            });
            if (response.ok) {
                const data = await response.json();
                setHistory(data || []);
            } else {
                console.error('Erreur lors de la récupération de l\'historique');
            }
        } catch (err) {
            console.error('Erreur réseau :', err);
        }
    };

    useEffect(() => {
        fetchHistory();
    }, []);

    const handleNewConversation = async () => {
        const token = localStorage.getItem('authToken');
        try {
            const response = await fetch('http://localhost:5000/conversations', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${token}`,
                },
                body: JSON.stringify({ title: 'Nouvelle conversation' }),
            });

            if (response.ok) {
                const data = await response.json();
                await fetchHistory();
                onSelectConversation(data.conversation_id);
            } else {
                console.error('Échec de la création de la conversation');
            }
        } catch (err) {
            console.error('Erreur réseau lors de la création :', err);
        }
    };

    const handleConversationClick = (conversationId) => {
        onSelectConversation(conversationId);
    };

    return (
        <aside className="w-80 flex-shrink-0 h-full bg-zinc-900 border-r shadow-lg p-5 space-y-6 overflow-y-auto">
      
          {/* Section Recherche */}
          <div>
            <h2 className="text-xl font-bold text-slate-200 mb-4 flex items-center gap-2">📚 Recherche</h2>
      
            <form onSubmit={handleSearchSubmit} className="space-y-3">
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                placeholder="Ex: deep learning, LLM..."
                className="w-full px-4 py-2 border rounded-lg shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
              />
      
              <select
                value={selectedSource}
                onChange={(e) => setSelectedSource(e.target.value)}
                className="w-full  px-4 py-2 border rounded-lg shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                <option value="">Toutes les sources</option>
                <option value="all">Toutes</option>
                <option value="arxiv">Arxiv</option>
                <option value="pubmed">PubMed</option>
                <option value="semantic_scholar">Semantic Scholar</option>
                <option value="openalex">OpenAlex</option>
              </select>
      
              <div className="flex flex-col sm:flex-row gap-2">
                <button
                  type="submit"
                  className="w-full bg-blue-600 hover:bg-blue-700 text-white py-2 rounded-lg transition"
                >
                  🔍 Rechercher
                </button>
                <button
                  type="button"
                  onClick={() => setShowResultsBubble(prev => !prev)}
                  className="w-full bg-gray-100 hover:bg-gray-200 text-gray-700 py-2 rounded-lg transition"
                >
                  {showResultsBubble ? '🔽 Replier' : '📄 Résultats'}
                </button>
              </div>
            </form>
          </div>
      
          {/* Section Conversations */}
          <div>
            <h2 className="text-xl font-bold text-slate-200 mb-3 flex items-center gap-2">💬 Discussions</h2>
      
            <button
              onClick={handleNewConversation}
              className="w-full mb-3 bg-green-600 hover:bg-green-700 text-white py-2 rounded-lg transition"
            >
              ➕ Nouvelle conversation
            </button>
      
            <ul className="flex flex-col overflow-y-auto max-h-[60vh]">
            {history.map((conversation) => (
               <li
               key={conversation.id}
               onClick={() => handleConversationClick(conversation.id)}
               className={`cursor-pointer p-3 rounded mb-1 transition-colors duration-200
                 ${
                   conversation.id === selectedConversationId
                     ? "bg-sky-600 text-white shadow-lg border border-blue-400"
                     : "hover:bg-sky-600 text-gray-800"
                 }`}
             >
                <div className="flex flex-col">
                    <strong className="truncate text-gray-100 font-semibold">{conversation.title}</strong>
                    <small className="text-xs text-gray-300">
                    {new Date(conversation.created_at).toLocaleString()}
                    </small>
                </div>
                </li>
            ))}
            </ul>
          </div>
      
          {/* Résultats de recherche flottants */}
          {showResultsBubble && (
            <div className="fixed inset-0 bg-black bg-opacity-30 z-40 flex justify-center items-start pt-24">
              <div className="bg-white w-full max-w-3xl max-h-[80vh] overflow-y-auto p-6 rounded-lg shadow-2xl relative animate-fadeIn">
                <button
                  className="absolute top-3 right-3 text-gray-600 hover:text-red-500 text-2xl font-bold"
                  onClick={() => setShowResultsBubble(false)}
                >
                  ×
                </button>
      
                <h3 className="text-xl font-semibold mb-4 flex items-center gap-2">📄 Résultats de recherche</h3>
      
                {searchResults.length > 0 ? (
                  <ul className="space-y-5">
                    {searchResults
                      .filter(doc => selectedSource === '' || selectedSource === 'all' || doc.source === selectedSource)
                      .map((doc, idx) => (
                        <li key={idx} className="border-b pb-3">
                          <div className="flex items-center gap-2 mb-1 text-sm">
                            <span title={doc.source}>
                              {sourceIcons[doc.source] || sourceIcons.unknown}
                            </span>
                            <span className="font-bold">{doc.title}</span>
                          </div>
                          <a
                            href={doc.url}
                            target="_blank"
                            rel="noopener noreferrer"
                            className="text-blue-600 hover:underline text-sm"
                          >
                            📥 Lire le PDF
                          </a>
                          <p className="text-sm text-gray-600 mt-1">{doc.summary?.slice(0, 100)}...</p>
                        </li>
                      ))}
                  </ul>
                ) : (
                  <p className="text-sm text-gray-600">Aucun résultat trouvé.</p>
                )}
              </div>
            </div>
          )}
        </aside>
      );
      
      
}

export default Sidebar;
