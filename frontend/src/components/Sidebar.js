import React, { useState, useEffect } from 'react';

function Sidebar({ searchQuery, setSearchQuery, handleSearchSubmit, showResultsBubble, setShowResultsBubble, searchResults, onSelectConversation }) {
    const [history, setHistory] = useState([]);

    // Récupérer l'historique des conversations depuis l'API
    useEffect(() => {
        const fetchHistory = async () => {
            const token = localStorage.getItem('authToken');
            try {
                const response = await fetch('http://localhost:5000/conversations', {
                    headers: {
                        Authorization: `Bearer ${token}`,
                    },
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

        fetchHistory();
    }, []);

    const handleConversationClick = (conversationId) => {
        onSelectConversation(conversationId); // Appeler la fonction pour sélectionner une conversation
    };

    return (
        <div className="sidebar">
            <h3 className="sidebar-title">📚 Rechercher</h3>
            <form onSubmit={handleSearchSubmit} className="search-form">
                <input
                    type="text"
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    placeholder="Ex : deep learning, LLM..."
                    className="search-input"
                />
                <button type="submit" className="search-btn">Rechercher</button>
                <button
                    className="toggle-bubble-btn"
                    onClick={() => setShowResultsBubble(prev => !prev)}
                >
                    {showResultsBubble ? '🔽 Replier' : '📄 Résultats'}
                </button>
            </form>

            <h3>Historique des Discussions</h3>
            <ul>
                {history.map((conversation, index) => (
                    <li key={index} className="conversation-item">
                        <div>
                            <strong>{conversation.title}</strong>
                        </div>
                        <div>
                            <small>{new Date(conversation.created_at).toLocaleString()}</small>
                        </div>
                        <button onClick={() => handleConversationClick(conversation.id)} className="select-btn">
                            Reprendre
                        </button>
                    </li>
                ))}
            </ul>

            {showResultsBubble && (
                <div className="floating-results">
                    <button className="close-btn" onClick={() => setShowResultsBubble(false)}>×</button>
                    <h4>📄 Résultats</h4>
                    {searchResults.length > 0 ? (
                        <ul className="results-list">
                            {searchResults.map((doc, idx) => (
                                <li key={idx}>
                                    <strong>{doc.title}</strong><br />
                                    <a href={doc.url} target="_blank" rel="noopener noreferrer">📥 Lire le PDF</a>
                                    <p className="summary">{doc.summary.slice(0, 100)}...</p>
                                </li>
                            ))}
                        </ul>
                    ) : (
                        <p>Aucun résultat trouvé.</p>
                    )}
                </div>
            )}
        </div>
    );
}

export default Sidebar;
