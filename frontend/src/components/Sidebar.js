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
        <div className="sidebar">
            {/* Barre de recherche */}
            <h3 className="sidebar-title">📚 Rechercher</h3>
            <form onSubmit={handleSearchSubmit} className="search-form">
                <input
                    type="text"
                    value={searchQuery}
                    onChange={(e) => setSearchQuery(e.target.value)}
                    placeholder="Ex : deep learning, LLM..."
                    className="search-input"
                />
                <select
                value={selectedSource}
                onChange={(e) => setSelectedSource(e.target.value)}
                className="custom-select"
            >
                <option value="">sélectionné une source</option>
                <option value="all">Toutes sources</option>
                <option value="arxiv">Arxiv</option>
                <option value="pubmed">PubMed</option>
                <option value="semantic_scholar">Semantic Scholar</option>
                <option value="openalex">OpenAlex</option>
            </select>
                <button type="submit" className="search-btn">Rechercher</button>
                <button
                    type="button"
                    className="toggle-bubble-btn"
                    onClick={() => setShowResultsBubble(prev => !prev)}
                >
                    {showResultsBubble ? '🔽 Replier' : '📄 Résultats'}
                </button>
            </form>

            {/* Historique des conversations */}
            <h3 className="section-title">💬 Discussions récentes</h3>
            <button onClick={handleNewConversation} className="btn-new-conversation">
                ➕ Nouvelle conversation
            </button>

            <ul className="conversation-history">
                {history.map((conversation) => (
                    <li
                        key={conversation.id}
                        className={`conversation-item ${conversation.id === selectedConversationId ? 'selected' : ''}`}
                        onClick={() => handleConversationClick(conversation.id)}
                        style={{ cursor: 'pointer' }}
                    >
                        <div className="conversation-content">
                            <strong className="conversation-title">{conversation.title}</strong>
                            <small className="conversation-date">
                                {new Date(conversation.created_at).toLocaleString()}
                            </small>
                        </div>
                    </li>
                ))}
            </ul>

            {/* Résultats de recherche */}
            {showResultsBubble && (
                <div className="floating-results">
                    <button className="close-btn" onClick={() => setShowResultsBubble(false)}>×</button>
                    <h4 className="result-title">📄 Résultats</h4>
                    {searchResults.length > 0 ? (
                    <ul className="results-list">
                        {searchResults
                        .filter(doc => selectedSource === '' || selectedSource === 'all' || doc.source === selectedSource)
                        .map((doc, idx) => (
                            <li key={idx}>
                           <span className="source-icon" title={doc.source.charAt(0).toUpperCase() + doc.source.slice(1).replace('_', ' ')}
                            style={{ marginRight: '6px' }}>{sourceIcons[doc.source] || sourceIcons.unknown}</span>
                            <strong>{doc.title}</strong><br />
                            <a href={doc.url} target="_blank"  rel="noopener noreferrer" className="pdf-link">
                            📥 Lire le PDF
                            </a>
                            <p className="summary">{doc.summary?.slice(0, 100)}...</p>
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
