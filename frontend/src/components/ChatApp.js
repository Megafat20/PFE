import React, { useState, useEffect, useRef } from 'react';
import { useNavigate } from 'react-router-dom';
import Sidebar from './Sidebar'; 
import io from 'socket.io-client';
import axios from 'axios';
import '../App.css';

function App({ user, setUser }) {
  const [socket, setSocket] = useState(null);
  const [userInput, setUserInput] = useState('');
  const [chat, setChat] = useState([]);
  const [currentResponse, setCurrentResponse] = useState('');
  const [files, setFiles] = useState([]);
  const [thumbnails, setThumbnails] = useState([]);
  const [isUploading, setIsUploading] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState([]);
  const [showResultsBubble, setShowResultsBubble] = useState(false);
  const [showWelcomeBubble, setShowWelcomeBubble] = useState(false);
  const [loading, setLoading] = useState(true);
  const [selectedConversationId, setSelectedConversationId] = useState(null);
  const [uploadedFiles, setUploadedFiles] = useState([]);
  const [selectedModel, setSelectedModel] = useState('llama3');
  const [selectedSource, setSelectedSource] = useState('all');
  const messagesEndRef = useRef(null);
  const [buffer, setBuffer] = useState('');
  const animationRef = useRef(null);
  const navigate = useNavigate();
 
  // --- Vérification token & récupération user ---
  useEffect(() => {
    const verifyToken = async () => {
      const token = localStorage.getItem('authToken');
      if (!token) {
        navigate('/');
        return;
      }
      try {
        const res = await fetch('http://localhost:5000/protected/me', {
          headers: { Authorization: `Bearer ${token}` },
        });
        if (!res.ok) throw new Error('Token invalide');
        const userData = await res.json();
        setUser(userData);
        setShowWelcomeBubble(true);
      } catch (err) {
        localStorage.removeItem('authToken');
        navigate('/');
      } finally {
        setLoading(false);
      }
    };
    verifyToken();
  }, [navigate, setUser]);

  // --- Initialisation socket.io ---
  useEffect(() => {
    if (loading) return; // éviter d'ouvrir socket avant user chargé
    const token = localStorage.getItem('authToken');
    if (!token) {
      navigate('/login');
      return;
    }
    const s = io('http://localhost:5000', {
      transports: ['websocket'],
      auth: { token }
    });
    setSocket(s);
    return () => s.disconnect();
  }, [loading, navigate]);

  // 🧠 Affichage lettre par lettre
  useEffect(() => {
    if (!buffer) return;

    const revealNextChar = () => {
      setCurrentResponse(prev => prev + buffer[0]);
      setBuffer(prev => prev.slice(1));

      if (buffer.length > 1) {
        animationRef.current = setTimeout(revealNextChar, 30); // vitesse ici (ms)
      }
    };

    animationRef.current = setTimeout(revealNextChar, 30);

    return () => clearTimeout(animationRef.current);
  }, [buffer]);

  // ✅ Connexion socket et réception token
  useEffect(() => {
    if (!socket) return;

    const onStreamResponse = ({ token }) => {
      setBuffer(prev => prev + token); // Ajout dans le buffer
    };

    socket.on('connect', () => console.log('✅ Connecté à SocketIO'));
    socket.on('stream_response', onStreamResponse);

    return () => {
      socket.off('connect');
      socket.off('stream_response', onStreamResponse);
    };
  }, [socket]);

  // ✉️ Envoi du message utilisateur
  const handleSubmit = e => {
    const token= localStorage.getItem('authToken');
    e.preventDefault();
    if (!userInput.trim() || !socket) return;

    setChat(prev => [...prev, { sender: 'user', text: userInput }]);
    setCurrentResponse('');
    setBuffer('');

    socket.emit('chat_message', {
      user_input: userInput,
      conversation_id: selectedConversationId,
      model: selectedModel,
      token:token
    });

    setUserInput('');
  };


  // --- Affichage progressif réponse bot (mot par mot) ---
  useEffect(() => {
    if (!currentResponse) return;

    let idx = 0;
    const tokens = currentResponse.trim().split(' ');
    let displayed = '';

    const interval = setInterval(() => {
      if (idx < tokens.length) {
        displayed += tokens[idx] + ' ';
        setChat(prev => {
          const copy = [...prev];
          if (copy.length && copy[copy.length - 1].sender === 'botTyping') {
            copy[copy.length - 1].text = displayed;
          } else {
            copy.push({ sender: 'botTyping', text: displayed });
          }
          return copy;
        });
        idx++;
      } else {
        clearInterval(interval);
        setChat(prev => prev.map(m =>
          m.sender === 'botTyping' ? { sender: 'bot', text: m.text.trim() } : m
        ));
        setCurrentResponse('');
      }
    }, 50);

    return () => clearInterval(interval);
  }, [currentResponse]);

  // --- Scroll auto sur nouveau message ---
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [chat]);

  // --- Gestion upload fichiers ---
  const handleFilesChange = (e) => setFiles(Array.from(e.target.files));

      
  const handleFileUpload = async e => {
    e.preventDefault();
    if (files.length === 0 || !selectedConversationId) return;

    setIsUploading(true);
    const formData = new FormData();
    files.forEach(file => formData.append('files', file));
    formData.append('conversation_id', selectedConversationId);

    try {
      const token = localStorage.getItem('authToken');
      const res = await axios.post('http://localhost:5000/upload', formData, {
        headers: {
          Authorization: `Bearer ${token}`,
          'Content-Type': 'multipart/form-data',
        }
      });
      const newFiles = res.data.files; // Fichiers renvoyés par le backend avec user_id et thumbnail
      setUploadedFiles((prev) => [...prev, ...newFiles]);
      alert('Fichiers uploadés avec succès');
      setFiles([]); // Reset fichiers sélectionnés
  } catch (error) {
    console.error('Échec de l\'upload:', error);
  } finally {
    setIsUploading(false);
  }
  };
  
  const handleDeleteFile = async (docId) => {
    const token = localStorage.getItem('authToken');
    try {
      await axios.delete(`http://localhost:5000/delete_document/${docId}`, {
        headers: { Authorization: `Bearer ${token}` }
      });
  
      // Mise à jour de l'état après suppression
      setUploadedFiles(prev => {
        const newFiles = prev.filter(file => file._id !== docId);
        setThumbnails(
          newFiles.map(file => 
            `http://localhost:5000/thumbnails/${file.user_id?.$oid || file.user_id}/${file.thumbnail}`
          )
        );
        return newFiles;
      });
      
    } catch (error) {
      console.error('Erreur suppression fichier :', error);
    }
  };

  // --- Récupération fichiers uploadés à chaque changement de conversation ---
  useEffect(() => {
    if (!selectedConversationId) {
      setChat([]);
      setUploadedFiles([]);
      setThumbnails([]);
      return;
    }
  
    const token = localStorage.getItem('authToken');
  
    const fetchFiles = async () => {
      try {
        const res = await axios.get(`http://localhost:5000/user_documents/${selectedConversationId}`, {
          headers: { Authorization: `Bearer ${token}` }
        });
        console.log("Réponse fichiers:", res.data);
        const files = res.data.files || [];
        setUploadedFiles(files);
  
        // Générer les URLs des miniatures
        const newThumbnails = files.map(file => 
          `http://localhost:5000/thumbnails/${file.user_id.$oid || file.user_id}/${file.thumbnail}`
        );
        setThumbnails(newThumbnails);
  
      } catch (error) {
        console.error("Erreur récupération fichiers:", error);
        setUploadedFiles([]);
        setThumbnails([]);
      }
    };
  
    const fetchConversation = async () => {
      try {
        const res = await axios.get(`http://localhost:5000/conversation/${selectedConversationId}`, {
          headers: { Authorization: `Bearer ${token}` }
        });
        setChat(res.data.messages || []);
        console.log("Messages récupérés :", res.data.messages);
      } catch (error) {
        console.error("Erreur récupération conversation:", error);
        setChat([]);
      }
    };
  
    fetchFiles();
    fetchConversation();
  }, [selectedConversationId]);

  // --- Recherche documents ---
  const fetchDocuments = (source, query) => {
    if (!query) return;
    const token = localStorage.getItem('authToken');
    const bar = document.getElementById("progress-bar");
    let width = 0;
  
    const startProgress = () => {
      if (bar) {
        bar.style.width = "0%";
        bar.style.opacity = "1";
        bar.interval = setInterval(() => {
          if (width < 90) {
            width += 1;
            bar.style.width = `${width}%`;
          }
        }, 30);
      }
    };
  
    const endProgress = () => {
      if (bar) {
        clearInterval(bar.interval);
        bar.style.width = "100%";
        setTimeout(() => {
          bar.style.opacity = "0";
          bar.style.width = "0%";
        }, 500);
      }
    };
  
    startProgress();
    fetch('http://localhost:5000/fetch_documents', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        Authorization: `Bearer ${token}`
      },
      body: JSON.stringify({ query, source, max_results: 5 }),
    })
      .then(res => res.json())
      .then(data => {
        setSearchResults(data.documents || []);
        setShowResultsBubble(true);
      })
      .catch(err => {
        console.error("Erreur recherche:", err);
        setSearchResults([{ title: "Erreur", summary: "Problème lors de la récupération.", url: "#" }]);
        setShowResultsBubble(true);
      }).finally(() => {
        endProgress();
      });
  };

  const handleSearchSubmit = (e) => {
    e.preventDefault();
    if (!searchQuery.trim()) return;
  
    setChat(prev => [...prev, { sender: 'user', text: searchQuery }]);
  
    fetchDocuments(selectedSource, searchQuery);
  };

  const handleClearChat = () => setChat([]);

  // --- Animation bulle bienvenue ---
  useEffect(() => {
    if (!showWelcomeBubble) return;
    const timer = setTimeout(() => setShowWelcomeBubble(false), 3000);
    return () => clearTimeout(timer);
  }, [showWelcomeBubble]);

  const handleSelectConversation = (id) => setSelectedConversationId(id);

  // --- JSX ---
  return (
    <div className="app-layout">
      <div id="progress-bar"></div>
      {showWelcomeBubble && (
        <div className="welcome-bubble">
          Bienvenue, {user ? user.name : 'utilisateur'} ! 😊
        </div>
      )}

      <Sidebar
        searchQuery={searchQuery}
        setSearchQuery={setSearchQuery}
        handleSearchSubmit={handleSearchSubmit}
        showResultsBubble={showResultsBubble}
        setShowResultsBubble={setShowResultsBubble}
        searchResults={searchResults}
        onSelectConversation={handleSelectConversation}
      />

      <div className="chat-container">
        <div className="chat-header-container">
          <h2 className="chat-header">🧠 Assistant IA</h2>
          <select
            onChange={e => setSelectedModel(e.target.value)}
            className="model-selector"
            value={selectedModel}
          >
            <option value="llama3">Llama 3</option>
            <option value="deepseek-coder">DeepSeek-coder</option>
            <option value="mistral">Mistral</option>
          </select>
        </div>

        <div className="chat-box">
        {chat.length === 0 ? (
    <div className="welcome-message" style={{ 
      padding: '20px', 
      color: '#666', 
      fontStyle: 'italic', 
      textAlign: 'center' 
    }}>
      Bienvenue sur MyAI, commencez une nouvelle discussion ou relancez les anciennes.
    </div>
  ) : (
    chat.map((msg, i) => (
      <div key={i} className={`message ${msg.role || msg.sender}`}>
        {msg.content || msg.text}
      </div>
    ))
  )}
          {currentResponse && (
            <div className="message bot typing-indicator">
              <span className="dot"></span><span className="dot"></span><span className="dot"></span>
            </div>
          )}
          <div ref={messagesEndRef} />
        </div>

        <form onSubmit={handleSubmit} className="input-form">
          <input
            type="text"
            value={userInput}
            onChange={e => setUserInput(e.target.value)}
            placeholder="Posez une question..."
            className="user-input"
            disabled={isUploading}
          />
          <button
            type="submit"
            className="submit-btn"
            disabled={!userInput.trim() || currentResponse || isUploading}
          >
            Envoyer
          </button>
          <button
            type="button"
            className="clear-btn"
            onClick={handleClearChat}
            disabled={isUploading}
          >🗑️ Effacer
                  </button>
    </form>

    <form
      onSubmit={handleFileUpload}
      className="upload-form"
      encType="multipart/form-data"
    >
      <input
        type="file"
        multiple
        onChange={handleFilesChange}
        disabled={isUploading}
        className="file-input"
      />
      <button
        type="submit"
        disabled={files.length === 0 || !selectedConversationId || isUploading}
        className="upload-btn"
      >
        {isUploading ? 'Chargement...' : 'Upload'}
      </button>
    </form>
    <h3>Documents Uploadés :</h3>
    <div style={{ marginTop: 10, display: 'flex', gap: 10, flexWrap: 'wrap' }}>
  {uploadedFiles.map((file, index) => (
    <div key={index} style={{ position: 'relative', display: 'inline-block' }}>
      <img
        src={`http://localhost:5000/thumbnails/${file.user_id.$oid || file.user_id}/${file.thumbnail}`}
        alt={`Thumbnail ${index}`}
        style={{ width: 100, height: 'auto', borderRadius: 8, border: '1px solid #ccc' }}
      />
      <button
        onClick={() => handleDeleteFile(file._id)}
        style={{
          position: 'absolute',
          top: 2,
          right: 2,
          background: 'red',
          color: 'white',
          border: 'none',
          borderRadius: '50%',
          width: 20,
          height: 20,
          cursor: 'pointer'
        }}
      >
        ×
      </button>
    </div>
  ))}
</div>
  </div>
</div>
);
}

export default App;
