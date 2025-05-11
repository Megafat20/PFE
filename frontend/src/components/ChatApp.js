import React, { useState, useEffect, useRef } from 'react';
import { useNavigate } from 'react-router-dom';
import Sidebar from './Sidebar'; 
import io from 'socket.io-client';
import '../App.css';

const socket = io('http://localhost:5000');

function App({ user, setUser }) {
  const [userInput, setUserInput] = useState('');
  const [chat, setChat] = useState([]);
  const [currentResponse, setCurrentResponse] = useState('');
  const [file, setFile] = useState(null);
  const [fileThumbnail, setFileThumbnail] = useState('');
  const [isUploading, setIsUploading] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState([]);
  const [showResultsBubble, setShowResultsBubble] = useState(false);
  const [showWelcomeBubble, setShowWelcomeBubble] = useState(false);  // Initialement masquée
  const [isUserLoggedIn, setIsUserLoggedIn] = useState(false); // Variable d'état pour vérifier si l'utilisateur est connecté
  const messagesEndRef = useRef(null);
  const [loading, setLoading] = useState(true);
  const navigate = useNavigate();

  useEffect(() => {
    const verifyToken = async () => {
      const token = localStorage.getItem('authToken');
      if (!token) return navigate('/');

      try {
        const response = await fetch('http://localhost:5000/protected/me', {
          headers: { Authorization: `Bearer ${token}` },
        });

        if (!response.ok) {
          localStorage.removeItem('authToken');
          navigate('/');
        } else {
          const userData = await response.json();
          setUser(userData);
          localStorage.setItem('user', JSON.stringify(userData));
          setShowWelcomeBubble(true);
          setLoading(false);
        }
      } catch (err) {
        console.error("Erreur de vérification du token:", err);
        localStorage.removeItem('authToken');
        navigate('/');
      }
    };

    verifyToken();
  }, [navigate, setUser]);

  useEffect(() => {
    const storedUser = localStorage.getItem('user');
    const isLoggedIn = localStorage.getItem('isLoggedIn') === 'true'; // Vérifier si l'utilisateur est marqué comme connecté

    if (storedUser && isLoggedIn) {
      setUser(JSON.parse(storedUser)); // Charger l'utilisateur si connecté
      setIsUserLoggedIn(true); // Marquer que l'utilisateur est connecté
      setShowWelcomeBubble(true); // Afficher la bulle de bienvenue
    }

    const token = localStorage.getItem('authToken');  // Adjust this if you're storing the token ailleurs

    if (!token) {
      console.error("Token manquant");
      return;
    }
    // Gérer les événements du socket
    socket.on('stream_response', ({ token }) => {
      setCurrentResponse(prev => prev + token + ' ');

      const storedUser = localStorage.getItem('user');
      if (storedUser) {
        setUser(JSON.parse(storedUser)); // Assurez-vous que l'utilisateur est récupéré
      }
    });

    socket.on('document_indexed', ({ title, thumbnail }) => {
      setChat(prev => [...prev, { sender: 'bot', text: `Le document "${title}" a été indexé avec succès.` }]);
      setFileThumbnail(thumbnail);
      setIsUploading(false);
    });

    socket.on('indexation_en_cours', () => {
      setIsUploading(true);
    });

    socket.on('connect', () => {
      console.log('✅ Connecté à Flask-SocketIO');
    });

    return () => {
      socket.off('stream_response');
      socket.off('document_indexed');
      socket.off('indexation_en_cours');
      socket.off('connect');
    };
  }, []);


  useEffect(() => {
    if (showWelcomeBubble) {
      const timer = setTimeout(() => {
        setShowWelcomeBubble(false);
        localStorage.setItem('isLoggedIn', 'false'); // Optionnel, pour que la bulle disparaisse après un certain temps
      }, 3000); // 3 secondes
      return () => clearTimeout(timer); // Nettoyage du timer
    }
  }, [showWelcomeBubble]);

  useEffect(() => {
    if (messagesEndRef.current) {
      messagesEndRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }, [chat]);

  useEffect(() => {
    if (!currentResponse) return;
    let index = 0;
    const tokens = currentResponse.trim().split(' ');
    let displayedText = '';

    const interval = setInterval(() => {
      if (index < tokens.length) {
        displayedText += tokens[index] + ' ';
        setChat(prev => {
          const newChat = [...prev];
          if (newChat.length > 0 && newChat[newChat.length - 1].sender === 'botTyping') {
            newChat[newChat.length - 1].text = displayedText;
          } else {
            newChat.push({ sender: 'botTyping', text: displayedText });
          }
          return newChat;
        });
        index++;
      } else {
        clearInterval(interval);
        setChat(prev =>
          prev.map(msg =>
            msg.sender === 'botTyping' ? { sender: 'bot', text: msg.text.trim() } : msg
          )
        );
        setCurrentResponse('');
      }
    }, 50);

    return () => clearInterval(interval);
  }, [currentResponse]);

  const handleSubmit = (e) => {
    
    e.preventDefault();
    if (userInput.trim() === '') return;
    const token = localStorage.getItem('authToken'); // ou autre nom selon ton app

    setChat([...chat, { sender: 'user', text: userInput }]);
    setCurrentResponse('');
    socket.emit('message', {
      user_input: userInput,
      token: token,
    });
    setUserInput('');
  };

  const handleLogout = () => {
    localStorage.removeItem('user'); // Supprimer l'utilisateur du localStorage
    setUser(null); // Réinitialiser l'état de l'utilisateur
    setShowWelcomeBubble(false); // Masquer la bulle de bienvenue après la déconnexion
  };

  const handleClearChat = () => {
    setChat([]);
  };

  const handleFileChange = (e) => {
    const selectedFile = e.target.files[0];
    if (selectedFile) {
      setFile(selectedFile);
    }
  };

  const handleFileUpload = (e) => {
    e.preventDefault();
    if (!file) return;

    const token = localStorage.getItem('authToken');
    const formData = new FormData();
    formData.append('file', file);
    setIsUploading(true);

    fetch('http://localhost:5000/upload', {
      method: 'POST',
      body: formData,
      headers: { Authorization: `Bearer ${token}` },
    })
      .then((res) => res.json())
      .then((data) => {
        setChat((prev) => [...prev, { sender: 'bot', text: data.message || data.error }]);
        setFile(null);
        setIsUploading(false);
      })
      .catch((err) => {
        console.error("Erreur d'upload:", err);
        setChat((prev) => [...prev, { sender: 'bot', text: "❌ Erreur lors du téléversement." }]);
        setIsUploading(false);
      });
  };

  const fetchDocuments = (source, query) => {
    if (!query) return;
    const token = localStorage.getItem('authToken');

    fetch('http://localhost:5000/fetch_documents', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${token}` },
      body: JSON.stringify({ query, source, max_results: 5 }),
    })
      .then((res) => res.json())
      .then((data) => {
        setSearchResults(data.documents || []);
        setShowResultsBubble(true);
      })
      .catch((err) => {
        console.error("Erreur recherche:", err);
        setSearchResults([{ title: "Erreur", summary: "Problème lors de la récupération.", url: "#" }]);
        setShowResultsBubble(true);
      });
  };


  const handleSearchSubmit = (e) => {
    e.preventDefault();
    if (!searchQuery.trim()) return;
    setChat([...chat, { sender: 'user', text: searchQuery }]); // Afficher la recherche de l'utilisateur
    fetchDocuments('arxiv', searchQuery); // Appel à la fonction fetchDocuments pour récupérer les articles
  };

  

  return (
    <div className="app-layout">
      {/* Bulle de bienvenue */}
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
      />

      <div className="chat-container">
        <h2 className="chat-header">🧠 Assistant IA</h2>
        <div className="chat-box">
          {chat.map((msg, i) => (
            <div key={i} className={`message ${msg.sender}`}>
              {msg.text}
            </div>
          ))}
          {currentResponse && (
            <div className="message bot typing-indicator">
              <span className="dot"></span>
              <span className="dot"></span>
              <span className="dot"></span>
            </div>
          )}
          <div ref={messagesEndRef} />
        </div>

        <form onSubmit={handleSubmit} className="input-form">
          <input
            type="text"
            value={userInput}
            onChange={(e) => setUserInput(e.target.value)}
            placeholder="Posez une question..."
            className="user-input"
          />
          <button type="submit" className="submit-btn" disabled={!userInput.trim() || currentResponse}>
            Envoyer
          </button>
          <button type="button" className="clear-btn" onClick={handleClearChat}>🗑️ Effacer</button>
        </form>

        {fileThumbnail && (
          <div className="file-preview">
            <img src={`http://localhost:5000/uploads/${fileThumbnail}`} alt="Aperçu du fichier" className="file-thumbnail" />
          </div>
        )}

        <form onSubmit={handleFileUpload} className="file-upload-form">
          <input type="file" onChange={handleFileChange} accept=".pdf, .txt" className="file-input" />
          <button type="submit" className="upload-btn" disabled={isUploading}>
            {isUploading ? <div className="spinner"></div> : 'Téléverser un fichier'}
          </button>
        </form>
      </div>
    </div>
  );
}

export default App;
