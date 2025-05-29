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
  const [sidebarOpen, setSidebarOpen] = useState(false);
  // --- Vérification token & récupération user ---
  useEffect(() => {
    const verifyToken = async () => {
      const token = localStorage.getItem('authToken');
      if (!token) {
        navigate('/');
        return;
      }
      try {
        const res = await fetch('http://localhost:5000/user/me', {
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
        const res = await axios.get(`http://localhost:5000/conversations/${selectedConversationId}`, {
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
    <div className="flex h-screen bg-gray-100  pt-16">
      <div id="progress-bar" />
  
      {/* Bulle de bienvenue */}
      {showWelcomeBubble && (
        <div className="fixed top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2 bg-blue-50 border border-blue-500 rounded p-5 shadow-lg text-lg font-bold z-50 animate-fadeIn">
          Bienvenue dans l'assistant IA 👋
        </div>
      )}
        <button
        className="md:hidden absolute top-4 left-4 z-50 p-2 bg-blue-600 text-white rounded"
        onClick={() => setSidebarOpen(prev => !prev)}
        aria-label="Toggle sidebar"
      >
        {sidebarOpen ? '✕' : '☰'}
      </button>
      {/* Sidebar */}
      <Sidebar
        searchQuery={searchQuery}
        setSearchQuery={setSearchQuery}
        handleSearchSubmit={handleSearchSubmit}
        showResultsBubble={showResultsBubble}
        setShowResultsBubble={setShowResultsBubble}
        searchResults={searchResults}
        onSelectConversation={handleSelectConversation}
        selectedConversationId={selectedConversationId}
      />
  
      {/* Chat container */}
      <div className="flex flex-col flex-grow p-6 bg-white overflow-y-auto">
        {/* Header */}
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-2xl font-semibold text-gray-800">🧠 Assistant IA</h2>
          <select
            onChange={e => setSelectedModel(e.target.value)}
            className="border border-gray-300 rounded p-2 text-sm"
            value={selectedModel}
          >
            <option value="llama3">Llama 3</option>
            <option value="deepseek-coder">DeepSeek-coder</option>
            <option value="mistral:7b-instruct">Mistral</option>
            <option value="gemini">Gemini</option>
          </select>
        </div>
  
        {/* Chat box */}
        <div className="flex flex-col flex-grow overflow-y-auto space-y-2 mb-4">
          {chat.length === 0 ? (
            <div className="p-5 text-center italic text-gray-500">
              Bienvenue sur MyAI, commencez une nouvelle discussion ou relancez les anciennes.
            </div>
          ) : (
            chat.map((msg, i) => (
              <div
                key={i}
                className={`max-w-[70%] p-3 rounded-lg ${
                  (msg.role || msg.sender) === 'user'
                    ? 'self-end bg-blue-100'
                    : 'self-start bg-gray-200'
                }`}
              >
                {msg.content || msg.text}
              </div>
            ))
          )}
  
          {/* Typing indicator */}
          {currentResponse && (
            <div className="flex items-center space-x-1 self-start px-4 py-2 bg-gray-200 rounded">
              <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce" />
              <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce delay-150" />
              <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce delay-300" />
            </div>
          )}
  
          <div ref={messagesEndRef} />
        </div>
  
        {/* Input form */}
        <form onSubmit={handleSubmit} className="flex items-center space-x-2 mb-4">
          <input
            type="text"
            value={userInput}
            onChange={e => setUserInput(e.target.value)}
            placeholder="Posez une question..."
            className="flex-grow p-2 border border-gray-300 rounded"
            disabled={isUploading}
          />
          <button
            type="submit"
            className="px-4 py-2 bg-green-600 text-white rounded hover:bg-green-700 disabled:opacity-50"
            disabled={!userInput.trim() || currentResponse || isUploading}
          >
            Envoyer
          </button>
          <button
            type="button"
            onClick={handleClearChat}
            className="px-3 py-2 bg-red-500 text-white rounded hover:bg-red-600 disabled:opacity-50"
            disabled={isUploading}
          >
            🗑️ Effacer
          </button>
        </form>
  
        {/* Upload form */}
        <form
          onSubmit={handleFileUpload}
          className="flex items-center space-x-2 mb-4"
          encType="multipart/form-data"
        >
          <input
            type="file"
            multiple
            onChange={handleFilesChange}
            disabled={isUploading}
            className="p-2 border border-gray-300 rounded"
          />
          <button
            type="submit"
            className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
            disabled={files.length === 0 || !selectedConversationId || isUploading}
          >
            {isUploading ? 'Chargement...' : 'Upload'}
          </button>
        </form>
  
        {/* Uploaded documents */}
        <h3 className="text-lg font-medium mb-2">Documents Uploadés :</h3>
        <div className="flex flex-wrap gap-3">
          {uploadedFiles.map((file, index) => (
            <div key={index} className="relative inline-block">
              <img
                src={`http://localhost:5000/thumbnails/${file.user_id.$oid || file.user_id}/${file.thumbnail}`}
                alt={`Thumbnail ${index}`}
                className="w-24 h-auto rounded border border-gray-300"
              />
              <button
                onClick={() => handleDeleteFile(file._id)}
                className="absolute top-1 right-1 w-5 h-5 bg-red-600 text-white text-xs rounded-full flex items-center justify-center hover:bg-red-700"
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
