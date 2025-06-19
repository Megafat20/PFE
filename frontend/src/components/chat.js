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
  const [files, setFiles] = useState([]);
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
  const [ragMode, setRagMode] = useState('standard');
  const messagesEndRef = useRef(null);
  const navigate = useNavigate();
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [conversations, setConversations] = useState([]);
  
  // --- Token verification & user retrieval ---
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

  // --- Fetch conversations ---
  const fetchConversations = async () => {
    const token = localStorage.getItem('authToken');
    try {
      const res = await axios.get('http://localhost:5000/conversations', {
        headers: { Authorization: `Bearer ${token}` }
      });
      setConversations(res.data || []);
    } catch (error) {
      console.error("Error fetching conversations:", error);
    }
  };

  // Load conversations on mount
  useEffect(() => {
    if (!loading) {
      fetchConversations();
    }
  }, [loading]);

  // --- Socket.io initialization ---
  useEffect(() => {
    if (loading) return;
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

  // Socket event handlers
  useEffect(() => {
    if (!socket) return;

    const onStreamResponse = ({ token }) => {
      setChat(prev => {
        const newChat = [...prev];
        
        if (newChat.length > 0 && newChat[newChat.length - 1].sender === 'botTyping') {
          newChat[newChat.length - 1].text += token;
        } 
        else if (newChat.length > 0 && newChat[newChat.length - 1].sender === 'user') {
          newChat.push({ sender: 'botTyping', text: token });
        }
        else {
          newChat.push({ sender: 'botTyping', text: token });
        }
        
        return newChat;
      });
    };

    const onStreamEnd = ({ full_response }) => {
      setChat(prev => {
        const newChat = [...prev];
        if (newChat.length > 0 && newChat[newChat.length - 1].sender === 'botTyping') {
          newChat[newChat.length - 1] = { 
            sender: 'bot', 
            text: full_response 
          };
        } else {
          newChat.push({ sender: 'bot', text: full_response });
        }
        return newChat;
      });
    };

    socket.on('connect', () => console.log('✅ Connecté à SocketIO'));
    socket.on('stream_response', onStreamResponse);
    socket.on('stream_end', onStreamEnd);

    return () => {
      socket.off('connect');
      socket.off('stream_response', onStreamResponse);
      socket.off('stream_end', onStreamEnd);
    };
  }, [socket]);

  // Send message handler
  const handleSubmit = e => {
    const token= localStorage.getItem('authToken');
    e.preventDefault();
    if (!userInput.trim() || !socket) return;

    setChat(prev => [...prev, { sender: 'user', text: userInput }]);
    setUserInput('');

    socket.emit('chat_message', {
      user_input: userInput,
      conversation_id: selectedConversationId,
      model: selectedModel,
      rag_mode: ragMode,
      token: token
    });
  };

  // Auto-scroll
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [chat]);

  // File upload management
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
      const newFiles = res.data.files;
      setUploadedFiles((prev) => [...prev, ...newFiles]);
      alert('Fichiers uploadés avec succès');
      setFiles([]);
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
      setUploadedFiles(prev => prev.filter(file => file._id !== docId));
    } catch (error) {
      console.error('Erreur suppression fichier :', error);
    }
  };

  // Retrieve documents on conversation change
  useEffect(() => {
    if (!selectedConversationId) {
      setChat([]);
      setUploadedFiles([]);
      return;
    }
  
    const token = localStorage.getItem('authToken');
  
    const fetchFiles = async () => {
      try {
        // Add ?only_masters=true to get only master documents
        const res = await axios.get(
          `http://localhost:5000/user_documents/${selectedConversationId}?only_masters=true`, 
          {
            headers: { Authorization: `Bearer ${token}` }
          }
        );
        setUploadedFiles(res.data.files || []);
      } catch (error) {
        console.error("Erreur récupération fichiers:", error);
        setUploadedFiles([]);
      }
    };
  
    const fetchConversation = async () => {
      try {
        const res = await axios.get(
          `http://localhost:5000/conversations/${selectedConversationId}`, 
          {
            headers: { Authorization: `Bearer ${token}` }
          }
        );
        setChat(res.data.messages || []);
      } catch (error) {
        console.error("Erreur récupération conversation:", error);
        setChat([]);
      }
    };
  
    fetchFiles();
    fetchConversation();
  }, [selectedConversationId]);

  // Document search
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

  // Welcome bubble animation
  useEffect(() => {
    if (!showWelcomeBubble) return;
    const timer = setTimeout(() => setShowWelcomeBubble(false), 3000);
    return () => clearTimeout(timer);
  }, [showWelcomeBubble]);

  const handleSelectConversation = (id) => setSelectedConversationId(id);
  
  // Scientific writing templates
  const handleTemplateClick = (templateType) => {
    if (templateType === 'literature_review') {
      setUserInput('Write a literature review about: ');  
    } else if (templateType === 'methodology') {
      setUserInput('Describe the methodology for: ');
    }
  };
  
  // Handle conversation deletion
  const handleDeleteConversation = async (conversationId) => {
    const token = localStorage.getItem('authToken');
    try {
      await axios.delete(`http://localhost:5000/conversations/${conversationId}`, {
        headers: { Authorization: `Bearer ${token}` }
      });
      
      // Refresh conversations list
      fetchConversations();
      
      // Clear current chat if deleted conversation was active
      if (selectedConversationId === conversationId) {
        setSelectedConversationId(null);
        setChat([]);
        setUploadedFiles([]);
      }
    } catch (error) {
      console.error('Error deleting conversation:', error);
    }
  };

  return (
    <div className="flex h-screen bg-gray-100 pt-16">
      <div id="progress-bar" />
  
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
      
      <Sidebar
        searchQuery={searchQuery}
        setSearchQuery={setSearchQuery}
        handleSearchSubmit={handleSearchSubmit}
        showResultsBubble={showResultsBubble}
        setShowResultsBubble={setShowResultsBubble}
        searchResults={searchResults}
        onSelectConversation={handleSelectConversation}
        selectedConversationId={selectedConversationId}
        onDeleteConversation={handleDeleteConversation}
        conversations={conversations}
      />
  
      <div className="flex flex-col flex-grow p-6 bg-white overflow-y-auto">
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-2xl font-semibold text-gray-800">🧠 Assistant IA</h2>
          <div className="flex items-center space-x-2">
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
            <select
              onChange={e => setRagMode(e.target.value)}
              className="border border-gray-300 rounded p-2 text-sm"
              value={ragMode}
            >
              <option value="standard">Standard RAG</option>
              <option value="deep">Deep Research</option>
              <option value="summary">Summary Focused</option>
              <option value="bm25">BM25 Only</option>
            </select>
          </div>
        </div>
        
        <div className="flex space-x-2 mb-4">
          <button 
            onClick={() => handleTemplateClick('literature_review')}
            className="px-3 py-1 bg-purple-100 text-purple-800 rounded text-sm hover:bg-purple-200"
          >
            Literature Review
          </button>
          <button 
            onClick={() => handleTemplateClick('methodology')}
            className="px-3 py-1 bg-blue-100 text-blue-800 rounded text-sm hover:bg-blue-200"
          >
            Methodology
          </button>
        </div>

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
                  msg.sender === 'user'
                    ? 'self-end bg-blue-100'
                    : msg.sender === 'botTyping'
                      ? 'self-start bg-gray-300'
                      : 'self-start bg-gray-200'
                }`}
              >
                {msg.text}
              </div>
            ))
          )}
          <div ref={messagesEndRef} />
        </div>
  
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
            disabled={!userInput.trim() || isUploading}
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
  
        <h3 className="text-lg font-medium mb-2">Documents Uploadés :</h3>
        <div className="flex flex-wrap gap-3">
          {uploadedFiles.map((file, index) => (
            <div key={index} className="relative inline-block">
              <img
                src={`http://localhost:5000/thumbnails/${file.user_id.$oid || file.user_id}/${file.thumbnail}`}
                alt={`Thumbnail ${index}`}
                className="w-24 h-auto rounded border border-gray-300"
              />
              <div className="text-xs mt-1 truncate w-24">
                {file.filename}
              </div>
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