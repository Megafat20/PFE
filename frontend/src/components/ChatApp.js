import React, { useState, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import Sidebar from "./Sidebar";
import io from "socket.io-client";
import axios from "axios";
import "../App.css";
import MarkdownRenderer from "./MarkdownRenderer";

function App({ user, setUser }) {
  const [socket, setSocket] = useState(null);
  const [userInput, setUserInput] = useState("");
  const [chat, setChat] = useState([]);
  
  const [currentResponse, setCurrentResponse] = useState("");
  const [files, setFiles] = useState([]);
  const [thumbnails, setThumbnails] = useState([]);
  const [isUploading, setIsUploading] = useState(false);
  const [searchQuery, setSearchQuery] = useState("");
  const [isTyping, setIsTyping] = useState(false);
  const [searchResults, setSearchResults] = useState([]);
  const [showResultsBubble, setShowResultsBubble] = useState(false);
  const [showWelcomeBubble, setShowWelcomeBubble] = useState(false);
  const [loading, setLoading] = useState(true);
  const [selectedConversationId, setSelectedConversationId] = useState(null);
  const [uploadedFiles, setUploadedFiles] = useState([]);
  const [selectedModel, setSelectedModel] = useState("llama3");
  const [selectedSource, setSelectedSource] = useState("all");
  const [ragMode, setRagMode] = useState("standard");
  const messagesEndRef = useRef(null);
  const [buffer, setBuffer] = useState("");
  const animationRef = useRef(null);
  const navigate = useNavigate();
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [ratings, setRatings] = React.useState({}); // { messageIndex: rating }
  const [submitteds, setSubmitteds] = React.useState({});
  const [forceLLM, setForceLLM] = React.useState(false);
  const [persona, setPersona] = useState("default");
  const [context, setContext] = useState("default");
  // --- Vérification token & récupération user ---
  useEffect(() => {
    const verifyToken = async () => {
      const token = localStorage.getItem("authToken");
      if (!token) {
        navigate("/login");
        return;
      }
      try {
        const res = await fetch("http://localhost:5000/protected/me", {
          headers: { Authorization: `Bearer ${token}` },
        });
        if (!res.ok) throw new Error("Token invalide");
        const userData = await res.json();
        setUser(userData);
      } catch (err) {
        console.error("Erreur auth:", err);
        localStorage.removeItem("authToken");
        navigate("/login");
      } finally {
        setLoading(false);
      }
    };
    verifyToken();
  }, [navigate]);

  // --- Initialisation socket.io ---
  useEffect(() => {
    if (loading) return; // éviter d'ouvrir socket avant user chargé
    const token = localStorage.getItem("authToken");
    if (!token) {
      navigate("/login");
      return;
    }
    const s = io("http://localhost:5000", {
      transports: ["websocket"],
      auth: { token },
    });
    setSocket(s);
    return () => s.disconnect();
  }, [loading, navigate]);


  
  
  useEffect(() => {
    if (!socket) return;
  
    const onStreamResponse = ({ token }) => {
      console.log("Token reçu:", token);
      setIsTyping(true);
      setBuffer((prev) => prev + token);
    };
    const onStreamEnd = () => {
      setIsTyping(false); // ⛔ fin forcée si pas de buffer
    };

    socket.on("connect", () => console.log("✅ Connecté à SocketIO"));
    socket.on("stream_response", onStreamResponse);
    socket.on("stream_end",onStreamEnd )
    return () => {
      socket.off("connect");
      socket.off("stream_response", onStreamResponse);
      socket.off("stream_end", onStreamEnd);
    };
  }, [socket]);
  
  useEffect(() => {
    if (!buffer) return;
  
    const tokens = buffer.trim().split(" ");
    let idx = 0;
  
    // Construire la phrase au fur et à mesure
    let displayed = "";
  
    const interval = setInterval(() => {
      if (idx < tokens.length) {
        displayed += tokens[idx] + " ";
  
        setChat(prev => {
          const copy = [...prev];
          if (copy.length && copy[copy.length - 1].sender === 'botTyping') {
            // Remplacer complètement le texte à chaque fois (pas de +=)
            copy[copy.length - 1].text = displayed;
          } else {
            copy.push({ sender: 'botTyping', text: displayed });
          }
          return copy;
        });
  
        idx++;
      } else {
        clearInterval(interval);
        setBuffer("");
        setChat(prev =>
          prev.map(m => m.sender === 'botTyping' ? { sender: 'bot', text: m.text.trim() } : m)
        );
        setIsTyping(false);
      }
    }, 50);
  
    return () => clearInterval(interval);
  }, [buffer]);
  
  
  

  // ✉️ Envoi du message utilisateur
  const handleSubmit = (e) => {
    e.preventDefault();
    const token = localStorage.getItem("authToken");
  
    if (!userInput.trim() || !socket) return;
  
    setChat((prev) => [
      ...prev,
      { sender: "user", text: userInput },
    ]);
    setIsTyping(true);
    setCurrentResponse("");
    setBuffer("");
  
    socket.emit("chat_message", {
      user_input: userInput,
      conversation_id: selectedConversationId,
      model: selectedModel,
      persona,
      context,
      force_llm: forceLLM,
      token: token,
    });
  
    setUserInput("");
  };
  

  // --- Scroll auto sur nouveau message ---
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [chat]);

  // --- Gestion upload fichiers ---
  const handleFilesChange = (e) => setFiles(Array.from(e.target.files));

  const handleFileUpload = async (e) => {
    e.preventDefault();
    if (files.length === 0 || !selectedConversationId) return;

    setIsUploading(true);
    const formData = new FormData();
    files.forEach((file) => formData.append("files", file));
    formData.append("conversation_id", selectedConversationId);

    try {
      const token = localStorage.getItem("authToken");
      const res = await axios.post("http://localhost:5000/upload", formData, {
        headers: {
          Authorization: `Bearer ${token}`,
          "Content-Type": "multipart/form-data",
        },
      });
      const newFiles = res.data.files;
      setUploadedFiles((prev) => [...prev, ...newFiles]);
      alert("Fichiers uploadés avec succès");
      setFiles([]);
    } catch (error) {
      console.error("Échec de l'upload:", error);
    } finally {
      setIsUploading(false);
    }
  };

  const handleDeleteFile = async (docId) => {
    const token = localStorage.getItem("authToken");
    try {
      await axios.delete(`http://localhost:5000/delete_document/${docId}`, {
        headers: { Authorization: `Bearer ${token}` },
      });

      // Mise à jour de l'état après suppression
      setUploadedFiles((prev) => {
        const newFiles = prev.filter((file) => file._id !== docId);
        setThumbnails(
          newFiles.map(
            (file) =>
              `http://localhost:5000/thumbnails/${
                file.user_id?.$oid || file.user_id
              }/${file.thumbnail}`
          )
        );
        return newFiles;
      });
    } catch (error) {
      console.error("Erreur suppression fichier :", error);
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

    const token = localStorage.getItem("authToken");

    const fetchFiles = async () => {
      try {
        const res = await axios.get(
          `http://localhost:5000/user_documents/${selectedConversationId}`,
          {
            headers: { Authorization: `Bearer ${token}` },
          }
        );
        console.log("Réponse fichiers:", res.data);
        const files = res.data.files || [];
        setUploadedFiles(files);

        // Générer les URLs des miniatures
        const newThumbnails = files.map(
          (file) =>
            `http://localhost:5000/thumbnails/${
              file.user_id.$oid || file.user_id
            }/${file.thumbnail}`
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
        const res = await axios.get(
          `http://localhost:5000/conversations/${selectedConversationId}`,
          {
            headers: { Authorization: `Bearer ${token}` },
          }
        );
        const messages = res.data.messages || [];
        setChat(messages);

        // Initialiser ratings et submitteds en fonction des messages
        const initialRatings = {};
        const initialSubmitteds = {};
        messages.forEach((msg, i) => {
          if (msg.rating) {
            initialRatings[i] = msg.rating;
            initialSubmitteds[i] = !!msg.validated;
          }
        });
        setRatings(initialRatings);
        setSubmitteds(initialSubmitteds);

        console.log("Messages récupérés :", messages);
      } catch (error) {
        console.error("Erreur récupération conversation:", error);
        setChat([]);
        setRatings({});
        setSubmitteds({});
      }
    };

    fetchFiles();
    fetchConversation();
  }, [selectedConversationId]);

  // --- Recherche documents ---
  const fetchDocuments = (source, query) => {
    if (!query) return;
    const token = localStorage.getItem("authToken");
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
    fetch("http://localhost:5000/fetch_documents", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${token}`,
      },
      body: JSON.stringify({ query, source, max_results: 5 }),
    })
      .then((res) => res.json())
      .then((data) => {
        setSearchResults(data.documents || []);
        setShowResultsBubble(true);
      })
      .catch((err) => {
        console.error("Erreur recherche:", err);
        setSearchResults([
          {
            title: "Erreur",
            summary: "Problème lors de la récupération.",
            url: "#",
          },
        ]);
        setShowResultsBubble(true);
      })
      .finally(() => {
        endProgress();
      });
  };

  const handleSearchSubmit = (e) => {
    e.preventDefault();
    if (!searchQuery.trim()) return;

    setChat((prev) => [...prev, { sender: "user", text: searchQuery }]);

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

  const handleTemplateClick = (templateType) => {
    if (templateType === "literature_review") {
      setUserInput("Write a literature review about: ");
    } else if (templateType === "methodology") {
      setUserInput("Describe the methodology for: ");
    }
  };

  const handlePersonaChange = async (e) => {
    const newPersona = e.target.value;
    console.log("Nouvelle personnalité sélectionnée :", newPersona);
    setPersona(newPersona);

    const token = localStorage.getItem("authToken");
    await fetch("http://localhost:5000/auth/update_profile", {
      method: "PUT",
      headers: {
        "Content-Type": "application/json",
        Authorization: `Bearer ${token}`,
      },
      body: JSON.stringify({
        persona: newPersona,
      }),
    });
  };

  const handleContextChange = async (e) => {
    const newContext = e.target.value;
    setContext(newContext); // met à jour le state React

    const token = localStorage.getItem("authToken");
    try {
      await fetch("http://localhost:5000/auth/update_profile", {
        method: "PUT",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          context: newContext,
        }),
      });
    } catch (error) {
      console.error("Erreur lors de la mise à jour du contexte :", error);
    }
  };

  // --- JSX ---
  return (
    <div className="flex h-screen bg-gray-100 pt-16">
      <div id="progress-bar" />

      {/* Bulle de bienvenue */}
      {showWelcomeBubble && (
        <div className="fixed top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2 bg-blue-50 border border-blue-500 rounded p-5 shadow-lg text-lg font-bold z-50 animate-fadeIn">
          Bienvenue dans l'assistant IA 👋
        </div>
      )}

      <button
        className="md:hidden absolute top-4 left-4 z-50 p-2 bg-blue-600 text-white rounded"
        onClick={() => setSidebarOpen((prev) => !prev)}
        aria-label="Toggle sidebar"
      >
        {sidebarOpen ? "✕" : "☰"}
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
          <h2 className="text-2xl font-semibold text-gray-800">
            🧠 Assistant IA
          </h2>
          <select
            onChange={(e) => setSelectedModel(e.target.value)}
            className="border border-gray-300 rounded p-2 text-sm"
            value={selectedModel}
          >
            <option value="llama3">Llama 3</option>
            <option value="deepseek-coder:6.7b">DeepSeek-coder</option>
            <option value="mistral:7b-instruct">Mistral</option>
            <option value="gemini">Gemini</option>
          </select>
        </div>

        {/* Personnalité et Contexte */}
        <div className="flex items-center space-x-4 mb-6">
          <select
            onChange={handlePersonaChange}
            className="border border-gray-300 rounded p-2 text-sm"
            value={persona}
          >
            <option value="default">🎭 Personnalité</option>
            <option value="formelle">Formelle</option>
            <option value="amicale">Amicale</option>
            <option value="concise">Concise</option>
          </select>

          <select
            onChange={handleContextChange}
            className="border border-gray-300 rounded p-2 text-sm"
            value={context}
          >
            <option value="default">⚙️ Contexte</option>
            <option value="scientifique">Scientifique</option>
            <option value="juridique">Juridique</option>
            <option value="général">Général</option>
          </select>
        </div>
        <div className="flex space-x-2 mb-4">
          <button
            onClick={() => handleTemplateClick("literature_review")}
            className="px-3 py-1 bg-purple-100 text-purple-800 rounded text-sm hover:bg-purple-200"
          >
            Literature Review
          </button>
          <button
            onClick={() => handleTemplateClick("methodology")}
            className="px-3 py-1 bg-blue-100 text-blue-800 rounded text-sm hover:bg-blue-200"
          >
            Methodology
          </button>
        </div>
        {/* Chat box */}
        <div className="flex flex-col flex-grow overflow-y-auto space-y-4 mb-4">
          {chat.length === 0 ? (
            <div className="p-5 text-center italic text-gray-500">
              Bienvenue sur MyAI, commencez une nouvelle discussion ou relancez
              les anciennes.
            </div>
          ) : (
            chat.map((msg, i) => {
              const isUser = (msg.role || msg.sender) === "user";
              const isTyping = msg.sender === "botTyping";
              // Fonction pour gérer la notation
              const handleRate = async (star, index) => {
                const token = localStorage.getItem("authToken");
              
                setRatings((prev) => ({ ...prev, [index]: star }));
              
                const payload = {
                  question: chat[index - 1]?.content || "",
                  answer: chat[index]?.content || chat[index]?.text,
                  rating: star,
                  validated: star >= 4,
                  conversation_id: selectedConversationId,
                };
              
                // Envoyer au backend si non encore soumis
                if (!submitteds[index]) {
                  try {
                    const res = await fetch("http://localhost:5000/rating/validate_answer", {
                      method: "POST",
                      headers: {
                        "Content-Type": "application/json",
                        Authorization: `Bearer ${token}`,
                      },
                      body: JSON.stringify(payload),
                    });
              
                    if (res.ok) {
                      setSubmitteds((prev) => ({ ...prev, [index]: true }));
                    } else {
                      alert("Erreur lors de la sauvegarde.");
                    }
                  } catch (err) {
                    alert("Erreur réseau.");
                  }
                }
              };

              // Composant notation étoiles inline
              const StarRating = ({ rating, onRate, disabled }) => (
                <div className="flex space-x-1 mt-2">
                  {[1, 2, 3, 4, 5].map((star) => {
                    const filled = star <= rating;
                    return (
                      <button
                        key={star}
                        type="button"
                        disabled={disabled}
                        onClick={() => !disabled && onRate(star)}
                        className={`text-2xl ${
                          filled ? "text-yellow-400" : "text-gray-300"
                        } hover:text-yellow-500 focus:outline-none`}
                        aria-label={`${star} étoile${star > 1 ? "s" : ""}`}
                      >
                        ★
                      </button>
                    );
                  })}
                </div>
              );

              return (
                <div
                  key={i}
                  className={`max-w-[70%] p-3 rounded-lg whitespace-pre-wrap ${
                    isUser ? "self-end bg-blue-100" : "self-start bg-gray-200"
                  } relative`}
                >
                  {isUser ? (
                    <p className="whitespace-pre-wrap">
                      {msg.content || msg.text}
                    </p>
                  ) : msg.sender === "botTyping" ? (
                    <div className="whitespace-pre-wrap">
                    <MarkdownRenderer content={msg.text} />
                    </div>
                  ) : (
                    <>
                      <div className="whitespace-pre-wrap">
                        <MarkdownRenderer content={msg.content || msg.text} />
                      </div>
                      <StarRating
                        rating={ratings[i] || 0}
                        onRate={(star) => handleRate(star, i)}
                        disabled={submitteds[i] || false}
                      />
                      {submitteds[i] && (
                        <div className="text-green-600 text-sm mt-1 font-semibold">
                          Merci pour votre avis !
                        </div>
                      )}
                    </>
                  )}
                </div>
              );
              
            })
          )}

        <div>
          {isTyping && (
  <div className="flex items-center space-x-1 self-start px-4 py-2 bg-gray-200 rounded">
    <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce" />
    <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce delay-150" />
    <span className="w-2 h-2 bg-gray-600 rounded-full animate-bounce delay-300" />
  </div>
)}
        </div>


          <div ref={messagesEndRef} />
        </div>
        <div className="flex items-center mb-2 space-x-2">
          <input
            id="forceLLM"
            type="checkbox"
            checked={forceLLM}
            onChange={() => setForceLLM((prev) => !prev)}
            className="w-4 h-4"
          />
          <label htmlFor="forceLLM" className="text-gray-700 select-none">
            Forcer génération LLM (ignorer documents indexés)
          </label>
        </div>
        {/* Input form */}
        <form
          onSubmit={handleSubmit}
          className="flex items-center space-x-2 mb-4"
        >
          <input
            type="text"
            value={userInput}
            onChange={(e) => setUserInput(e.target.value)}
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
            disabled={
              files.length === 0 || !selectedConversationId || isUploading
            }
          >
            {isUploading ? "Chargement..." : "Upload"}
          </button>
        </form>

        {/* Uploaded documents */}
        <h3 className="text-lg font-medium mb-2">Documents Uploadés :</h3>
        <div className="flex flex-wrap gap-3">
          {uploadedFiles
            .filter((file) => file.is_master) // Garde uniquement les masters
            .map((file, index) => (
              <div key={index} className="relative inline-block">
                <img
                  src={`http://localhost:5000/thumbnails/${
                    file.user_id.$oid || file.user_id
                  }/${file.thumbnail}`}
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
