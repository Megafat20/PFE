// App.js
import React , { useState, useEffect } from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Navbar from './components/Navbar';
import Home from './components/Home';
import ChatApp from './components/ChatApp';
import Login from './components/Login';
import Register from './components/Register';
import './App.css';



function App() {

  const [user, setUser] = useState(null);  // État pour stocker l'utilisateur


  useEffect(() => {
    const fetchUserProfile = async () => {
      const token = localStorage.getItem('authToken');
  
      if (!token) return;
  
      try {
        const res = await fetch('http://localhost:5000/protected/me', {
          method: 'GET',
          headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${token}`
          }
        });
  
        if (!res.ok) {
          throw new Error('Non autorisé');
        }
  
        const data = await res.json();
        console.log('Profil utilisateur récupéré au chargement de l\'app :', data);
        setUser(data);
      } catch (err) {
        console.error('Erreur de récupération du profil :', err.message);
        setUser(null);  // Optionnel : réinitialiser si token expiré
        localStorage.removeItem('authToken');
      }
    };
  
    fetchUserProfile();
  }, []);
  


  const handleRegister = (newUser) => {
    setUser(newUser); // Met à jour l'utilisateur après l'inscription
  };

   // Fonction pour gérer la connexion de l'utilisateur
   const handleLogin = (userData) => {
    setUser(userData); // Mettre à jour l'état utilisateur
    localStorage.setItem('user', JSON.stringify(userData)); // Sauvegarder dans localStorage
    console.log("Utilisateur connecté :", userData);
  };

  const handleLogout = () => {
    localStorage.removeItem('user');
    localStorage.removeItem('authToken');
    setUser(null);
    console.log("Utilisateur déconnecté :", user);
  };
  return (
    <Router>
      <div className="App">
        <Navbar user={user} onLogout={handleLogout} />
        <Routes>
          <Route path="/" element={<Home />} />
          <Route 
            path="/login" 
            element={<Login onLogin={handleLogin} />} 
          />
          <Route 
            path="/ChatApp" 
            element={<ChatApp user={user} setUser={setUser} />} 
          />
          <Route path="/register" element={<Register onRegister={handleRegister} />} />
        </Routes>
      </div>
    </Router>
  );
}

export default App;
