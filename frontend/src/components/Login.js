import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom'; // Importer useNavigate pour la redirection
import './AuthStyles.css';

const Login = ({ onLogin }) => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false); // Ajout d'un état de chargement
  const navigate = useNavigate(); // Utiliser useNavigate pour la redirection

  const handleLogin = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError('');
  
    try {
      const res = await fetch('http://localhost:5000/auth/login', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password }),
      });
  
      if (res.ok) {
        const data = await res.json();
        console.log("Token reçu : ", data.token);
        console.log("Utilisateur reçu : ", data.user);  // Ajout du log pour vérifier les données
  
        // Stockage des données utilisateur et du token dans localStorage
        localStorage.setItem('authToken', data.token);
        localStorage.setItem('user', JSON.stringify(data.user));
        localStorage.setItem('isLoggedIn', 'true');
  
        // Mise à jour de l'état de l'application avec les informations utilisateur
        onLogin(data.user);
  
        // Redirection vers la page ChatApp
        navigate('/ChatApp');
      } else {
        const data = await res.json();
        setError(data.error || 'Erreur de connexion');
      }
    } catch (err) {
      console.error('Erreur réseau:', err);
      setError('Erreur réseau');
    } finally {
      setLoading(false);
    }
  };
  


  return (
    <div className="auth-container">
      <h2 className="auth-heading">Connexion</h2>
      {error && <p className="error">{error}</p>}
      <form onSubmit={handleLogin} className="auth-form">
        
        <input
          type="email"
          placeholder="Adresse email"
          value={email}
          onChange={(e) => setEmail(e.target.value)}
          required
          className="auth-input"
        />
        <input
          type="password"
          placeholder="Mot de passe"
          value={password}
          onChange={(e) => setPassword(e.target.value)}
          required
          className="auth-input"
        />
        <button 
          type="submit" 
          className="auth-button" 
          disabled={loading} // Désactiver le bouton pendant le chargement
        >
          {loading ? 'Connexion...' : 'Se connecter'}
        </button>
      </form>
      <div className="auth-footer">
        <p>Pas encore de compte? <a href="/register">S'inscrire</a></p>
      </div>
    </div>
  );
};

export default Login;
