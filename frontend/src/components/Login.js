import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom'; // Importer useNavigate pour la redirection
import './AuthStyles.css';

const Login = ({ onLogin }) => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false); // Ajout d'un état de chargement
  const navigate = useNavigate(); // Utiliser useNavigate pour la redirection

   const handleLogin = async e => {
    e.preventDefault();
    try {
      const res = await fetch('http://localhost:5000/auth/login', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password }),

      });
      if (!res.ok) throw new Error('Identifiants invalides');
      const { token, user } = await res.json();
      // 1) Stocke le token
      localStorage.setItem('authToken', token);
      // 2) Passe le user en haut pour instancier le socket après
      onLogin(user);
      // 3) Redirige vers le chat
      navigate('/ChatApp');
    } catch (err) {
      console.error('Erreur réseau:', err);
      alert(err.message);
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
