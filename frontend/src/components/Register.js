import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import './AuthStyles.css';

const Register = ({ onRegister }) => {
  const [email, setEmail] = useState('');
  const [name, setName] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState('');
  const navigate = useNavigate();

  const handleRegister = async (e) => {
    e.preventDefault();
    try {
      const res = await fetch('http://localhost:5000/auth/register', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ name, email, password }),
      });

      const data = await res.json();

      if (res.ok) {
        // Appeler la fonction onRegister pour mettre à jour l'état de l'utilisateur dans App.js
        onRegister({ name, email }); // Passez l'objet utilisateur
        navigate('/ChatApp'); // Redirige vers ChatApp
      } else {
        setError(data.error || 'Erreur lors de l\'inscription');
      }
    } catch (err) {
      setError('Erreur réseau');
    }
  };

  return (
    <div className="auth-container">
      <h2 className="auth-heading">Inscription</h2>
      {error && <p className="error">{error}</p>}
      <form onSubmit={handleRegister} className="auth-form">
        <input
          type="text"
          placeholder="Nom"
          value={name}
          onChange={(e) => setName(e.target.value)}
          required
          className="auth-input"
        />
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
        <button type="submit" className="auth-button">S'inscrire</button>
      </form>
      <div className="auth-footer">
        <p>Déjà un compte? <a href="/login">Se connecter</a></p>
      </div>
    </div>
  );
};

export default Register;
