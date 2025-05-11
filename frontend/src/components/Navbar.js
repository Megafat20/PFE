import React from 'react';
import { Link,useNavigate  } from 'react-router-dom';
import './Navbar.css'; // Ajoutez vos styles ici

function Navbar({ user, onLogout }) {
  const navigate = useNavigate(); 
  const handleLogout = () => {
    
    onLogout();
    localStorage.removeItem("user");
    localStorage.removeItem("authToken");  
    navigate('/'); // Redirige vers la page d'accueil après la déconnexion
  };
  return (
    <nav className="navbar">
      <h2 className="navbar-logo">🤖 Mon Assistant IA</h2>
      <ul className="navbar-links">
        <li><Link to="/">Accueil</Link></li>
        {user ? (
          <>
            <li className="welcome-message">Bienvenue, {user.name}!</li> {/* Afficher le nom de l'utilisateur */}
            <li>
              <button onClick={handleLogout} className="logout-button">Déconnexion</button>
            </li>
          </>
        ) : (
          <>
            <li><Link to="/login">Connexion</Link></li>
            <li><Link to="/register">Inscription</Link></li>
          </>
        )}
      </ul>
    </nav>
  );
}

export default Navbar;
