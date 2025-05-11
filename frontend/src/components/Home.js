import React from 'react';
import './HomeStyles.css'; // Assurez-vous de lier votre fichier CSS

const Home = () => {
  return (
    <div className="home-container">
      <header className="home-header">
        <h1>Bienvenue sur notre Application</h1>
        <p className="home-subheading">Votre espace pour gérer votre expérience</p>
      </header>

      <div className="home-content">
        <div className="home-card">
          <h3>Accédez à votre compte</h3>
          <p>Déjà inscrit? Connectez-vous et continuez à utiliser nos services.</p>
          <a href="/login" className="home-button">Se connecter</a>
        </div>

        <div className="home-card">
          <h3>Pas encore de compte?</h3>
          <p>Créez un compte pour commencer à profiter de toutes nos fonctionnalités.</p>
          <a href="/register" className="home-button">S'inscrire</a>
        </div>
      </div>

      <footer className="home-footer">
        <p>© 2025 Application - Tous droits réservés</p>
      </footer>
    </div>
  );
};

export default Home;
