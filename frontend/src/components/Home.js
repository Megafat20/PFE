import React from 'react';
import { Link } from 'react-router-dom';


const Home = () => {
  return (
    <div className="min-h-screen pt-24 flex flex-col justify-between bg-gradient-to-br from-indigo-50 to-blue-100 p-6">
      {/* Header */}
      <header className="text-center mt-8">
        <h1 className="text-4xl md:text-5xl font-bold text-gray-800">
          Bienvenue sur notre Application
        </h1>
        <p className="mt-2 text-lg text-gray-600">
          Votre espace pour gérer votre expérience
        </p>
      </header>

      {/* Contenu principal */}
      <div className="flex flex-col md:flex-row justify-center items-center gap-8 flex-grow">
        {/* Carte Connexion */}
        <div className="bg-white shadow-xl rounded-xl p-6 max-w-sm w-full text-center">
          <h3 className="text-2xl font-semibold text-gray-800">Accédez à votre compte</h3>
          <p className="mt-2 text-gray-600">
            Déjà inscrit ? Connectez-vous et continuez à utiliser nos services.
          </p>
          <Link
            to="/login"
            className="inline-block mt-4 px-6 py-2 bg-indigo-600 hover:bg-indigo-700 text-white rounded-full font-medium transition"
          >
            Se connecter
          </Link>
        </div>
        {/* Carte Inscription */}
        <div className="bg-white shadow-xl rounded-xl p-6 max-w-sm w-full text-center">
          <h3 className="text-2xl font-semibold text-gray-800">Pas encore de compte ?</h3>
          <p className="mt-2 text-gray-600">
            Créez un compte pour commencer à profiter de toutes nos fonctionnalités.
          </p>
          <Link
            to="/register"
            className="inline-block mt-4 px-6 py-2 bg-green-600 hover:bg-green-700 text-white rounded-full font-medium transition"
          >
            S'inscrire
          </Link>
        </div>
      </div>

      {/* Footer */}
      <footer className="text-center py-4 border-t mt-8">
        <p className="text-gray-500 text-sm">© 2025 Application - Tous droits réservés</p>
      </footer>
    </div>
  );
};

export default Home;
