import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';
import './AuthStyles.css';
import { useAuth } from './AuthProvider'; // ✅ utilise le contexte

const Login = () => {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);
  const { setUser } = useAuth(); // ✅ récupère setUser du contexte
  const navigate = useNavigate();

  const handleLogin = async (e) => {
    e.preventDefault();
    setLoading(true);
    try {
      const res = await fetch('http://localhost:5000/auth/login', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email, password }),
      });
      if (!res.ok) throw new Error('Identifiants invalides');

      const { token, user } = await res.json();

      localStorage.setItem('authToken', token);
      setUser(user); // ✅ enregistre dans le contexte
      localStorage.setItem("just_logged_in", "true");
      await new Promise(resolve => setTimeout(resolve, 100)); 
      navigate('/chatapp');
    } catch (err) {
      console.error('Erreur login :', err);
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen flex items-center justify-center bg-gradient-to-br from-blue-50 to-indigo-100 px-4">
      <div className="w-full max-w-md bg-white rounded-2xl shadow-lg p-8">
        <h2 className="text-2xl font-bold text-center text-gray-800 mb-6">Connexion à votre compte</h2>

        {error && (
          <div className="mb-4 text-sm text-red-600 bg-red-100 border border-red-300 rounded-md px-4 py-2">
            {error}
          </div>
        )}

        <form onSubmit={handleLogin} className="space-y-4">
          <input
            type="email"
            placeholder="Adresse email"
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            required
            className="w-full px-4 py-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-indigo-400"
          />
          <input
            type="password"
            placeholder="Mot de passe"
            value={password}
            onChange={(e) => setPassword(e.target.value)}
            required
            className="w-full px-4 py-3 border border-gray-300 rounded-lg focus:outline-none focus:ring-2 focus:ring-indigo-400"
          />
          <button
            type="submit"
            disabled={loading}
            className={`w-full py-3 rounded-lg text-white font-semibold transition duration-300 ${
              loading
                ? 'bg-indigo-300 cursor-not-allowed'
                : 'bg-indigo-600 hover:bg-indigo-700'
            }`}
          >
            {loading ? 'Connexion...' : 'Se connecter'}
          </button>
        </form>

        <p className="mt-6 text-center text-sm text-gray-600">
          Pas encore de compte ?{' '}
          <a href="/register" className="text-indigo-600 hover:underline">
            S'inscrire
          </a>
        </p>
      </div>
    </div>
  );
};

export default Login;
