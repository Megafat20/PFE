import React, { useState,useEffect } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { Menu, X } from 'lucide-react'; // Icônes modernes
import './Navbar.css';
import NotificationsDropdown from "./features/NotificationsDropdown";

function Navbar({ user, onLogout }) {
  const [isOpen, setIsOpen] = useState(false);
  const navigate = useNavigate();

  const handleLogout = () => {
    onLogout();
    localStorage.removeItem("user");
    localStorage.removeItem("authToken");
    navigate('/');
  };

  const toggleMenu = () => setIsOpen(!isOpen);


  const [notifications, setNotifications] = useState([]);

  useEffect(() => {
    // Charger les notifications ici depuis ton API
    const fetchNotifications = async () => {
      const token = localStorage.getItem("authToken");
      const res = await fetch("http://localhost:5000/notifications", {
        headers: { Authorization: `Bearer ${token}` },
      });
      if (res.ok) {
        const data = await res.json();
        setNotifications(data);
      }
    };
    fetchNotifications();
  }, []);

  const markAsRead = async (id) => {
    const token = localStorage.getItem("authToken");
    await fetch(`http://localhost:5000/notifications/${id}/read`, {
      method: "POST",
      headers: { Authorization: `Bearer ${token}` },
    });
    setNotifications((prev) => prev.map(n => n._id === id ? { ...n, read: true } : n));
  };
  return (
<nav className="bg-white  shadow-md border-b border-gray-200 fixed top-0 left-0 right-0 z-50">
  <div className="max-w-9xl mx-auto px-6 py-3 flex items-center justify-between">
    
    {/* Logo à gauche */}
    <div className="flex items-center flex-1">
      <h2 className="text-2xl font-bold text-blue-600 flex items-center gap-2 select-none">
        <span className="text-5xl">🤖</span> Mon Assistant IA
      </h2>
    </div>

    {/* Menu Desktop & Bouton Mobile à droite */}
    <div className="flex items-center gap-4">
      {/* Menu Desktop */}
      <ul className="hidden md:flex items-center gap-6 text-gray-700 font-medium">
        <li><Link to="/" className="hover:text-blue-600 transition">Accueil</Link></li>

        {user ? (
          <>
            <li><Link to="/ChatApp" className="bg-blue-600 text-white px-3 py-1 rounded hover:bg-blue-700 transition">Chat</Link></li>
            <li className="text-blue-600 font-semibold">Bienvenue, {user.name}!</li>
            <li><Link to="/search" className="hover:text-blue-600 transition">Multi-fonctions</Link></li>
            <li>
              <button onClick={handleLogout} className="bg-red-500 hover:bg-red-600 text-white px-3 py-1 rounded transition">
                Déconnexion
              </button>
            </li>
            <NotificationsDropdown notifications={notifications} onMarkAsRead={markAsRead} />
          </>
        ) : (
          <>
            <li><Link to="/login" className="hover:text-blue-600 transition">Connexion</Link></li>
            <li><Link to="/register" className="hover:text-blue-600 transition">Inscription</Link></li>
          </>
        )}
      </ul>

      {/* Bouton Menu Mobile */}
      <div className="md:hidden">
        <button onClick={toggleMenu} className="text-gray-700">
          {isOpen ? <X size={28} /> : <Menu size={28} />}
        </button>
      </div>
    </div>
  </div>

  {/* Menu Mobile */}
  {isOpen && (
    <ul className="md:hidden bg-white px-6 pb-4 pt-2 space-y-2 text-gray-700 font-medium shadow-md">
      <li><Link to="/" onClick={toggleMenu} className="block hover:text-blue-600">Accueil</Link></li>

      {user ? (
        <>
          <li><Link to="/ChatApp" onClick={toggleMenu} className="block text-blue-600">Chat</Link></li>
          <li className="block text-blue-600 font-semibold">Bienvenue, {user.name}!</li>
          <li><Link to="/multifeature" onClick={toggleMenu} className="block hover:text-blue-600">Multi-fonctions</Link></li>
          <li>
            <button onClick={() => { handleLogout(); toggleMenu(); }} className="w-full text-left text-red-500 hover:text-red-600">
              Déconnexion
            </button>
          </li>
        </>
      ) : (
        <>
          <li><Link to="/login" onClick={toggleMenu} className="block hover:text-blue-600">Connexion</Link></li>
          <li><Link to="/register" onClick={toggleMenu} className="block hover:text-blue-600">Inscription</Link></li>
        </>
      )}
    </ul>
  )}
</nav>

  );
}

export default Navbar;
