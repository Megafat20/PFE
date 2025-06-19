// App.js
import React, { useState, useEffect } from 'react';
import { BrowserRouter as Router, Routes, Route, useLocation } from 'react-router-dom';

import Navbar from './components/Navbar';
import Home from './components/Home';
import ChatApp from './components/ChatApp';
import Login from './components/Login';
import Register from './components/Register';
import MultiFeatureApp from './components/features/MultiFeatureApp';
import Summary from './components/features/Summary';
import Translate from './components/features/Translate';
import Paraphrase from './components/features/Paraphrase';
import Favorites from './components/features/Favorites';

import LayoutWithSidebar from './components/layouts/LayoutWithSidebar';
import PlainLayout from './components/layouts/PlainLayout';
import Search from './components/features/Search';
import { AuthProvider } from './components/AuthProvider';
import { ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';


import './i18n/i18n';

function App() {
  const [user, setUser] = useState(null);

  useEffect(() => {
    const token = localStorage.getItem('authToken');
    if (!token) return;

    const fetchUserProfile = async () => {
      try {
        const res = await fetch('http://localhost:5000/protected/me', {
          headers: { 'Authorization': `Bearer ${token}` }
        });
        if (!res.ok) throw new Error('Non autorisé');
        const data = await res.json();
        setUser(data);
      } catch (err) {
        setUser(null);
        localStorage.removeItem('authToken');
      }
    };

    fetchUserProfile();
  }, []);

  const handleRegister = (newUser) => setUser(newUser);
  const handleLogin = (userData) => {
    setUser(userData);
    localStorage.setItem('user', JSON.stringify(userData));
  };
  const handleLogout = () => {
    localStorage.removeItem('user');
    localStorage.removeItem('authToken');
    setUser(null);
  };

  return (
    
    <Router>
    <AuthProvider>
    <ToastContainer />
      <Routes>
        {/* Layout sans sidebar */}
        <Route path="/login" element={<PlainLayout user={user} onLogout={handleLogout}><Login onLogin={handleLogin} /></PlainLayout>} />
        <Route path="/register" element={<PlainLayout user={user} onLogout={handleLogout}><Register onRegister={handleRegister} /></PlainLayout>} />
        <Route path="/ChatApp" element={<PlainLayout user={user} onLogout={handleLogout}><ChatApp user={user} setUser={setUser} /></PlainLayout>} />
        <Route path="/" element={<PlainLayout user={user} onLogout={handleLogout}><Home /></PlainLayout>} />

        {/* Layout avec sidebar */}
        <Route path="/search" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Search /></LayoutWithSidebar>} />
        <Route path="/summary" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Summary /></LayoutWithSidebar>} />
        <Route path="/translate" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Translate /></LayoutWithSidebar>} />
        <Route path="/paraphrase" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Paraphrase /></LayoutWithSidebar>} />
        <Route path="/favorites" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Favorites /></LayoutWithSidebar>} />
        <Route path="/multifeature" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><MultiFeatureApp /></LayoutWithSidebar>} />
        
      </Routes>
    </AuthProvider>
  </Router>
  );
}

export default App;
