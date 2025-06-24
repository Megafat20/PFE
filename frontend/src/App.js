import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';

import { AuthProvider } from './components/AuthProvider';
import AuthInitializer from './components/AuthInitializer';

import Home from './components/Home';
import Login from './components/Login';
import Register from './components/Register';
import ChatApp from './components/ChatApp';
import Summary from './components/features/Summary';
import Translate from './components/features/Translate';
import Paraphrase from './components/features/Paraphrase';
import ReportGenerator from './components/features/ReportGenerator';

import Favorites from './components/features/Favorites';
import MultiFeatureApp from './components/features/MultiFeatureApp';
import Search from './components/features/Search';

import LayoutWithSidebar from './components/layouts/LayoutWithSidebar';
import PlainLayout from './components/layouts/PlainLayout';
import { ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';
import './i18n/i18n';

import { useAuth } from './components/AuthProvider';

function AppRoutes() {
  const { user, loading } = useAuth();

  const handleLogout = () => {
    localStorage.removeItem('authToken');
    window.location.href = '/login';
  };

  if (loading) return <div>Chargement...</div>;

  return (
    <Routes>
      <Route path="/" element={<PlainLayout user={user} onLogout={handleLogout}><Home /></PlainLayout>} />
      <Route path="/login" element={<PlainLayout user={user} onLogout={handleLogout}><Login /></PlainLayout>} />
      <Route path="/register" element={<PlainLayout user={user} onLogout={handleLogout}><Register /></PlainLayout>} />
      <Route path="/chatapp" element={<PlainLayout user={user} onLogout={handleLogout}><ChatApp /></PlainLayout>} />

      {/* Authenticated pages with sidebar */}
      <Route path="/search" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Search /></LayoutWithSidebar>} />
      <Route path="/summary" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Summary /></LayoutWithSidebar>} />
      <Route path="/translate" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Translate /></LayoutWithSidebar>} />
      <Route path="/paraphrase" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Paraphrase /></LayoutWithSidebar>} />
      <Route path="/ReportGenerator" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><ReportGenerator /></LayoutWithSidebar>} />
      <Route path="/favorites" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><Favorites /></LayoutWithSidebar>} />
      <Route path="/multifeature" element={<LayoutWithSidebar user={user} onLogout={handleLogout}><MultiFeatureApp /></LayoutWithSidebar>} />
    </Routes>
  );
}

function App() {
  return (
    <Router>
      <AuthProvider>
        <AuthInitializer />
        <ToastContainer />
        <AppRoutes />
      </AuthProvider>
    </Router>
  );
}

export default App;
