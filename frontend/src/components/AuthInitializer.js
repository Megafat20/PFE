// AuthInitializer.js
import { useEffect } from 'react';
import { useNavigate, useLocation } from 'react-router-dom';
import { io } from 'socket.io-client';
import { useAuth } from './AuthProvider';

const publicPaths = ['/', '/login', '/register'];

export default function AuthInitializer() {
  const { setUser, setSocket, setLoading } = useAuth();
  const navigate = useNavigate();
  const location = useLocation();

  useEffect(() => {
    const init = async () => {
      const token = localStorage.getItem('authToken');

      if (!token && !publicPaths.includes(location.pathname)) {
        navigate('/login');
        setLoading(false);
        return;
      }

      if (!token && publicPaths.includes(location.pathname)) {
        setLoading(false);
        return;
      }

      try {
        const res = await fetch('http://localhost:5000/protected/me', {
          headers: { Authorization: `Bearer ${token}` }
        });
        if (!res.ok) throw new Error('Token invalide');
        const userData = await res.json();
        setUser(userData);

        const socket = io('http://localhost:5000', {
          transports: ['websocket'],
          auth: { token }
        });
        setSocket(socket);
      } catch (err) {
        console.error('Erreur auth:', err);
        localStorage.removeItem('authToken');
        if (!publicPaths.includes(location.pathname)) {
          navigate('/login');
        }
      } finally {
        setLoading(false);
      }
    };

    init();
  }, [location.pathname]);

  return null; // Rien à afficher
}
