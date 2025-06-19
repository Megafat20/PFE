import React from "react";
import { Link, useLocation } from "react-router-dom";
import { ScrollText, Languages, RefreshCcw, Star } from "lucide-react";

const Sidebar = () => {
  const location = useLocation();

  const links = [
    { to: "/search", label: "Recherche", icon: <ScrollText className="w-6 h-6" /> },
    { to: "/summary", label: "Résumé", icon: <ScrollText className="w-6 h-6" /> },
    { to: "/translate", label: "Traduction", icon: <Languages className="w-6 h-6" /> },
    { to: "/paraphrase", label: "Paraphrase", icon: <RefreshCcw className="w-6 h-6" /> },
  ];

  return (
    <aside className="fixed top-20 left-0 w-60 h-[calc(100vh-64px)] bg-white border-r shadow-md flex flex-col p-6 overflow-y-auto z-50">
      <nav className="flex flex-col gap-4 flex-grow">
        {links.map(({ to, label, icon }) => {
          const isActive = location.pathname === to;
          return (
            <Link
              key={to}
              to={to}
              className={`flex items-center gap-3 px-4 py-3 rounded-lg font-medium transition-colors duration-200
                ${isActive
                  ? "bg-blue-600 text-white shadow-lg"
                  : "text-gray-700 hover:bg-blue-100 hover:text-blue-600"}
              `}
            >
              {icon}
              <span>{label}</span>
            </Link>
          );
        })}
      </nav>

      {/* Bouton Favoris tout en bas, avant le footer */}
      <nav className="mt-6">
        <Link
          to="/favorites"
          className={`flex items-center gap-3 px-4 py-3 rounded-lg font-medium transition-colors duration-200
            ${location.pathname === "/favorites"
              ? "bg-yellow-500 text-white shadow-lg"
              : "text-yellow-600 hover:bg-yellow-100 hover:text-yellow-700"}
          `}
        >
          <Star className="w-6 h-6" />
          <span>Favoris</span>
        </Link>
      </nav>

      <footer className="text-sm text-gray-400 mt-auto pt-6 border-t select-none">
        &copy; {new Date().getFullYear()} MonAssistant IA
      </footer>
    </aside>
  );
};

export default Sidebar;
