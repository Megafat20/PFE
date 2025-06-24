import React from "react";
import { Link, useLocation } from "react-router-dom";
import { ScrollText, Languages, RefreshCcw, Star } from "lucide-react";

const Sidebar = () => {
  const location = useLocation();

  const links = [
    { to: "/search", label: "Recherche", icon: <ScrollText className="w-5 h-5" /> },
    { to: "/summary", label: "Résumé", icon: <ScrollText className="w-5 h-5" /> },
    { to: "/translate", label: "Traduction", icon: <Languages className="w-5 h-5" /> },
    { to: "/paraphrase", label: "Paraphrase", icon: <RefreshCcw className="w-5 h-5" /> },
    { to: "/ReportGenerator", label: "Report", icon: <RefreshCcw className="w-5 h-5" /> },

  ];

  return (
    <aside className="fixed top-20 left-0 w-80 h-[calc(100vh-64px)] text-white bg-gray-900 dark:bg-gray-900 border-r border-gray-800 p-6 overflow-y-auto z-50">
      <h2 className="text-xl font-semibold text-white dark:text-gray-100 mb-6 tracking-tight">
        🧠 Assistant IA
      </h2>

      <nav className="flex flex-col gap-2">
        {links.map(({ to, label, icon }) => {
          const isActive = location.pathname === to;
          return (
            <Link
              key={to}
              to={to}
              className={`flex items-center gap-3 px-4 py-3 rounded-lg font-medium transition-all duration-200
                ${
                  isActive
                    ? "bg-blue-600 text-white shadow-md"
                    : "text-white dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-800 hover:text-blue-600"
                }
              `}
            >
              {icon}
              <span>{label}</span>
            </Link>
          );
        })}
      </nav>

      {/* 📁 Favoris */}
      <div className="mt-8 pt-4 border-t border-gray-200 dark:border-gray-700">
        <Link
          to="/favorites"
          className={`flex items-center gap-3 px-4 py-3 rounded-lg font-medium transition-all duration-200 mt-2
            ${
              location.pathname === "/favorites"
                ? "bg-yellow-500 text-white shadow-md"
                : "text-yellow-700 dark:text-yellow-400 hover:bg-yellow-100 dark:hover:bg-yellow-900"
            }
          `}
        >
          <Star className="w-5 h-5" />
          <span>Favoris</span>
        </Link>
      </div>

      <footer className="text-xs text-gray-400 dark:text-gray-500 mt-auto pt-6 border-t border-gray-200 dark:border-gray-800 select-none">
        &copy; {new Date().getFullYear()} MonAssistant IA
      </footer>
    </aside>
  );
};

export default Sidebar;
