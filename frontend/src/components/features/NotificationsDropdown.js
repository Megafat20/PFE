import React, { useState } from "react";

const NotificationsDropdown = ({ notifications, onMarkAsRead }) => {
  const [open, setOpen] = useState(false);

  // Compte des notifications non lues
  const unreadCount = notifications.filter((n) => !n.read).length;

  return (
    <div className="relative">
      {/* Bouton cloche */}
      <button
        onClick={() => setOpen((o) => !o)}
        className="relative focus:outline-none"
        aria-label="Afficher notifications"
      >
        <svg
          className="w-6 h-6 text-gray-700 hover:text-gray-900"
          fill="none"
          stroke="currentColor"
          strokeWidth={2}
          viewBox="0 0 24 24"
          strokeLinecap="round"
          strokeLinejoin="round"
        >
          <path d="M18 8a6 6 0 0 0-12 0c0 7-3 9-3 9h18s-3-2-3-9" />
          <path d="M13.73 21a2 2 0 0 1-3.46 0" />
        </svg>

        {/* Badge */}
        {unreadCount > 0 && (
          <span className="absolute -top-1 -right-1 inline-flex items-center justify-center px-2 py-1 text-xs font-bold leading-none text-white bg-red-600 rounded-full">
            {unreadCount}
          </span>
        )}
      </button>

      {/* Dropdown */}
      {open && (
        <div className="absolute right-0 mt-2 w-80 max-h-96 overflow-y-auto bg-white border border-gray-200 rounded shadow-lg z-50">
          {notifications.length === 0 && (
            <p className="p-4 text-center text-gray-500">Aucune notification</p>
          )}

          <ul>
            {notifications.map((n) => (
              <li
                key={n._id}
                className={`p-3 border-b cursor-pointer hover:bg-gray-100 ${
                  !n.read ? "bg-indigo-50 font-semibold" : ""
                }`}
              >
                <a
                  href={n.link}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="block truncate"
                  title={n.title}
                  onClick={() => {
                    if (!n.read) {
                      onMarkAsRead(n._id);
                    }
                  }}
                >
                  {n.title}
                </a>

                <div className="text-xs text-gray-500 truncate">
                  {n.message}
                </div>

                <time className="text-xs text-gray-400 italic">
                  {new Date(n.created_at).toLocaleString()}
                </time>

                <button
                  onClick={() => onMarkAsRead(n._id)}
                  className="text-indigo-600 hover:text-indigo-900 font-medium underline text-sm mt-1"
                >
                  Marquer comme lu
                </button>
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
};

export default NotificationsDropdown;
