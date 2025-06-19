import React from "react";
import Sidebar2 from "../features/Sidebar2";
import Navbar from "../Navbar";

const LayoutWithSidebar = ({ children, user, onLogout }) => {
  return (
    <>
    <div className="flex flex-col min-h-screen">
    {/* Navbar fixe en haut */}
    <Navbar user={user} onLogout={onLogout} />

    <div className="flex flex-1 min-h-0  pt-16">
      {/* Sidebar à gauche */}
      <div className="hidden md:block w-64 bg-white border-r">
        <Sidebar2 />
      </div>

      {/* Contenu principal */}
      <main
        className="flex-1 bg-gray-50 p-4 sm:p-6 md:p-8 overflow-auto"
        style={{ minHeight: "calc(100vh - 64px)" }}
      >
        {children}
      </main>
    </div>
  </div>
  </>
  );
};

export default LayoutWithSidebar;
