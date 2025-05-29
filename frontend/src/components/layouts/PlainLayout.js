
import React from "react";
import Navbar from "../Navbar";

const PlainLayout = ({ children, user, onLogout }) => {
    return (
      <div className="min-h-screen flex flex-col">
        <Navbar user={user} onLogout={onLogout} />
        {/* padding-top si navbar fixe */}
        <main className="flex-grow p-6 pt-16 bg-gray-50">
          {children}
        </main>
      </div>
    );
  };

export default PlainLayout;
