import React from "react";
import Sidebar2 from "../features/Sidebar2";
import Navbar from "../Navbar";

const LayoutWithSidebar = ({ children, user, onLogout }) => {
  return (
    <>
      <Navbar user={user} onLogout={onLogout} />
      <Sidebar2 />

      <main
        className="pt-16 px-0 min-h-screen bg-gray-50 p-8"
        style={{ minHeight: "calc(100vh - 64px)" }}
      >
        {children}
      </main>
    </>
  );
};

export default LayoutWithSidebar;
