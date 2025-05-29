
import React from "react";
import Sidebar2 from "../features/Sidebar2";
import Navbar from "../Navbar";

const MainLayout = ({ children, user, onLogout }) => {
  return (
    <div className="flex">
      <Sidebar2 />
      <div className="flex-1">
        <Navbar user={user} onLogout={onLogout} />
        <main className="p-4">{children}</main>
      </div>
    </div>
  );
};

export default MainLayout;
