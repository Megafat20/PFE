import React from "react";

// ProgressBar.jsx (ou .js)
const ProgressBar = ({ progress }) => {
    return (
      <div
        style={{
          position: "fixed",
          top: 0,
          left: 0,
          height: "3px",
          width: `${progress}%`,
          backgroundColor: "#2563EB",
          transition: "width 0.2s ease",
          zIndex: 9999,
          pointerEvents: "none",
        }}
      />
    );
  };
  

export default ProgressBar;
