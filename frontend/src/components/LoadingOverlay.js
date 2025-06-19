import '../App.css';
import React from 'react';

export default function LoadingOverlay({ progress }) {
    return (
      <div className="fixed inset-0 bg-black bg-opacity-50 flex flex-col justify-center items-center z-50">
        <div className="bg-white rounded p-6 shadow-lg flex flex-col items-center">
          <div className="loader mb-4" /> {/* tu peux faire un spinner CSS ou emoji */}
          <p className="text-lg font-semibold mb-2">Indexation en cours...</p>
          <p>{progress}%</p>
        </div>
      </div>
    );
  }
  