import React, { useState } from 'react';

export default function ChatModal({ document, onClose }) {
  const [question, setQuestion] = useState('');
  const [messages, setMessages] = useState([]); // { from: 'user'|'bot', text: '' }
  const [loading, setLoading] = useState(false);
  const [answer, setAnswer] = useState("");
  const [sources, setSources] = useState([]);
  const sendQuestion = async () => {
    if (!question.trim()) return;

  
    setLoading(true);
    try {
      const token = localStorage.getItem('AuthToken'); // ou sessionStorage

      const res = await fetch('http://localhost:5000/chat_document', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          document_id: document._id,
          question: question,
        }),
      });
      
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error || 'Erreur serveur');
      }
      const data = await res.json();
      setMessages((msgs) => [...msgs, { from: 'bot', text: data.answer }]);
      setSources(data.sources || []);
   
    } catch (err) {
      setMessages((msgs) => [
        ...msgs,
        { from: 'bot', text: `Erreur : ${err.message}` },
      ]);
    } finally {
      setLoading(false);
      setQuestion('');
    }
  };

  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex justify-center items-center z-50">
      <div className="bg-white rounded-lg w-full max-w-xl max-h-[80vh] flex flex-col p-4 shadow-lg">
        <div className="flex justify-between items-center mb-4">
          <h2 className="text-xl font-bold text-gray-800">💬 Chat sur : {document.filename}</h2>
          <button
            onClick={onClose}
            className="text-gray-500 hover:text-gray-700 font-bold text-xl"
          >
            &times;
          </button>
        </div>

        <div className="flex-1 overflow-y-auto mb-4 p-2 border rounded bg-gray-50">
          {messages.length === 0 && (
            <p className="text-gray-500 italic">Pose ta question sur ce document...</p>
          )}
          {messages.map((msg, idx) => (
            <div
              key={idx}
              className={`mb-2 p-2 rounded ${
                msg.from === 'user'
                  ? 'bg-blue-200 text-blue-900 self-end ml-auto'
                  : 'bg-gray-300 text-gray-900 self-start mr-auto'
              }`}
              style={{ maxWidth: '75%' }}
            >
              {msg.text}
            </div>
          ))}
        </div>
        {sources.length > 0 && (
        <div className="mt-2 p-2 border-t text-sm text-gray-600">
          <h3 className="font-semibold mb-1">Sources utilisées :</h3>
          <ul className="list-disc pl-5 space-y-1">
            {sources.map((src, idx) => (
              <li key={idx}>
                {src.content.slice(0, 200)}...{" "}
                {src.metadata?.source && (
                  <span className="italic text-xs">({src.metadata.source})</span>
                )}
              </li>
            ))}
          </ul>
        </div>
         )}
        <div className="flex gap-2">
          <input
            type="text"
            value={question}
            onChange={(e) => setQuestion(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === 'Enter') sendQuestion();
            }}
            placeholder="Tape ta question..."
            className="flex-1 p-2 border rounded"
            disabled={loading}
          />
          <button
            onClick={sendQuestion}
            disabled={loading}
            className="bg-blue-600 text-white px-4 rounded hover:bg-blue-700 disabled:opacity-50"
          >
            {loading ? '...' : 'Envoyer'}
          </button>
        </div>
      </div>
    </div>
  );
}
