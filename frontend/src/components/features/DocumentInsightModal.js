import React, { useState , useRef, useEffect  } from "react";
import { Dialog } from "@headlessui/react";

const DocumentInsightModal = ({
  isOpen,
  onClose,
  document,
  chunks = [],
  fullText = "",
}) => {
  const [question, setQuestion] = useState("");
  const [messages, setMessages] = useState([]);
  const [loading, setLoading] = useState(false);
  const [sources, setSources] = useState([]);
  const messagesEndRef = useRef(null);
  const sendQuestion = async () => {
    if (!question.trim()) return;

    const newMessage = { from: "user", text: question };
    setMessages((prev) => [...prev, newMessage]);
    setLoading(true);

    try {
      const token = localStorage.getItem("authToken");

      const res = await fetch("http://localhost:5000/chat_document", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`,
        },
        body: JSON.stringify({
          document_id: document._id,
          question: question,
        }),
      });

      const data = await res.json();
      const answer = data.answer || "Aucune réponse disponible.";
      setMessages((prev) => [...prev, { from: "bot", text: answer }]);
      setSources(data.sources || []);
    } catch (err) {
      setMessages((prev) => [
        ...prev,
        { from: "bot", text: `Erreur : ${err.message}` },
      ]);
    } finally {
      setLoading(false);
      setQuestion("");
    }
  };

  useEffect(() => {
    if (messagesEndRef.current) {
      messagesEndRef.current.scrollIntoView({ behavior: "smooth" });
    }
  }, [messages]);



  return (
    <Dialog open={isOpen} onClose={onClose} className="fixed inset-0 z-50">
      <div className="flex items-center justify-center min-h-screen bg-black bg-opacity-40 p-4">
        <Dialog.Panel className="bg-white w-full max-w-6xl rounded-lg shadow-lg p-6 overflow-y-auto max-h-[90vh]">
          <Dialog.Title className="text-xl font-bold mb-4 text-center">
          📌 Passages pertinents
          </Dialog.Title>
          <div className="space-y-4">
            {chunks.length === 0 && (
              <p className="text-gray-600">Aucun passage pertinent trouvé.</p>
            )}
            {chunks.map((chunk, index) => (
              <div
                key={index}
                className="border-l-4 border-blue-500 pl-4 text-sm bg-gray-50 rounded p-2"
              >
                <strong>Passage {index + 1}</strong>
                <p>{chunk.text}</p>
              </div>
            ))}
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 h-[65vh] overflow-hidden mt-2">
            {/* Zone Document */}
            <div className="border rounded p-4 overflow-y-auto bg-gray-50 text-sm leading-relaxed">
              <h3 className="text-lg font-semibold mb-2">📄 Aperçu du document</h3>
              {highlightChunks(fullText, chunks)}
            </div>

            {/* Zone Chat */}
            <div className="flex flex-col border rounded p-4 bg-white">
              <h3 className="text-lg font-semibold mb-2">💬 Chat avec le document</h3>

              <div className="flex-1 overflow-y-auto mb-2 space-y-2" style={{ maxHeight: "calc(65vh - 100px)" }} >
                {messages.length === 0 && (
                  <p className="text-gray-400 italic">Pose ta question...</p>
                )}
                {messages.map((msg, idx) => (
                  <div
                    key={idx}
                    className={`p-2 rounded ${
                      msg.from === "user"
                        ? "bg-blue-100 text-blue-800 ml-auto text-right"
                        : "bg-gray-200 text-gray-900 mr-auto text-left"
                    }`}
                    style={{ maxWidth: "80%" }}
                  >
                    {msg.text}
                  </div>
                ))}
              </div>

              {sources.length > 0 && (
                <div className="mt-2 border-t pt-2 text-sm text-gray-600">
                  <h4 className="font-semibold">🔍 Sources utilisées :</h4>
                  <ul className="list-disc pl-5 space-y-1">
                    {sources.map((src, idx) => (
                      <li key={idx}>
                        {src.content.slice(0, 150)}...
                        {src.metadata?.source && (
                          <span className="italic text-xs">
                            {" "}
                            ({src.metadata.source})
                          </span>
                        )}
                      </li>
                    ))}
                  </ul>
                </div>
              )}

              {/* Input */}
              <div className="mt-3 flex gap-2">
                <input
                  type="text"
                  value={question}
                  onChange={(e) => setQuestion(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && sendQuestion()}
                  placeholder="Ex: Quelle est la conclusion ?"
                  className="flex-1 p-2 border rounded text-sm"
                  disabled={loading}
                />
                <button
                  onClick={sendQuestion}
                  disabled={loading}
                  className="bg-blue-600 hover:bg-blue-700 text-white px-4 rounded"
                >
                  {loading ? "..." : "Envoyer"}
                </button>
              </div>
            </div>
          </div>

          <button
            onClick={onClose}
            className="mt-6 bg-gray-700 hover:bg-gray-800 text-white px-4 py-2 rounded w-full"
          >
            Fermer
          </button>
        </Dialog.Panel>
      </div>
    </Dialog>
  );
};

// Helper : surlignage des chunks dans le texte
const highlightChunks = (text, chunks) => {
  if (!text) return null;

  let lastIndex = 0;
  const elements = [];

  chunks.forEach((chunk, index) => {
    const { start, end } = chunk;

    if (start > lastIndex) {
      elements.push(
        <span key={`text-${index}`}>{text.slice(lastIndex, start)}</span>
      );
    }

    elements.push(
      <span
        key={`highlight-${index}`}
        className="bg-yellow-200 font-medium px-1 rounded"
        id={`highlight-${index}`}
      >
        {text.slice(start, end)}
      </span>
    );

    lastIndex = end;
  });

  if (lastIndex < text.length) {
    elements.push(<span key="end">{text.slice(lastIndex)}</span>);
  }

  return elements;
};

export default DocumentInsightModal;
