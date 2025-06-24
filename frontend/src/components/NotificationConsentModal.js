import { useState } from "react";


function NotificationConsentModal({ onClose }) {
  const [enabled, setEnabled] = useState(null);
  const [interests, setInterests] = useState([]);
  const token = localStorage.getItem("authToken");

  const allInterests = ["IA", "Médecine", "Physique", "Biologie", "Mathématiques"];

  const submit = async () => {
    try {
      await fetch("http://localhost:5000/enable_notifications", {
        method: "POST",
        headers: { 
          "Content-Type": "application/json",
          Authorization: `Bearer ${token}`
        },
        body: JSON.stringify({ enabled, interests }),
      });
      onClose();
    } catch (err) {
      alert("Erreur lors de la sauvegarde des préférences.");
    }
  };

  // Quand l'utilisateur refuse les notifications
  if (enabled === false) {
    onClose();
    return null;
  }

  return (
    <div
      className="fixed inset-0 bg-black bg-opacity-40 flex items-center justify-center z-50"
      role="dialog"
      aria-modal="true"
      aria-labelledby="modal-title"
    >
      <div className="bg-white rounded-xl shadow-lg max-w-md w-full p-6 text-center">
        {enabled === null ? (
          <>
            <h2 id="modal-title" className="mb-4 text-xl font-bold">
              Notifications
            </h2>
            <p className="mb-6">
              Voulez-vous recevoir des notifications sur les nouveaux articles ?
            </p>
            <div className="flex justify-center space-x-4">
              <button
                onClick={() => setEnabled(true)}
                className="px-5 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition"
              >
                Oui
              </button>
              <button
                onClick={() => setEnabled(false)}
                className="px-5 py-2 bg-red-600 text-white rounded-md hover:bg-red-700 transition"
              >
                Non
              </button>
            </div>
          </>
        ) : (
          <>
            <h2 id="modal-title" className="mb-4 text-xl font-bold">
              Choisissez vos centres d’intérêt
            </h2>
            <form
              onSubmit={(e) => {
                e.preventDefault();
                submit();
              }}
              className="text-left max-h-72 overflow-y-auto"
            >
              {allInterests.map((interest) => (
                <label
                  key={interest}
                  className="block mb-3 cursor-pointer select-none"
                >
                  <input
                    type="checkbox"
                    value={interest}
                    className="mr-2"
                    onChange={(e) => {
                      const checked = e.target.checked;
                      setInterests((prev) =>
                        checked ? [...prev, interest] : prev.filter((i) => i !== interest)
                      );
                    }}
                  />
                  {interest}
                </label>
              ))}

              <div className="mt-6 text-center">
                <button
                  type="submit"
                  className="px-6 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 transition font-semibold"
                >
                  Valider
                </button>
              </div>
            </form>
          </>
        )}
      </div>
    </div>
  );
}

export default NotificationConsentModal;
