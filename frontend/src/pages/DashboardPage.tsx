import { useAuth } from "../contexts/AuthContext";

export default function DashboardPage() {
  const { user, logout } = useAuth();

  return (
    <div className="p-8">
      <h1 className="text-2xl font-bold mb-6">Dashboard</h1>
      {user && (
        <div className="mb-6">
          <p>Signed in as: {user.email || user.displayName || "Unknown"}</p>
          <img src={user.photoURL} alt="profile" className="w-16 h-16 rounded-full" />
        </div>
      )}
      <button onClick={logout} className="px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700">
        Sign Out
      </button>
    </div>
  );
}
