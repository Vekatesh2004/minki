import { useAuth } from "../../contexts/AuthContext";

export default function LoginPage() {
  const { login, isLoading } = useAuth();

  return (
    <div className="p-8">
      <h1 className="text-2xl font-bold mb-6">Sign In</h1>
      <button
        onClick={login}
        disabled={isLoading}
        className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-50"
      >
        {isLoading ? "Signing in..." : "Sign in with Google"}
      </button>
    </div>
  );
}
