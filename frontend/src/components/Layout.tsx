import { Outlet } from "react-router-dom";

export default function Layout() {
  return (
    <div className="p-4">
      <nav className="mb-4 flex gap-4">
        <a href="/" className="text-blue-600 hover:underline">Home</a>
        <a href="/demo" className="text-blue-600 hover:underline">Demo</a>
        <a href="/login" className="text-blue-600 hover:underline">Login</a>
        <a href="/dashboard" className="text-blue-600 hover:underline">Dashboard</a>
      </nav>
      <Outlet />
    </div>
  );
}
