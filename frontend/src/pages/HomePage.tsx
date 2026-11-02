export default function HomePage() {
  return (
    <div className="p-8 text-center">
      <h1 className="text-3xl font-bold mb-4">Pharmacogenomics Pipeline</h1>
      <p className="mb-6">Analyze genetic variants for drug interaction predictions.</p>
      <a href="/login" className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700">
        Sign in with Google
      </a>
    </div>
  );
}
