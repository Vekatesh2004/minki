#!/bin/bash

# Local Installation Script for Pharmacogenomics Pipeline
# This script sets up everything you need for local testing

echo "🧬 Pharmacogenomics Pipeline - Local Setup"
echo "=========================================="

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is required but not installed"
    echo "Please install Python 3.8+ and try again"
    exit 1
fi

echo "✅ Python 3 found: $(python3 --version)"

# Create virtual environment (optional but recommended)
read -p "🐍 Create virtual environment? (recommended) [y/n]: " create_venv

if [[ $create_venv =~ ^[Yy]$ ]]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    
    # Activate virtual environment
    if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
        source venv/Scripts/activate
    else
        source venv/bin/activate
    fi
    
    echo "✅ Virtual environment created and activated"
fi

# Install minimal requirements for local testing
echo "📦 Installing minimal requirements..."

# Create a minimal requirements file for local testing
cat > requirements_local.txt << EOF
# Minimal requirements for local testing
fastapi>=0.104.0
uvicorn[standard]>=0.24.0
python-multipart>=0.0.6
requests>=2.31.0
pandas>=2.0.0
python-dotenv>=1.0.0

# Optional: for full functionality
# cyvcf2>=0.30.0
# biopython>=1.81
# aiohttp>=3.8.0
EOF

# Install requirements
pip install -r requirements_local.txt

if [ $? -eq 0 ]; then
    echo "✅ Requirements installed successfully"
else
    echo "❌ Failed to install requirements"
    echo "Trying alternative installation method..."
    
    # Try installing one by one
    pip install fastapi uvicorn python-multipart requests pandas python-dotenv
fi

# Run the local setup
echo "🚀 Running local setup..."
python3 run_local.py

# Check if setup was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Local setup completed successfully!"
    echo ""
    echo "📋 What was set up:"
    echo "   ✅ Configuration files"
    echo "   ✅ Directory structure"
    echo "   ✅ Sample data"
    echo "   ✅ Simplified web interface"
    echo ""
    echo "🚀 To start the server:"
    echo "   python3 simple_backend.py"
    echo ""
    echo "🌐 Then visit: http://localhost:8000"
    echo ""
    echo "📚 API docs: http://localhost:8000/docs"
    echo ""
else
    echo "❌ Setup encountered issues"
    echo "💡 You can try running 'python3 run_local.py' manually"
fi