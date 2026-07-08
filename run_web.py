#!/usr/bin/env python3
"""
Easy launcher for the Pharmacogenomics Web Interface
"""

import os
import sys
import subprocess
from pathlib import Path

def check_requirements():
    """Check if basic requirements are met"""
    
    print("🔍 Checking requirements...")
    
    # Check if we're in the right directory
    if not Path("config.json").exists():
        print("❌ config.json not found. Please run from the project directory.")
        return False
        
    if not Path("web_app.py").exists():
        print("❌ web_app.py not found. Please run from the project directory.")
        return False
        
    # Check if modules directory exists
    if not Path("modules").exists():
        print("❌ modules directory not found.")
        return False
        
    print("✅ Basic requirements met")
    return True

def install_dependencies():
    """Ask user if they want to install dependencies"""
    
    try:
        import flask
        print("✅ Flask is installed")
        return True
    except ImportError:
        print("⚠️  Flask is not installed")
        
        response = input("Would you like to install dependencies now? (y/n): ").strip().lower()
        if response in ['y', 'yes']:
            print("📦 Installing dependencies...")
            try:
                subprocess.run([sys.executable, "-m", "pip", "install", "flask", "flask-cors"], 
                             check=True)
                print("✅ Dependencies installed successfully")
                return True
            except subprocess.CalledProcessError:
                print("❌ Failed to install dependencies")
                print("Please install manually: pip install flask flask-cors")
                return False
        else:
            print("Please install dependencies manually: pip install flask flask-cors")
            return False

def start_web_interface():
    """Start the web interface"""
    
    print("\n🚀 Starting Pharmacogenomics Web Interface...")
    print("📊 The interface will be available at: http://localhost:5000")
    print("🌐 You can also access it from other devices at: http://YOUR_IP:5000")
    print("\n💡 To stop the server, press Ctrl+C")
    print("=" * 60)
    
    try:
        # Import and run the web app
        import web_app
        # This will start the Flask development server
        
    except ImportError as e:
        print(f"❌ Error importing web app: {e}")
        print("Make sure all dependencies are installed")
        return False
    except KeyboardInterrupt:
        print("\n👋 Web interface stopped by user")
        return True
    except Exception as e:
        print(f"❌ Error starting web interface: {e}")
        return False

def main():
    """Main function"""
    
    print("Pharmacogenomics Pipeline - Web Interface Launcher")
    print("=" * 55)
    
    # Check basic requirements
    if not check_requirements():
        print("\n❌ Requirements not met. Please check the setup.")
        sys.exit(1)
    
    # Check and install dependencies
    if not install_dependencies():
        print("\n❌ Dependencies not available. Cannot start web interface.")
        sys.exit(1)
    
    # Start the web interface
    success = start_web_interface()
    
    if success:
        print("\n✅ Web interface session completed")
    else:
        print("\n❌ Failed to start web interface")
        sys.exit(1)

if __name__ == "__main__":
    main()