#!/usr/bin/env python3
"""
Easy launcher for the Pharmacogenomics Streamlit Web Interface
"""

import os
import sys
import subprocess
from pathlib import Path


def check_requirements():
    """Check if basic requirements are met."""
    print("🔍 Checking requirements...")

    if not Path("config.json").exists():
        print("❌ config.json not found. Please run from the project directory.")
        return False

    if not Path("streamlit_app.py").exists():
        print("❌ streamlit_app.py not found. Please run from the project directory.")
        return False

    if not Path("modules").exists():
        print("❌ modules directory not found.")
        return False

    print("✅ Basic requirements met")
    return True


def install_dependencies():
    """Check if streamlit is installed."""
    try:
        import streamlit  # noqa: F401
        print("✅ Streamlit is installed")
        return True
    except ImportError:
        print("⚠️  Streamlit is not installed")
        response = input("Would you like to install dependencies now? (y/n): ").strip().lower()
        if response in ["y", "yes"]:
            print("📦 Installing dependencies...")
            try:
                subprocess.run(
                    [sys.executable, "-m", "pip", "install", "-r", "requirements_streamlit.txt"],
                    check=True,
                )
                print("✅ Dependencies installed successfully")
                return True
            except subprocess.CalledProcessError:
                print("❌ Failed to install dependencies")
                print("Please install manually: pip install -r requirements_streamlit.txt")
                return False
        else:
            print("Please install dependencies manually: pip install -r requirements_streamlit.txt")
            return False


def start_streamlit_interface():
    """Start the Streamlit web interface."""
    print("\n🚀 Starting Pharmacogenomics Streamlit Interface...")
    print("📊 The interface will be available at: http://localhost:8501")
    print("💡 To stop the server, press Ctrl+C")
    print("=" * 60)

    try:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "streamlit",
                "run",
                "streamlit_app.py",
                "--server.port=8501",
                "--server.address=0.0.0.0",
            ],
            check=True,
        )
    except KeyboardInterrupt:
        print("\n👋 Streamlit interface stopped by user")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error starting Streamlit interface: {e}")
        return False


def main():
    """Main function."""
    print("Pharmacogenomics Pipeline - Streamlit Interface Launcher")
    print("=" * 60)

    if not check_requirements():
        print("\n❌ Requirements not met. Please check the setup.")
        sys.exit(1)

    if not install_dependencies():
        print("\n❌ Dependencies not available. Cannot start Streamlit interface.")
        sys.exit(1)

    success = start_streamlit_interface()

    if success:
        print("\n✅ Streamlit interface session completed")
    else:
        print("\n❌ Failed to start Streamlit interface")
        sys.exit(1)


if __name__ == "__main__":
    main()
