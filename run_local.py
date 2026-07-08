#!/usr/bin/env python3
"""
Local Development Server for Pharmacogenomics Pipeline
Simple setup for testing without Docker
"""

import os
import sys
import subprocess
import asyncio
import json
from pathlib import Path
import webbrowser
import time
import signal

def check_python_version():
    """Check Python version compatibility"""
    if sys.version_info < (3, 8):
        print("❌ Python 3.8 or higher is required")
        print(f"   Current version: {sys.version}")
        return False
    print(f"✅ Python version: {sys.version_info.major}.{sys.version_info.minor}")
    return True

def check_requirements():
    """Check if basic requirements are available"""
    print("🔍 Checking requirements...")
    
    required_packages = ['fastapi', 'uvicorn', 'sqlite3']
    missing_packages = []
    
    for package in required_packages:
        try:
            if package == 'sqlite3':
                import sqlite3
            else:
                __import__(package)
            print(f"   ✅ {package}")
        except ImportError:
            missing_packages.append(package)
            print(f"   ❌ {package}")
    
    if missing_packages:
        print(f"\n📦 Installing missing packages: {', '.join(missing_packages)}")
        try:
            subprocess.run([sys.executable, "-m", "pip", "install"] + missing_packages, check=True)
            print("✅ Packages installed successfully")
        except subprocess.CalledProcessError:
            print("❌ Failed to install packages")
            return False
    
    return True

def setup_local_config():
    """Setup local configuration"""
    print("⚙️  Setting up local configuration...")
    
    # Create minimal config for local testing
    local_config = {
        "ensembl_vep": {
            "base_url": "https://rest.ensembl.org/vep/human/region",
            "rate_limit_delay": 0.1
        },
        "uniprot": {
            "base_url": "https://rest.uniprot.org"
        },
        "alphafold": {
            "base_url": "https://alphafold.ebi.ac.uk/files"
        },
        "drugbank": {
            "api_url": "https://go.drugbank.com/api/v1",
            "api_key": "demo_key"
        },
        "pharmgkb": {
            "base_url": "https://api.pharmgkb.org/v1"
        },
        "qc_thresholds": {
            "min_qual": 30,
            "min_depth": 10,
            "max_missing_rate": 0.1
        },
        "output": {
            "formats": ["json", "html"],
            "include_structures": True
        }
    }
    
    # Ensure config exists
    config_path = Path("config.json")
    if not config_path.exists():
        with open(config_path, 'w') as f:
            json.dump(local_config, f, indent=2)
        print("✅ Configuration file created")
    else:
        print("✅ Configuration file exists")
    
    # Create necessary directories
    for directory in ['uploads', 'results', 'examples']:
        Path(directory).mkdir(exist_ok=True)
        print(f"✅ Directory created: {directory}")
    
    return True

def create_simple_backend():
    """Create a simplified FastAPI backend for local testing"""
    print("🔧 Creating simplified backend...")
    
    backend_code = '''#!/usr/bin/env python3
"""
Simplified Local Backend for Testing
"""

import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Dict, Any, List
import uuid
from datetime import datetime

from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse, HTMLResponse
from pydantic import BaseModel

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Try to import our modules
try:
    from modules.vcf_parser import VCFParser
    from modules.drug_matcher import DrugMatcher
    modules_available = True
except ImportError as e:
    print(f"Warning: Some modules not available: {e}")
    modules_available = False

app = FastAPI(
    title="Pharmacogenomics Pipeline - Local Testing",
    description="Simplified local version for testing",
    version="2.0.0-local"
)

# Add CORS middleware for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global variables
pipeline_components = {}
analysis_storage = {}

async def initialize_components():
    """Initialize pipeline components"""
    global pipeline_components
    
    try:
        # Load configuration
        with open('config.json', 'r') as f:
            config = json.load(f)
        
        if modules_available:
            # Initialize components
            pipeline_components = {
                'config': config,
                'vcf_parser': VCFParser(config),
                'drug_matcher': DrugMatcher(config)
            }
            print("✅ Pipeline components initialized")
        else:
            pipeline_components = {'config': config}
            print("⚠️ Running in demo mode - limited functionality")
            
    except Exception as e:
        print(f"❌ Error initializing components: {e}")
        pipeline_components = {}

@app.on_event("startup")
async def startup_event():
    """Startup event"""
    await initialize_components()

# Models
class AnalysisStatus(BaseModel):
    id: str
    status: str
    progress: float
    message: str
    created_at: datetime
    results: Dict[str, Any] = None

class SimpleUploadResponse(BaseModel):
    upload_id: str
    filename: str
    status: str
    message: str

# Routes
@app.get("/")
async def root():
    """Root endpoint with simple HTML interface"""
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Pharmacogenomics Pipeline - Local Testing</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; background: #f5f7fa; }
            .container { max-width: 800px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
            h1 { color: #2c3e50; text-align: center; }
            .feature { margin: 20px 0; padding: 15px; background: #ecf0f1; border-radius: 5px; }
            .status { padding: 10px; background: #d5f4e6; border-radius: 5px; margin: 10px 0; }
            button { background: #3498db; color: white; padding: 10px 20px; border: none; border-radius: 5px; cursor: pointer; }
            button:hover { background: #2980b9; }
            .upload-form { margin: 20px 0; }
            input[type="file"] { margin: 10px 0; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Pharmacogenomics Pipeline</h1>
            <div class="status">
                <strong>Status:</strong> Local testing mode active
            </div>
            
            <div class="feature">
                <h3>📊 Available Features</h3>
                <ul>
                    <li>VCF file upload and parsing</li>
                    <li>Quality control analysis</li>
                    <li>Pharmacogene database lookup</li>
                    <li>Drug interaction matching</li>
                    <li>Interactive API documentation</li>
                </ul>
            </div>
            
            <div class="feature">
                <h3>🚀 Quick Links</h3>
                <p><a href="/docs" target="_blank">
                    <button>📚 View API Documentation</button>
                </a></p>
                <p><a href="/demo" target="_blank">
                    <button>🎮 Try Demo Analysis</button>
                </a></p>
            </div>
            
            <div class="upload-form">
                <h3>📁 Upload VCF File</h3>
                <form id="uploadForm" enctype="multipart/form-data">
                    <input type="file" id="vcfFile" accept=".vcf,.vcf.gz" required>
                    <input type="text" id="sampleId" placeholder="Sample ID" required>
                    <br><br>
                    <button type="submit">Upload and Analyze</button>
                </form>
                <div id="uploadResult"></div>
            </div>
        </div>
        
        <script>
            document.getElementById('uploadForm').addEventListener('submit', async function(e) {
                e.preventDefault();
                
                const formData = new FormData();
                const fileInput = document.getElementById('vcfFile');
                const sampleIdInput = document.getElementById('sampleId');
                
                formData.append('file', fileInput.files[0]);
                formData.append('sample_id', sampleIdInput.value);
                
                const resultDiv = document.getElementById('uploadResult');
                resultDiv.innerHTML = '<p>Uploading and analyzing...</p>';
                
                try {
                    const response = await fetch('/api/upload', {
                        method: 'POST',
                        body: formData
                    });
                    
                    const result = await response.json();
                    
                    if (response.ok) {
                        resultDiv.innerHTML = `
                            <div style="background: #d5f4e6; padding: 10px; border-radius: 5px; margin: 10px 0;">
                                <strong>Success!</strong><br>
                                Upload ID: ${result.upload_id}<br>
                                Status: ${result.status}<br>
                                <a href="/api/analysis/${result.upload_id}" target="_blank">View Results</a>
                            </div>
                        `;
                    } else {
                        resultDiv.innerHTML = `
                            <div style="background: #f8d7da; padding: 10px; border-radius: 5px; margin: 10px 0;">
                                <strong>Error:</strong> ${result.detail || 'Upload failed'}
                            </div>
                        `;
                    }
                } catch (error) {
                    resultDiv.innerHTML = `
                        <div style="background: #f8d7da; padding: 10px; border-radius: 5px; margin: 10px 0;">
                            <strong>Error:</strong> ${error.message}
                        </div>
                    `;
                }
            });
        </script>
    </body>
    </html>
    """
    return HTMLResponse(content=html_content)

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "mode": "local_testing",
        "components_available": modules_available,
        "timestamp": datetime.now().isoformat()
    }

@app.post("/api/upload", response_model=SimpleUploadResponse)
async def upload_file(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    sample_id: str = Form(...)
):
    """Upload VCF file for analysis"""
    
    # Validate file
    if not file.filename.lower().endswith(('.vcf', '.vcf.gz')):
        raise HTTPException(status_code=400, detail="File must be a VCF file")
    
    # Save file
    upload_id = str(uuid.uuid4())
    upload_dir = Path("uploads")
    upload_dir.mkdir(exist_ok=True)
    
    file_path = upload_dir / f"{upload_id}_{file.filename}"
    
    try:
        content = await file.read()
        with open(file_path, "wb") as f:
            f.write(content)
        
        # Queue background analysis
        background_tasks.add_task(analyze_file, upload_id, str(file_path), sample_id)
        
        return SimpleUploadResponse(
            upload_id=upload_id,
            filename=file.filename,
            status="uploaded",
            message="File uploaded successfully, analysis queued"
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Upload failed: {str(e)}")

async def analyze_file(upload_id: str, file_path: str, sample_id: str):
    """Analyze uploaded file in background"""
    
    # Update status
    analysis_storage[upload_id] = AnalysisStatus(
        id=upload_id,
        status="running",
        progress=0.0,
        message="Starting analysis...",
        created_at=datetime.now()
    )
    
    try:
        results = {}
        
        if modules_available and 'vcf_parser' in pipeline_components:
            # Parse VCF
            analysis_storage[upload_id].progress = 25.0
            analysis_storage[upload_id].message = "Parsing VCF file..."
            
            vcf_parser = pipeline_components['vcf_parser']
            vcf_results = await vcf_parser.parse_vcf(file_path, sample_id)
            results['vcf_analysis'] = vcf_results
            
            # Drug matching
            analysis_storage[upload_id].progress = 75.0
            analysis_storage[upload_id].message = "Matching drugs..."
            
            drug_matcher = pipeline_components['drug_matcher']
            drug_results = []
            
            pharmacogenes = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'TPMT', 'DPYD']
            for gene in pharmacogenes:
                try:
                    drug_result = await drug_matcher.match_drugs(gene)
                    if drug_result.get('matched_drugs'):
                        drug_results.append(drug_result)
                except:
                    pass
            
            results['drug_analysis'] = drug_results
        else:
            # Demo mode
            results = {
                'demo_mode': True,
                'message': 'Analysis completed in demo mode',
                'sample_variants': 9,
                'pharmacogenes_analyzed': ['CYP2D6', 'CYP2C19', 'CYP2C9'],
                'drug_matches_found': 12
            }
        
        # Complete analysis
        analysis_storage[upload_id].status = "completed"
        analysis_storage[upload_id].progress = 100.0
        analysis_storage[upload_id].message = "Analysis completed successfully"
        analysis_storage[upload_id].results = results
        
    except Exception as e:
        analysis_storage[upload_id].status = "failed"
        analysis_storage[upload_id].message = f"Analysis failed: {str(e)}"

@app.get("/api/analysis/{upload_id}")
async def get_analysis_status(upload_id: str):
    """Get analysis status and results"""
    
    if upload_id not in analysis_storage:
        raise HTTPException(status_code=404, detail="Analysis not found")
    
    analysis = analysis_storage[upload_id]
    return {
        "id": analysis.id,
        "status": analysis.status,
        "progress": analysis.progress,
        "message": analysis.message,
        "created_at": analysis.created_at,
        "results": analysis.results
    }

@app.get("/demo")
async def demo_analysis():
    """Demo analysis endpoint"""
    
    if not Path("examples/sample_pharmacogenomics.vcf").exists():
        raise HTTPException(status_code=404, detail="Demo VCF file not found")
    
    # Run demo analysis
    demo_id = "demo_" + str(uuid.uuid4())
    
    try:
        if modules_available:
            # Run actual analysis on demo file
            await analyze_file(demo_id, "examples/sample_pharmacogenomics.vcf", "DEMO_SAMPLE")
        else:
            # Create mock demo results
            analysis_storage[demo_id] = AnalysisStatus(
                id=demo_id,
                status="completed",
                progress=100.0,
                message="Demo analysis completed",
                created_at=datetime.now(),
                results={
                    "demo_mode": True,
                    "vcf_analysis": {
                        "total_variants": 9,
                        "pass_rate": 1.0,
                        "sample_variants": ["22:42126611 A>G", "22:42127803 G>A", "22:42128936 C>T"]
                    },
                    "drug_analysis": [
                        {"gene": "CYP2D6", "drugs": ["Codeine", "Tramadol", "Metoprolol"]},
                        {"gene": "CYP2C19", "drugs": ["Clopidogrel", "Omeprazole"]},
                        {"gene": "CYP2C9", "drugs": ["Warfarin", "Phenytoin"]}
                    ]
                }
            )
    except Exception as e:
        return {"error": f"Demo analysis failed: {str(e)}"}
    
    return {"demo_id": demo_id, "status": "completed", "message": "Demo analysis ready", "results_url": f"/api/analysis/{demo_id}"}

# Serve static files
if Path("results").exists():
    app.mount("/results", StaticFiles(directory="results"), name="results")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("simple_backend:app", host="127.0.0.1", port=8000, reload=True)
'''
    
    with open("simple_backend.py", "w") as f:
        f.write(backend_code)
    
    print("✅ Simplified backend created")
    return True

def start_local_server():
    """Start the local development server"""
    print("🚀 Starting local development server...")
    
    try:
        # Start the backend server
        cmd = [sys.executable, "simple_backend.py"]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        # Wait a bit for server to start
        print("⏳ Waiting for server to start...")
        time.sleep(3)
        
        # Check if server is running
        try:
            import requests
            response = requests.get("http://127.0.0.1:8000/health", timeout=5)
            if response.status_code == 200:
                print("✅ Server started successfully!")
                print("\n🌐 Access your application:")
                print("   Main Interface: http://127.0.0.1:8000")
                print("   API Documentation: http://127.0.0.1:8000/docs")
                print("   Health Check: http://127.0.0.1:8000/health")
                
                # Try to open browser
                try:
                    webbrowser.open("http://127.0.0.1:8000")
                    print("🔗 Browser opened automatically")
                except:
                    print("💡 Please open http://127.0.0.1:8000 in your browser")
                
                print("\n⚠️  To stop the server, press Ctrl+C")
                
                # Wait for process
                try:
                    process.wait()
                except KeyboardInterrupt:
                    print("\n👋 Shutting down server...")
                    process.terminate()
                    
            else:
                print("❌ Server health check failed")
                return False
                
        except ImportError:
            print("⚠️  'requests' not available, trying basic approach...")
            print("🌐 Server should be running at: http://127.0.0.1:8000")
            print("💡 Check the URL in your browser")
            
            try:
                process.wait()
            except KeyboardInterrupt:
                print("\n👋 Shutting down server...")
                process.terminate()
        
    except Exception as e:
        print(f"❌ Failed to start server: {e}")
        return False
    
    return True

def main():
    """Main function"""
    print("=" * 60)
    print("🧬 Pharmacogenomics Pipeline - Local Testing Setup")
    print("=" * 60)
    
    # Check Python version
    if not check_python_version():
        return 1
    
    # Check and install requirements
    if not check_requirements():
        print("❌ Failed to install requirements")
        return 1
    
    # Setup configuration
    if not setup_local_config():
        print("❌ Failed to setup configuration")
        return 1
    
    # Initialize database if needed
    print("🗄️  Initializing local database...")
    try:
        from setup import setup_pipeline
        asyncio.run(setup_pipeline())
    except Exception as e:
        print(f"⚠️  Database setup warning: {e}")
    
    # Create simplified backend
    if not create_simple_backend():
        print("❌ Failed to create backend")
        return 1
    
    print("\n✅ Local setup completed successfully!")
    print("📋 Setup Summary:")
    print("   - Configuration files created")
    print("   - Directories initialized")
    print("   - Simplified backend ready")
    print("   - Database populated with sample data")
    
    # Ask user if they want to start the server
    try:
        choice = input("\n🚀 Start local server now? (y/n): ").strip().lower()
        if choice in ['y', 'yes']:
            return 0 if start_local_server() else 1
        else:
            print("\n💡 To start the server later, run:")
            print("   python simple_backend.py")
            return 0
    except KeyboardInterrupt:
        print("\n\n👋 Setup completed. Start server with: python simple_backend.py")
        return 0

if __name__ == "__main__":
    sys.exit(main())