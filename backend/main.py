#!/usr/bin/env python3
"""
Production FastAPI Backend for Pharmacogenomics Pipeline
"""

import asyncio
import logging
import os
import sys
from contextlib import asynccontextmanager
from typing import Dict, List, Optional
import uuid
from datetime import datetime, timedelta
import json

from fastapi import FastAPI, HTTPException, Depends, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel, Field
import uvicorn
from prometheus_client import Counter, Histogram, generate_latest
import structlog

# Add parent directory to path for module imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.vcf_parser import VCFParser
from modules.vep_annotator import VEPAnnotator
from modules.protein_structure import ProteinStructureAnalyzer
from modules.drug_matcher import DrugMatcher
from modules.report_generator import ReportGenerator
from backend.models import *
from backend.database import DatabaseManager
from backend.auth import AuthManager
from backend.task_manager import TaskManager

# Configure structured logging
structlog.configure(
    processors=[
        structlog.stdlib.filter_by_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
        structlog.processors.UnicodeDecoder(),
        structlog.processors.JSONRenderer()
    ],
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    wrapper_class=structlog.stdlib.BoundLogger,
    cache_logger_on_first_use=True,
)

logger = structlog.get_logger()

# Metrics
request_count = Counter('http_requests_total', 'Total HTTP requests', ['method', 'endpoint'])
request_duration = Histogram('http_request_duration_seconds', 'HTTP request duration')
analysis_duration = Histogram('analysis_duration_seconds', 'Analysis duration', ['analysis_type'])
upload_size = Histogram('upload_file_size_bytes', 'Uploaded file sizes')

# Global components
pipeline_components = {}
db_manager = None
auth_manager = None
task_manager = None

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager"""
    logger.info("Starting Pharmacogenomics Pipeline API")
    
    # Initialize components
    await initialize_components()
    
    yield
    
    # Cleanup
    logger.info("Shutting down Pharmacogenomics Pipeline API")
    if db_manager:
        await db_manager.close()

async def initialize_components():
    """Initialize all pipeline components"""
    global pipeline_components, db_manager, auth_manager, task_manager
    
    try:
        # Load configuration
        config_path = os.path.join(os.path.dirname(__file__), '..', 'config.json')
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        # Initialize database
        db_manager = DatabaseManager()
        await db_manager.initialize()
        
        # Initialize auth manager
        auth_manager = AuthManager(db_manager)
        
        # Initialize task manager
        task_manager = TaskManager()
        
        # Initialize pipeline components
        pipeline_components = {
            'config': config,
            'vcf_parser': VCFParser(config),
            'vep_annotator': VEPAnnotator(config),
            'structure_analyzer': ProteinStructureAnalyzer(config),
            'drug_matcher': DrugMatcher(config),
            'report_generator': ReportGenerator(config)
        }
        
        logger.info("All components initialized successfully")
        
    except Exception as e:
        logger.error("Failed to initialize components", error=str(e))
        raise

# Create FastAPI app
app = FastAPI(
    title="Pharmacogenomics Analysis Pipeline",
    description="Production API for pharmacogenomics variant analysis and drug interaction prediction",
    version="2.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
    lifespan=lifespan
)

# Add middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.add_middleware(GZipMiddleware, minimum_size=1000)

# Security
security = HTTPBearer()

async def get_current_user(credentials: HTTPAuthorizationCredentials = Depends(security)) -> UserModel:
    """Get current authenticated user"""
    if not auth_manager:
        raise HTTPException(status_code=500, detail="Authentication not initialized")
    
    user = await auth_manager.verify_token(credentials.credentials)
    if not user:
        raise HTTPException(status_code=401, detail="Invalid authentication token")
    
    return user

# Middleware for metrics
@app.middleware("http")
async def metrics_middleware(request, call_next):
    start_time = datetime.now()
    
    response = await call_next(request)
    
    # Record metrics
    duration = (datetime.now() - start_time).total_seconds()
    request_count.labels(method=request.method, endpoint=request.url.path).inc()
    request_duration.observe(duration)
    
    return response

# Health check endpoints
@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "version": "2.0.0",
        "components": {
            "database": "connected" if db_manager and db_manager.is_connected() else "disconnected",
            "pipeline": "initialized" if pipeline_components else "not_initialized"
        }
    }

@app.get("/metrics")
async def get_metrics():
    """Prometheus metrics endpoint"""
    return Response(generate_latest(), media_type="text/plain")

# Authentication endpoints
@app.post("/api/auth/register", response_model=UserResponse)
async def register_user(user_data: UserCreate):
    """Register a new user"""
    if not auth_manager:
        raise HTTPException(status_code=500, detail="Authentication not available")
    
    try:
        user = await auth_manager.create_user(user_data)
        return UserResponse.from_model(user)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/auth/login", response_model=TokenResponse)
async def login_user(login_data: UserLogin):
    """Login user and return JWT token"""
    if not auth_manager:
        raise HTTPException(status_code=500, detail="Authentication not available")
    
    try:
        token = await auth_manager.authenticate_user(login_data.email, login_data.password)
        return TokenResponse(access_token=token, token_type="bearer")
    except ValueError as e:
        raise HTTPException(status_code=401, detail=str(e))

@app.get("/api/auth/me", response_model=UserResponse)
async def get_current_user_info(current_user: UserModel = Depends(get_current_user)):
    """Get current user information"""
    return UserResponse.from_model(current_user)

# File upload endpoints
@app.post("/api/upload", response_model=FileUploadResponse)
async def upload_vcf_file(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    sample_id: str = Form(...),
    analysis_type: str = Form("basic"),
    current_user: UserModel = Depends(get_current_user)
):
    """Upload VCF file for analysis"""
    
    # Validate file
    if not file.filename.lower().endswith(('.vcf', '.vcf.gz')):
        raise HTTPException(status_code=400, detail="File must be a VCF file (.vcf or .vcf.gz)")
    
    # Check file size (50MB limit)
    content = await file.read()
    if len(content) > 50 * 1024 * 1024:
        raise HTTPException(status_code=400, detail="File size exceeds 50MB limit")
    
    upload_size.observe(len(content))
    
    try:
        # Save file
        upload_id = str(uuid.uuid4())
        upload_dir = f"uploads/{current_user.id}"
        os.makedirs(upload_dir, exist_ok=True)
        
        file_path = f"{upload_dir}/{upload_id}_{file.filename}"
        with open(file_path, "wb") as f:
            f.write(content)
        
        # Create analysis record
        analysis_record = AnalysisCreate(
            upload_id=upload_id,
            sample_id=sample_id,
            filename=file.filename,
            file_path=file_path,
            analysis_type=analysis_type,
            user_id=current_user.id
        )
        
        # Save to database
        analysis = await db_manager.create_analysis(analysis_record)
        
        # Queue background analysis
        if analysis_type != "upload_only":
            background_tasks.add_task(run_background_analysis, analysis.id)
        
        return FileUploadResponse(
            upload_id=upload_id,
            analysis_id=analysis.id,
            message="File uploaded successfully",
            queued_for_analysis=analysis_type != "upload_only"
        )
        
    except Exception as e:
        logger.error("File upload failed", error=str(e), user_id=current_user.id)
        raise HTTPException(status_code=500, detail=f"Upload failed: {str(e)}")

# Analysis endpoints
@app.post("/api/analysis/{analysis_id}/start", response_model=AnalysisResponse)
async def start_analysis(
    analysis_id: int,
    background_tasks: BackgroundTasks,
    current_user: UserModel = Depends(get_current_user)
):
    """Start analysis for uploaded file"""
    
    # Get analysis record
    analysis = await db_manager.get_analysis(analysis_id)
    if not analysis:
        raise HTTPException(status_code=404, detail="Analysis not found")
    
    if analysis.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    if analysis.status != "uploaded":
        raise HTTPException(status_code=400, detail="Analysis already started or completed")
    
    # Update status and queue analysis
    await db_manager.update_analysis_status(analysis_id, "queued")
    background_tasks.add_task(run_background_analysis, analysis_id)
    
    updated_analysis = await db_manager.get_analysis(analysis_id)
    return AnalysisResponse.from_model(updated_analysis)

@app.get("/api/analysis/{analysis_id}", response_model=AnalysisResponse)
async def get_analysis(
    analysis_id: int,
    current_user: UserModel = Depends(get_current_user)
):
    """Get analysis status and results"""
    
    analysis = await db_manager.get_analysis(analysis_id)
    if not analysis:
        raise HTTPException(status_code=404, detail="Analysis not found")
    
    if analysis.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    return AnalysisResponse.from_model(analysis)

@app.get("/api/analysis/{analysis_id}/results")
async def get_analysis_results(
    analysis_id: int,
    format: str = "json",
    current_user: UserModel = Depends(get_current_user)
):
    """Download analysis results"""
    
    analysis = await db_manager.get_analysis(analysis_id)
    if not analysis:
        raise HTTPException(status_code=404, detail="Analysis not found")
    
    if analysis.user_id != current_user.id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    if analysis.status != "completed":
        raise HTTPException(status_code=400, detail="Analysis not completed")
    
    # Return results based on format
    if format == "json":
        return JSONResponse(content=analysis.results)
    elif format == "html" and analysis.report_html_path:
        return FileResponse(
            analysis.report_html_path,
            media_type="text/html",
            filename=f"{analysis.sample_id}_report.html"
        )
    else:
        raise HTTPException(status_code=400, detail="Invalid format or format not available")

@app.get("/api/analyses", response_model=List[AnalysisResponse])
async def list_user_analyses(
    skip: int = 0,
    limit: int = 50,
    current_user: UserModel = Depends(get_current_user)
):
    """List user's analyses"""
    
    analyses = await db_manager.get_user_analyses(current_user.id, skip, limit)
    return [AnalysisResponse.from_model(analysis) for analysis in analyses]

# Pipeline-specific endpoints
@app.post("/api/pipeline/parse-vcf", response_model=VCFParseResponse)
async def parse_vcf_endpoint(
    request: VCFParseRequest,
    current_user: UserModel = Depends(get_current_user)
):
    """Parse VCF file and return QC metrics"""
    
    try:
        parser = pipeline_components['vcf_parser']
        results = await parser.parse_vcf(request.file_path, request.sample_id)
        
        return VCFParseResponse(
            sample_id=request.sample_id,
            total_variants=len(results.get('variants', [])),
            qc_summary=results.get('qc_summary', {}),
            variants=results.get('variants', [])[:100]  # Limit for API response
        )
    except Exception as e:
        logger.error("VCF parsing failed", error=str(e))
        raise HTTPException(status_code=500, detail=f"VCF parsing failed: {str(e)}")

@app.post("/api/pipeline/match-drugs", response_model=DrugMatchResponse)
async def match_drugs_endpoint(
    request: DrugMatchRequest,
    current_user: UserModel = Depends(get_current_user)
):
    """Match drugs for specified gene"""
    
    try:
        drug_matcher = pipeline_components['drug_matcher']
        results = await drug_matcher.match_drugs(
            request.gene_symbol,
            request.uniprot_id,
            request.hgvsp
        )
        
        return DrugMatchResponse(
            gene_symbol=request.gene_symbol,
            matched_drugs=results.get('matched_drugs', []),
            total_matches=len(results.get('matched_drugs', [])),
            pharmgkb_annotations=results.get('pharmgkb_annotations', [])
        )
    except Exception as e:
        logger.error("Drug matching failed", error=str(e))
        raise HTTPException(status_code=500, detail=f"Drug matching failed: {str(e)}")

# Admin endpoints
@app.get("/api/admin/stats", dependencies=[Depends(get_current_user)])
async def get_system_stats():
    """Get system statistics (admin only)"""
    
    stats = {
        "total_analyses": await db_manager.count_analyses(),
        "active_analyses": await db_manager.count_analyses_by_status("running"),
        "completed_analyses": await db_manager.count_analyses_by_status("completed"),
        "failed_analyses": await db_manager.count_analyses_by_status("failed"),
        "total_users": await db_manager.count_users(),
        "system_info": {
            "version": "2.0.0",
            "uptime": str(datetime.now() - app.state.start_time) if hasattr(app.state, 'start_time') else "unknown"
        }
    }
    
    return stats

# Background task functions
async def run_background_analysis(analysis_id: int):
    """Run analysis in background"""
    
    start_time = datetime.now()
    
    try:
        # Get analysis record
        analysis = await db_manager.get_analysis(analysis_id)
        if not analysis:
            logger.error("Analysis not found", analysis_id=analysis_id)
            return
        
        # Update status
        await db_manager.update_analysis_status(analysis_id, "running")
        
        logger.info("Starting background analysis", analysis_id=analysis_id, sample_id=analysis.sample_id)
        
        # Run analysis based on type
        if analysis.analysis_type == "basic":
            results = await run_basic_analysis_pipeline(analysis)
        elif analysis.analysis_type == "full":
            results = await run_full_analysis_pipeline(analysis)
        else:
            raise ValueError(f"Unknown analysis type: {analysis.analysis_type}")
        
        # Save results
        await db_manager.update_analysis_results(analysis_id, results, "completed")
        
        duration = (datetime.now() - start_time).total_seconds()
        analysis_duration.labels(analysis_type=analysis.analysis_type).observe(duration)
        
        logger.info("Analysis completed", analysis_id=analysis_id, duration=duration)
        
    except Exception as e:
        # Update status to failed
        await db_manager.update_analysis_status(analysis_id, "failed")
        await db_manager.update_analysis_error(analysis_id, str(e))
        
        duration = (datetime.now() - start_time).total_seconds()
        
        logger.error("Analysis failed", 
                    analysis_id=analysis_id, 
                    error=str(e), 
                    duration=duration)

async def run_basic_analysis_pipeline(analysis):
    """Run basic analysis pipeline"""
    
    vcf_parser = pipeline_components['vcf_parser']
    drug_matcher = pipeline_components['drug_matcher']
    
    # Parse VCF
    vcf_results = await vcf_parser.parse_vcf(analysis.file_path, analysis.sample_id)
    
    # Match drugs for known pharmacogenes
    pharmacogenes = ['CYP2D6', 'CYP2C19', 'CYP2C9', 'TPMT', 'DPYD', 'VKORC1']
    drug_results = []
    
    for gene in pharmacogenes:
        try:
            drug_result = await drug_matcher.match_drugs(gene)
            if drug_result.get('matched_drugs'):
                drug_results.append(drug_result)
        except Exception as e:
            logger.warning("Drug matching failed for gene", gene=gene, error=str(e))
    
    return {
        "analysis_type": "basic",
        "vcf_results": vcf_results,
        "drug_results": drug_results,
        "summary": {
            "total_variants": len(vcf_results.get('variants', [])),
            "pharmacogenes_analyzed": len(pharmacogenes),
            "drug_matches": sum(len(r.get('matched_drugs', [])) for r in drug_results)
        }
    }

async def run_full_analysis_pipeline(analysis):
    """Run full analysis pipeline"""
    
    # This would include VEP annotation, structure analysis, etc.
    # For now, same as basic but could be extended
    return await run_basic_analysis_pipeline(analysis)

# Static file serving for production
app.mount("/static", StaticFiles(directory="static"), name="static")

if __name__ == "__main__":
    # For development
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_config={
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "default": {
                    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                },
            },
            "handlers": {
                "default": {
                    "formatter": "default",
                    "class": "logging.StreamHandler",
                    "stream": "ext://sys.stdout",
                },
            },
            "root": {
                "level": "INFO",
                "handlers": ["default"],
            },
        }
    )