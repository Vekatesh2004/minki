"""
Pydantic models for the Pharmacogenomics Pipeline API
"""

from datetime import datetime
from typing import List, Optional, Dict, Any, Union
from enum import Enum

from pydantic import BaseModel, Field, EmailStr, validator


class AnalysisStatus(str, Enum):
    """Analysis status enumeration"""
    UPLOADED = "uploaded"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class AnalysisType(str, Enum):
    """Analysis type enumeration"""
    BASIC = "basic"
    FULL = "full"
    CUSTOM = "custom"
    UPLOAD_ONLY = "upload_only"


# User models
class UserBase(BaseModel):
    """Base user model"""
    email: EmailStr
    full_name: str
    is_active: bool = True


class UserCreate(UserBase):
    """User creation model"""
    password: str = Field(..., min_length=8)
    
    @validator('password')
    def validate_password(cls, v):
        if len(v) < 8:
            raise ValueError('Password must be at least 8 characters long')
        return v


class UserLogin(BaseModel):
    """User login model"""
    email: EmailStr
    password: str


class UserModel(UserBase):
    """User model with ID"""
    id: int
    created_at: datetime
    last_login: Optional[datetime] = None
    
    class Config:
        from_attributes = True


class UserResponse(BaseModel):
    """User response model"""
    id: int
    email: str
    full_name: str
    is_active: bool
    created_at: datetime
    last_login: Optional[datetime] = None
    
    @classmethod
    def from_model(cls, user: UserModel):
        return cls(**user.dict())


class TokenResponse(BaseModel):
    """Token response model"""
    access_token: str
    token_type: str
    expires_in: int = 3600


# Analysis models
class AnalysisCreate(BaseModel):
    """Analysis creation model"""
    upload_id: str
    sample_id: str
    filename: str
    file_path: str
    analysis_type: AnalysisType
    user_id: int
    description: Optional[str] = None


class AnalysisModel(BaseModel):
    """Analysis model"""
    id: int
    upload_id: str
    sample_id: str
    filename: str
    file_path: str
    analysis_type: AnalysisType
    status: AnalysisStatus
    user_id: int
    description: Optional[str] = None
    results: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    report_html_path: Optional[str] = None
    report_json_path: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    
    class Config:
        from_attributes = True


class AnalysisResponse(BaseModel):
    """Analysis response model"""
    id: int
    upload_id: str
    sample_id: str
    filename: str
    analysis_type: str
    status: str
    description: Optional[str] = None
    progress: Optional[float] = None
    results_available: bool = False
    error_message: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    duration_seconds: Optional[float] = None
    
    @classmethod
    def from_model(cls, analysis: AnalysisModel):
        duration = None
        if analysis.started_at and analysis.completed_at:
            duration = (analysis.completed_at - analysis.started_at).total_seconds()
        
        return cls(
            id=analysis.id,
            upload_id=analysis.upload_id,
            sample_id=analysis.sample_id,
            filename=analysis.filename,
            analysis_type=analysis.analysis_type.value,
            status=analysis.status.value,
            description=analysis.description,
            results_available=analysis.status == AnalysisStatus.COMPLETED and analysis.results is not None,
            error_message=analysis.error_message,
            created_at=analysis.created_at,
            started_at=analysis.started_at,
            completed_at=analysis.completed_at,
            duration_seconds=duration
        )


# File upload models
class FileUploadResponse(BaseModel):
    """File upload response"""
    upload_id: str
    analysis_id: int
    message: str
    queued_for_analysis: bool = False


# Pipeline-specific models
class VCFParseRequest(BaseModel):
    """VCF parsing request"""
    file_path: str
    sample_id: str


class VCFParseResponse(BaseModel):
    """VCF parsing response"""
    sample_id: str
    total_variants: int
    qc_summary: Dict[str, Any]
    variants: List[Dict[str, Any]]


class DrugMatchRequest(BaseModel):
    """Drug matching request"""
    gene_symbol: str
    uniprot_id: Optional[str] = None
    hgvsp: Optional[str] = None


class DrugMatchResponse(BaseModel):
    """Drug matching response"""
    gene_symbol: str
    matched_drugs: List[Dict[str, Any]]
    total_matches: int
    pharmgkb_annotations: List[Dict[str, Any]]


class ProteinStructureRequest(BaseModel):
    """Protein structure analysis request"""
    uniprot_id: str
    residue_position: int
    amino_acid_change: Optional[str] = None


class ProteinStructureResponse(BaseModel):
    """Protein structure analysis response"""
    uniprot_id: str
    residue_position: int
    amino_acid_change: Optional[str] = None
    structure_analysis: Dict[str, Any]
    visualization_available: bool = False


# Report models
class ReportGenerationRequest(BaseModel):
    """Report generation request"""
    sample_id: str
    analysis_results: Dict[str, Any]
    format: str = Field(default="html", regex="^(html|json|pdf)$")


class ReportGenerationResponse(BaseModel):
    """Report generation response"""
    report_id: str
    format: str
    file_path: str
    download_url: str


# Batch analysis models
class BatchAnalysisRequest(BaseModel):
    """Batch analysis request"""
    sample_ids: List[str]
    analysis_type: AnalysisType = AnalysisType.BASIC
    description: Optional[str] = None


class BatchAnalysisResponse(BaseModel):
    """Batch analysis response"""
    batch_id: str
    analysis_ids: List[int]
    total_samples: int
    queued_count: int


# System models
class SystemStatus(BaseModel):
    """System status model"""
    status: str
    version: str
    uptime: str
    components: Dict[str, str]
    active_analyses: int
    queue_size: int


class SystemStats(BaseModel):
    """System statistics model"""
    total_analyses: int
    active_analyses: int
    completed_analyses: int
    failed_analyses: int
    total_users: int
    avg_analysis_time: float
    system_load: Dict[str, float]


# Configuration models
class DatabaseConfig(BaseModel):
    """Database configuration"""
    host: str = "localhost"
    port: int = 5432
    database: str = "pharmacogenomics"
    username: str
    password: str
    pool_size: int = 10
    max_overflow: int = 20


class RedisConfig(BaseModel):
    """Redis configuration"""
    host: str = "localhost"
    port: int = 6379
    db: int = 0
    password: Optional[str] = None


class AppConfig(BaseModel):
    """Application configuration"""
    debug: bool = False
    secret_key: str
    access_token_expire_minutes: int = 1440  # 24 hours
    upload_max_size: int = 50 * 1024 * 1024  # 50MB
    allowed_file_types: List[str] = [".vcf", ".vcf.gz"]
    cors_origins: List[str] = ["*"]


# Pagination models
class PaginatedResponse(BaseModel):
    """Paginated response model"""
    items: List[Any]
    total: int
    page: int
    size: int
    pages: int


class PaginationParams(BaseModel):
    """Pagination parameters"""
    page: int = Field(1, ge=1)
    size: int = Field(20, ge=1, le=100)
    
    @property
    def skip(self) -> int:
        return (self.page - 1) * self.size


# WebSocket models
class WebSocketMessage(BaseModel):
    """WebSocket message model"""
    type: str
    data: Dict[str, Any]
    timestamp: datetime = Field(default_factory=datetime.now)


class AnalysisProgressUpdate(BaseModel):
    """Analysis progress update"""
    analysis_id: int
    status: AnalysisStatus
    progress: float = Field(ge=0.0, le=100.0)
    current_step: str
    message: Optional[str] = None
    estimated_remaining: Optional[int] = None  # seconds


# Error models
class APIError(BaseModel):
    """API error response"""
    error: str
    message: str
    details: Optional[Dict[str, Any]] = None
    timestamp: datetime = Field(default_factory=datetime.now)


class ValidationError(BaseModel):
    """Validation error model"""
    field: str
    message: str
    value: Any


# Export models for external integrations
class ExternalAnalysisRequest(BaseModel):
    """External analysis request model"""
    vcf_data: str  # Base64 encoded VCF content
    sample_id: str
    analysis_type: AnalysisType = AnalysisType.BASIC
    callback_url: Optional[str] = None
    api_key: str


class ExternalAnalysisResponse(BaseModel):
    """External analysis response model"""
    request_id: str
    status: str
    results_url: Optional[str] = None
    estimated_completion: Optional[datetime] = None