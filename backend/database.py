"""
Database management for Pharmacogenomics Pipeline
"""

import asyncio
import os
from datetime import datetime
from typing import List, Optional, Dict, Any

import asyncpg
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, DateTime, Boolean, JSON, Text, ForeignKey, Float
from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession
from sqlalchemy.orm import sessionmaker, declarative_base
from sqlalchemy.sql import select, insert, update, delete
import structlog

from .models import UserCreate, UserModel, AnalysisCreate, AnalysisModel, AnalysisStatus

logger = structlog.get_logger()

Base = declarative_base()

# Database tables
class User(Base):
    __tablename__ = "users"
    
    id = Column(Integer, primary_key=True, index=True)
    email = Column(String, unique=True, index=True)
    full_name = Column(String)
    hashed_password = Column(String)
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=datetime.now)
    last_login = Column(DateTime, nullable=True)

class Analysis(Base):
    __tablename__ = "analyses"
    
    id = Column(Integer, primary_key=True, index=True)
    upload_id = Column(String, unique=True, index=True)
    sample_id = Column(String, index=True)
    filename = Column(String)
    file_path = Column(String)
    analysis_type = Column(String)
    status = Column(String, default="uploaded")
    user_id = Column(Integer, ForeignKey("users.id"))
    description = Column(Text, nullable=True)
    results = Column(JSON, nullable=True)
    error_message = Column(Text, nullable=True)
    report_html_path = Column(String, nullable=True)
    report_json_path = Column(String, nullable=True)
    progress = Column(Float, default=0.0)
    created_at = Column(DateTime, default=datetime.now)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

class APIKey(Base):
    __tablename__ = "api_keys"
    
    id = Column(Integer, primary_key=True, index=True)
    key_hash = Column(String, unique=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"))
    name = Column(String)
    is_active = Column(Boolean, default=True)
    created_at = Column(DateTime, default=datetime.now)
    last_used = Column(DateTime, nullable=True)

class SystemMetrics(Base):
    __tablename__ = "system_metrics"
    
    id = Column(Integer, primary_key=True, index=True)
    metric_name = Column(String, index=True)
    metric_value = Column(Float)
    recorded_at = Column(DateTime, default=datetime.now)

class DatabaseManager:
    """Database manager for the application"""
    
    def __init__(self):
        self.engine = None
        self.async_session = None
        self.pool = None
        self._initialized = False
    
    async def initialize(self):
        """Initialize database connection"""
        try:
            # Database configuration from environment
            db_url = os.getenv(
                "DATABASE_URL",
                "postgresql+asyncpg://postgres:password@localhost/pharmacogenomics"
            )
            
            # Create async engine
            self.engine = create_async_engine(
                db_url,
                pool_size=10,
                max_overflow=20,
                pool_timeout=30,
                echo=False
            )
            
            # Create session factory
            self.async_session = sessionmaker(
                self.engine,
                class_=AsyncSession,
                expire_on_commit=False
            )
            
            # Create tables if they don't exist
            async with self.engine.begin() as conn:
                await conn.run_sync(Base.metadata.create_all)
            
            # Create connection pool for raw queries
            self.pool = await asyncpg.create_pool(
                os.getenv("DATABASE_URL", "postgresql://postgres:password@localhost/pharmacogenomics").replace("+asyncpg", ""),
                min_size=5,
                max_size=20
            )
            
            self._initialized = True
            logger.info("Database initialized successfully")
            
        except Exception as e:
            logger.error("Failed to initialize database", error=str(e))
            raise
    
    async def close(self):
        """Close database connections"""
        if self.pool:
            await self.pool.close()
        
        if self.engine:
            await self.engine.dispose()
        
        self._initialized = False
        logger.info("Database connections closed")
    
    def is_connected(self) -> bool:
        """Check if database is connected"""
        return self._initialized and self.engine is not None
    
    # User management
    async def create_user(self, user_data: UserCreate) -> UserModel:
        """Create a new user"""
        async with self.async_session() as session:
            # Check if user already exists
            result = await session.execute(
                select(User).where(User.email == user_data.email)
            )
            existing_user = result.scalar_one_or_none()
            
            if existing_user:
                raise ValueError("User with this email already exists")
            
            # Create new user
            new_user = User(
                email=user_data.email,
                full_name=user_data.full_name,
                hashed_password=user_data.password,  # This should be hashed in auth manager
                is_active=user_data.is_active
            )
            
            session.add(new_user)
            await session.commit()
            await session.refresh(new_user)
            
            return UserModel(
                id=new_user.id,
                email=new_user.email,
                full_name=new_user.full_name,
                is_active=new_user.is_active,
                created_at=new_user.created_at,
                last_login=new_user.last_login
            )
    
    async def get_user_by_email(self, email: str) -> Optional[UserModel]:
        """Get user by email"""
        async with self.async_session() as session:
            result = await session.execute(
                select(User).where(User.email == email)
            )
            user = result.scalar_one_or_none()
            
            if not user:
                return None
            
            return UserModel(
                id=user.id,
                email=user.email,
                full_name=user.full_name,
                is_active=user.is_active,
                created_at=user.created_at,
                last_login=user.last_login
            )
    
    async def get_user_by_id(self, user_id: int) -> Optional[UserModel]:
        """Get user by ID"""
        async with self.async_session() as session:
            result = await session.execute(
                select(User).where(User.id == user_id)
            )
            user = result.scalar_one_or_none()
            
            if not user:
                return None
            
            return UserModel(
                id=user.id,
                email=user.email,
                full_name=user.full_name,
                is_active=user.is_active,
                created_at=user.created_at,
                last_login=user.last_login
            )
    
    async def update_user_login(self, user_id: int):
        """Update user's last login time"""
        async with self.async_session() as session:
            await session.execute(
                update(User)
                .where(User.id == user_id)
                .values(last_login=datetime.now())
            )
            await session.commit()
    
    # Analysis management
    async def create_analysis(self, analysis_data: AnalysisCreate) -> AnalysisModel:
        """Create a new analysis record"""
        async with self.async_session() as session:
            new_analysis = Analysis(
                upload_id=analysis_data.upload_id,
                sample_id=analysis_data.sample_id,
                filename=analysis_data.filename,
                file_path=analysis_data.file_path,
                analysis_type=analysis_data.analysis_type.value,
                user_id=analysis_data.user_id,
                description=analysis_data.description
            )
            
            session.add(new_analysis)
            await session.commit()
            await session.refresh(new_analysis)
            
            return self._analysis_to_model(new_analysis)
    
    async def get_analysis(self, analysis_id: int) -> Optional[AnalysisModel]:
        """Get analysis by ID"""
        async with self.async_session() as session:
            result = await session.execute(
                select(Analysis).where(Analysis.id == analysis_id)
            )
            analysis = result.scalar_one_or_none()
            
            if not analysis:
                return None
            
            return self._analysis_to_model(analysis)
    
    async def get_user_analyses(self, user_id: int, skip: int = 0, limit: int = 50) -> List[AnalysisModel]:
        """Get user's analyses with pagination"""
        async with self.async_session() as session:
            result = await session.execute(
                select(Analysis)
                .where(Analysis.user_id == user_id)
                .order_by(Analysis.created_at.desc())
                .offset(skip)
                .limit(limit)
            )
            analyses = result.scalars().all()
            
            return [self._analysis_to_model(analysis) for analysis in analyses]
    
    async def update_analysis_status(self, analysis_id: int, status: str):
        """Update analysis status"""
        async with self.async_session() as session:
            update_data = {"status": status}
            
            if status == "running":
                update_data["started_at"] = datetime.now()
            elif status == "completed":
                update_data["completed_at"] = datetime.now()
            
            await session.execute(
                update(Analysis)
                .where(Analysis.id == analysis_id)
                .values(**update_data)
            )
            await session.commit()
    
    async def update_analysis_progress(self, analysis_id: int, progress: float):
        """Update analysis progress"""
        async with self.async_session() as session:
            await session.execute(
                update(Analysis)
                .where(Analysis.id == analysis_id)
                .values(progress=progress)
            )
            await session.commit()
    
    async def update_analysis_results(self, analysis_id: int, results: Dict[str, Any], status: str):
        """Update analysis results"""
        async with self.async_session() as session:
            update_data = {
                "results": results,
                "status": status,
                "progress": 100.0
            }
            
            if status == "completed":
                update_data["completed_at"] = datetime.now()
            
            await session.execute(
                update(Analysis)
                .where(Analysis.id == analysis_id)
                .values(**update_data)
            )
            await session.commit()
    
    async def update_analysis_error(self, analysis_id: int, error_message: str):
        """Update analysis error"""
        async with self.async_session() as session:
            await session.execute(
                update(Analysis)
                .where(Analysis.id == analysis_id)
                .values(
                    error_message=error_message,
                    status="failed",
                    completed_at=datetime.now()
                )
            )
            await session.commit()
    
    # Statistics methods
    async def count_analyses(self) -> int:
        """Count total analyses"""
        async with self.async_session() as session:
            result = await session.execute(select(Analysis.id).count())
            return result.scalar()
    
    async def count_analyses_by_status(self, status: str) -> int:
        """Count analyses by status"""
        async with self.async_session() as session:
            result = await session.execute(
                select(Analysis.id)
                .where(Analysis.status == status)
                .count()
            )
            return result.scalar()
    
    async def count_users(self) -> int:
        """Count total users"""
        async with self.async_session() as session:
            result = await session.execute(select(User.id).count())
            return result.scalar()
    
    async def get_recent_analyses(self, limit: int = 10) -> List[AnalysisModel]:
        """Get recent analyses"""
        async with self.async_session() as session:
            result = await session.execute(
                select(Analysis)
                .order_by(Analysis.created_at.desc())
                .limit(limit)
            )
            analyses = result.scalars().all()
            
            return [self._analysis_to_model(analysis) for analysis in analyses]
    
    # System metrics
    async def record_metric(self, metric_name: str, metric_value: float):
        """Record a system metric"""
        async with self.async_session() as session:
            metric = SystemMetrics(
                metric_name=metric_name,
                metric_value=metric_value
            )
            session.add(metric)
            await session.commit()
    
    async def get_metrics(self, metric_name: str, hours: int = 24) -> List[Dict[str, Any]]:
        """Get metrics for the last N hours"""
        async with self.async_session() as session:
            cutoff_time = datetime.now() - timedelta(hours=hours)
            
            result = await session.execute(
                select(SystemMetrics)
                .where(
                    SystemMetrics.metric_name == metric_name,
                    SystemMetrics.recorded_at >= cutoff_time
                )
                .order_by(SystemMetrics.recorded_at)
            )
            metrics = result.scalars().all()
            
            return [
                {
                    "value": metric.metric_value,
                    "timestamp": metric.recorded_at
                }
                for metric in metrics
            ]
    
    # Helper methods
    def _analysis_to_model(self, analysis: Analysis) -> AnalysisModel:
        """Convert Analysis SQLAlchemy model to Pydantic model"""
        return AnalysisModel(
            id=analysis.id,
            upload_id=analysis.upload_id,
            sample_id=analysis.sample_id,
            filename=analysis.filename,
            file_path=analysis.file_path,
            analysis_type=analysis.analysis_type,
            status=analysis.status,
            user_id=analysis.user_id,
            description=analysis.description,
            results=analysis.results,
            error_message=analysis.error_message,
            report_html_path=analysis.report_html_path,
            report_json_path=analysis.report_json_path,
            created_at=analysis.created_at,
            started_at=analysis.started_at,
            completed_at=analysis.completed_at
        )
    
    # Raw SQL methods for complex queries
    async def execute_raw_query(self, query: str, params: Dict[str, Any] = None) -> List[Dict[str, Any]]:
        """Execute raw SQL query"""
        async with self.pool.acquire() as conn:
            if params:
                result = await conn.fetch(query, **params)
            else:
                result = await conn.fetch(query)
            
            return [dict(row) for row in result]
    
    async def get_analysis_stats_by_user(self, user_id: int) -> Dict[str, Any]:
        """Get analysis statistics for a user"""
        query = """
        SELECT 
            COUNT(*) as total_analyses,
            COUNT(CASE WHEN status = 'completed' THEN 1 END) as completed_analyses,
            COUNT(CASE WHEN status = 'failed' THEN 1 END) as failed_analyses,
            COUNT(CASE WHEN status = 'running' THEN 1 END) as running_analyses,
            AVG(EXTRACT(EPOCH FROM (completed_at - started_at))) as avg_duration_seconds
        FROM analyses 
        WHERE user_id = $1
        """
        
        async with self.pool.acquire() as conn:
            result = await conn.fetchrow(query, user_id)
            return dict(result) if result else {}