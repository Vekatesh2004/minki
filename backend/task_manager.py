"""
Task management for background processing
"""

import asyncio
import json
import os
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List
from enum import Enum
import uuid

import redis.asyncio as redis
from celery import Celery
import structlog

logger = structlog.get_logger()

class TaskStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"

class TaskManager:
    """Manages background tasks and job queues"""
    
    def __init__(self):
        self.redis_client = None
        self.celery_app = None
        self.task_registry = {}
        
    async def initialize(self):
        """Initialize task manager"""
        try:
            # Initialize Redis connection
            redis_url = os.getenv("REDIS_URL", "redis://localhost:6379/0")
            self.redis_client = redis.from_url(redis_url)
            
            # Test Redis connection
            await self.redis_client.ping()
            
            # Initialize Celery if available
            try:
                broker_url = os.getenv("CELERY_BROKER_URL", redis_url)
                result_backend = os.getenv("CELERY_RESULT_BACKEND", redis_url)
                
                self.celery_app = Celery(
                    'pharmacogenomics_tasks',
                    broker=broker_url,
                    backend=result_backend
                )
                
                # Configure Celery
                self.celery_app.conf.update(
                    task_serializer='json',
                    accept_content=['json'],
                    result_serializer='json',
                    timezone='UTC',
                    enable_utc=True,
                    task_track_started=True,
                    task_time_limit=3600,  # 1 hour
                    worker_prefetch_multiplier=1,
                    result_expires=86400,  # 24 hours
                )
                
            except Exception as e:
                logger.warning("Celery not available, using simple task queue", error=str(e))
                self.celery_app = None
            
            logger.info("Task manager initialized successfully")
            
        except Exception as e:
            logger.error("Failed to initialize task manager", error=str(e))
            raise
    
    async def close(self):
        """Close task manager connections"""
        if self.redis_client:
            await self.redis_client.close()
    
    async def queue_task(self, task_name: str, task_data: Dict[str, Any], 
                        priority: int = 5, delay: Optional[int] = None) -> str:
        """Queue a task for background processing"""
        
        task_id = str(uuid.uuid4())
        
        task_info = {
            "task_id": task_id,
            "task_name": task_name,
            "task_data": task_data,
            "status": TaskStatus.PENDING,
            "priority": priority,
            "created_at": datetime.now().isoformat(),
            "scheduled_at": (datetime.now() + timedelta(seconds=delay)).isoformat() if delay else None,
            "started_at": None,
            "completed_at": None,
            "progress": 0.0,
            "result": None,
            "error": None
        }
        
        try:
            if self.celery_app:
                # Use Celery for task management
                celery_task = self.celery_app.send_task(
                    task_name,
                    args=[task_data],
                    task_id=task_id,
                    countdown=delay,
                    priority=priority
                )
                
                # Store task info in Redis
                await self.redis_client.set(
                    f"task:{task_id}",
                    json.dumps(task_info),
                    ex=86400  # Expire after 24 hours
                )
                
            else:
                # Use simple Redis-based queue
                await self._queue_simple_task(task_info)
            
            logger.info("Task queued", task_id=task_id, task_name=task_name)
            return task_id
            
        except Exception as e:
            logger.error("Failed to queue task", error=str(e), task_name=task_name)
            raise
    
    async def _queue_simple_task(self, task_info: Dict[str, Any]):
        """Queue task in simple Redis queue"""
        
        # Add to priority queue
        priority = task_info.get("priority", 5)
        queue_key = f"queue:priority:{priority}"
        
        await self.redis_client.lpush(queue_key, json.dumps(task_info))
        
        # Add to task tracking
        await self.redis_client.set(
            f"task:{task_info['task_id']}",
            json.dumps(task_info),
            ex=86400
        )
    
    async def get_task_status(self, task_id: str) -> Optional[Dict[str, Any]]:
        """Get task status and information"""
        
        try:
            task_data = await self.redis_client.get(f"task:{task_id}")
            
            if not task_data:
                return None
            
            task_info = json.loads(task_data)
            
            # If using Celery, get additional status info
            if self.celery_app:
                celery_result = self.celery_app.AsyncResult(task_id)
                
                if celery_result.state == "PENDING":
                    task_info["status"] = TaskStatus.PENDING
                elif celery_result.state == "STARTED":
                    task_info["status"] = TaskStatus.RUNNING
                    if not task_info.get("started_at"):
                        task_info["started_at"] = datetime.now().isoformat()
                elif celery_result.state == "SUCCESS":
                    task_info["status"] = TaskStatus.COMPLETED
                    task_info["result"] = celery_result.result
                    if not task_info.get("completed_at"):
                        task_info["completed_at"] = datetime.now().isoformat()
                elif celery_result.state == "FAILURE":
                    task_info["status"] = TaskStatus.FAILED
                    task_info["error"] = str(celery_result.info)
                    if not task_info.get("completed_at"):
                        task_info["completed_at"] = datetime.now().isoformat()
                
                # Update stored task info
                await self.redis_client.set(
                    f"task:{task_id}",
                    json.dumps(task_info),
                    ex=86400
                )
            
            return task_info
            
        except Exception as e:
            logger.error("Failed to get task status", error=str(e), task_id=task_id)
            return None
    
    async def update_task_progress(self, task_id: str, progress: float, 
                                  status: Optional[TaskStatus] = None, 
                                  message: Optional[str] = None):
        """Update task progress"""
        
        try:
            task_data = await self.redis_client.get(f"task:{task_id}")
            
            if not task_data:
                logger.warning("Task not found for progress update", task_id=task_id)
                return
            
            task_info = json.loads(task_data)
            task_info["progress"] = progress
            
            if status:
                task_info["status"] = status
                
                if status == TaskStatus.RUNNING and not task_info.get("started_at"):
                    task_info["started_at"] = datetime.now().isoformat()
                elif status in [TaskStatus.COMPLETED, TaskStatus.FAILED] and not task_info.get("completed_at"):
                    task_info["completed_at"] = datetime.now().isoformat()
            
            if message:
                task_info["message"] = message
            
            await self.redis_client.set(
                f"task:{task_id}",
                json.dumps(task_info),
                ex=86400
            )
            
            # Publish progress update for real-time notifications
            await self.redis_client.publish(
                f"task_progress:{task_id}",
                json.dumps({
                    "task_id": task_id,
                    "progress": progress,
                    "status": status,
                    "message": message
                })
            )
            
        except Exception as e:
            logger.error("Failed to update task progress", error=str(e), task_id=task_id)
    
    async def complete_task(self, task_id: str, result: Any, 
                           status: TaskStatus = TaskStatus.COMPLETED):
        """Mark task as completed with result"""
        
        try:
            task_data = await self.redis_client.get(f"task:{task_id}")
            
            if not task_data:
                logger.warning("Task not found for completion", task_id=task_id)
                return
            
            task_info = json.loads(task_data)
            task_info["status"] = status
            task_info["result"] = result
            task_info["progress"] = 100.0
            task_info["completed_at"] = datetime.now().isoformat()
            
            await self.redis_client.set(
                f"task:{task_id}",
                json.dumps(task_info),
                ex=86400
            )
            
            # Publish completion notification
            await self.redis_client.publish(
                f"task_complete:{task_id}",
                json.dumps({
                    "task_id": task_id,
                    "status": status,
                    "result": result
                })
            )
            
            logger.info("Task completed", task_id=task_id, status=status)
            
        except Exception as e:
            logger.error("Failed to complete task", error=str(e), task_id=task_id)
    
    async def fail_task(self, task_id: str, error: str):
        """Mark task as failed with error"""
        
        try:
            task_data = await self.redis_client.get(f"task:{task_id}")
            
            if not task_data:
                logger.warning("Task not found for failure", task_id=task_id)
                return
            
            task_info = json.loads(task_data)
            task_info["status"] = TaskStatus.FAILED
            task_info["error"] = error
            task_info["completed_at"] = datetime.now().isoformat()
            
            await self.redis_client.set(
                f"task:{task_id}",
                json.dumps(task_info),
                ex=86400
            )
            
            # Publish failure notification
            await self.redis_client.publish(
                f"task_failed:{task_id}",
                json.dumps({
                    "task_id": task_id,
                    "error": error
                })
            )
            
            logger.error("Task failed", task_id=task_id, error=error)
            
        except Exception as e:
            logger.error("Failed to mark task as failed", error=str(e), task_id=task_id)
    
    async def cancel_task(self, task_id: str) -> bool:
        """Cancel a pending or running task"""
        
        try:
            if self.celery_app:
                # Revoke Celery task
                self.celery_app.control.revoke(task_id, terminate=True)
            
            # Update task status
            task_data = await self.redis_client.get(f"task:{task_id}")
            
            if task_data:
                task_info = json.loads(task_data)
                task_info["status"] = TaskStatus.CANCELLED
                task_info["completed_at"] = datetime.now().isoformat()
                
                await self.redis_client.set(
                    f"task:{task_id}",
                    json.dumps(task_info),
                    ex=86400
                )
            
            logger.info("Task cancelled", task_id=task_id)
            return True
            
        except Exception as e:
            logger.error("Failed to cancel task", error=str(e), task_id=task_id)
            return False
    
    async def get_queue_stats(self) -> Dict[str, Any]:
        """Get queue statistics"""
        
        try:
            stats = {
                "total_tasks": 0,
                "pending_tasks": 0,
                "running_tasks": 0,
                "completed_tasks": 0,
                "failed_tasks": 0,
                "queue_lengths": {}
            }
            
            # Get all task keys
            task_keys = await self.redis_client.keys("task:*")
            stats["total_tasks"] = len(task_keys)
            
            # Count tasks by status
            for key in task_keys:
                task_data = await self.redis_client.get(key)
                if task_data:
                    task_info = json.loads(task_data)
                    status = task_info.get("status", TaskStatus.PENDING)
                    
                    if status == TaskStatus.PENDING:
                        stats["pending_tasks"] += 1
                    elif status == TaskStatus.RUNNING:
                        stats["running_tasks"] += 1
                    elif status == TaskStatus.COMPLETED:
                        stats["completed_tasks"] += 1
                    elif status == TaskStatus.FAILED:
                        stats["failed_tasks"] += 1
            
            # Get queue lengths
            for priority in range(1, 11):  # Priority 1-10
                queue_key = f"queue:priority:{priority}"
                queue_length = await self.redis_client.llen(queue_key)
                if queue_length > 0:
                    stats["queue_lengths"][f"priority_{priority}"] = queue_length
            
            return stats
            
        except Exception as e:
            logger.error("Failed to get queue stats", error=str(e))
            return {}
    
    async def cleanup_old_tasks(self, max_age_hours: int = 48):
        """Clean up old completed/failed tasks"""
        
        try:
            cutoff_time = datetime.now() - timedelta(hours=max_age_hours)
            
            # Get all task keys
            task_keys = await self.redis_client.keys("task:*")
            deleted_count = 0
            
            for key in task_keys:
                task_data = await self.redis_client.get(key)
                if task_data:
                    task_info = json.loads(task_data)
                    
                    # Check if task is old and completed/failed
                    created_at = datetime.fromisoformat(task_info.get("created_at", ""))
                    status = task_info.get("status")
                    
                    if (created_at < cutoff_time and 
                        status in [TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED]):
                        
                        await self.redis_client.delete(key)
                        deleted_count += 1
            
            logger.info("Cleaned up old tasks", deleted_count=deleted_count)
            return deleted_count
            
        except Exception as e:
            logger.error("Failed to cleanup old tasks", error=str(e))
            return 0
    
    async def get_user_tasks(self, user_id: int, limit: int = 50) -> List[Dict[str, Any]]:
        """Get tasks for a specific user"""
        
        try:
            # This is a simplified implementation
            # In production, you'd want to index tasks by user_id
            
            task_keys = await self.redis_client.keys("task:*")
            user_tasks = []
            
            for key in task_keys:
                if len(user_tasks) >= limit:
                    break
                    
                task_data = await self.redis_client.get(key)
                if task_data:
                    task_info = json.loads(task_data)
                    
                    # Check if task belongs to user (stored in task_data)
                    task_user_id = task_info.get("task_data", {}).get("user_id")
                    if task_user_id == user_id:
                        user_tasks.append(task_info)
            
            # Sort by created_at (newest first)
            user_tasks.sort(
                key=lambda x: x.get("created_at", ""),
                reverse=True
            )
            
            return user_tasks
            
        except Exception as e:
            logger.error("Failed to get user tasks", error=str(e), user_id=user_id)
            return []