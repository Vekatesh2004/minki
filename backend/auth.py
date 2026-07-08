"""
Authentication and authorization for Pharmacogenomics Pipeline
"""

import os
import secrets
from datetime import datetime, timedelta
from typing import Optional

import bcrypt
from jose import JWTError, jwt
from passlib.context import CryptContext
import structlog

from .models import UserCreate, UserModel, UserLogin
from .database import DatabaseManager

logger = structlog.get_logger()

class AuthManager:
    """Authentication manager"""
    
    def __init__(self, db_manager: DatabaseManager):
        self.db_manager = db_manager
        self.pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
        
        # JWT configuration
        self.secret_key = os.getenv("SECRET_KEY", self._generate_secret_key())
        self.algorithm = "HS256"
        self.access_token_expire_minutes = int(os.getenv("ACCESS_TOKEN_EXPIRE_MINUTES", "1440"))  # 24 hours
    
    def _generate_secret_key(self) -> str:
        """Generate a secure secret key"""
        return secrets.token_urlsafe(32)
    
    def verify_password(self, plain_password: str, hashed_password: str) -> bool:
        """Verify password against hash"""
        try:
            return self.pwd_context.verify(plain_password, hashed_password)
        except Exception as e:
            logger.error("Password verification failed", error=str(e))
            return False
    
    def hash_password(self, password: str) -> str:
        """Hash password"""
        return self.pwd_context.hash(password)
    
    def create_access_token(self, data: dict, expires_delta: Optional[timedelta] = None):
        """Create JWT access token"""
        to_encode = data.copy()
        
        if expires_delta:
            expire = datetime.utcnow() + expires_delta
        else:
            expire = datetime.utcnow() + timedelta(minutes=self.access_token_expire_minutes)
        
        to_encode.update({"exp": expire})
        encoded_jwt = jwt.encode(to_encode, self.secret_key, algorithm=self.algorithm)
        
        return encoded_jwt
    
    async def create_user(self, user_data: UserCreate) -> UserModel:
        """Create a new user"""
        try:
            # Hash password
            hashed_password = self.hash_password(user_data.password)
            
            # Create user data with hashed password
            user_create = UserCreate(
                email=user_data.email,
                full_name=user_data.full_name,
                password=hashed_password,
                is_active=user_data.is_active
            )
            
            # Create user in database
            user = await self.db_manager.create_user(user_create)
            
            logger.info("User created successfully", user_id=user.id, email=user.email)
            return user
            
        except ValueError as e:
            logger.warning("User creation failed", error=str(e), email=user_data.email)
            raise
        except Exception as e:
            logger.error("User creation failed", error=str(e), email=user_data.email)
            raise ValueError("Failed to create user")
    
    async def authenticate_user(self, email: str, password: str) -> str:
        """Authenticate user and return JWT token"""
        try:
            # Get user from database
            user = await self.db_manager.get_user_by_email(email)
            
            if not user:
                logger.warning("Authentication failed - user not found", email=email)
                raise ValueError("Invalid email or password")
            
            if not user.is_active:
                logger.warning("Authentication failed - user inactive", email=email)
                raise ValueError("Account is deactivated")
            
            # Get user with password hash for verification
            # Note: In a real implementation, you'd need to modify the database 
            # model to include the hashed password in a secure way
            async with self.db_manager.async_session() as session:
                from .database import User
                from sqlalchemy import select
                
                result = await session.execute(
                    select(User).where(User.email == email)
                )
                db_user = result.scalar_one_or_none()
                
                if not db_user:
                    raise ValueError("Invalid email or password")
                
                # Verify password
                if not self.verify_password(password, db_user.hashed_password):
                    logger.warning("Authentication failed - invalid password", email=email)
                    raise ValueError("Invalid email or password")
            
            # Update last login
            await self.db_manager.update_user_login(user.id)
            
            # Create access token
            access_token_expires = timedelta(minutes=self.access_token_expire_minutes)
            access_token = self.create_access_token(
                data={"sub": str(user.id), "email": user.email},
                expires_delta=access_token_expires
            )
            
            logger.info("User authenticated successfully", user_id=user.id, email=email)
            return access_token
            
        except ValueError:
            raise
        except Exception as e:
            logger.error("Authentication failed", error=str(e), email=email)
            raise ValueError("Authentication failed")
    
    async def verify_token(self, token: str) -> Optional[UserModel]:
        """Verify JWT token and return user"""
        try:
            # Decode token
            payload = jwt.decode(token, self.secret_key, algorithms=[self.algorithm])
            user_id: str = payload.get("sub")
            
            if user_id is None:
                return None
            
            # Get user from database
            user = await self.db_manager.get_user_by_id(int(user_id))
            
            if not user or not user.is_active:
                return None
            
            return user
            
        except JWTError as e:
            logger.warning("Token verification failed", error=str(e))
            return None
        except Exception as e:
            logger.error("Token verification error", error=str(e))
            return None
    
    async def refresh_token(self, token: str) -> Optional[str]:
        """Refresh JWT token"""
        try:
            # Verify current token
            user = await self.verify_token(token)
            
            if not user:
                return None
            
            # Create new token
            access_token_expires = timedelta(minutes=self.access_token_expire_minutes)
            new_token = self.create_access_token(
                data={"sub": str(user.id), "email": user.email},
                expires_delta=access_token_expires
            )
            
            logger.info("Token refreshed successfully", user_id=user.id)
            return new_token
            
        except Exception as e:
            logger.error("Token refresh failed", error=str(e))
            return None
    
    async def invalidate_token(self, token: str) -> bool:
        """Invalidate JWT token (add to blacklist)"""
        # In a production system, you'd maintain a token blacklist
        # For now, we'll just log the action
        try:
            payload = jwt.decode(token, self.secret_key, algorithms=[self.algorithm])
            user_id = payload.get("sub")
            
            logger.info("Token invalidated", user_id=user_id)
            return True
            
        except Exception as e:
            logger.error("Token invalidation failed", error=str(e))
            return False
    
    async def change_password(self, user_id: int, current_password: str, new_password: str) -> bool:
        """Change user password"""
        try:
            # Get user
            user = await self.db_manager.get_user_by_id(user_id)
            if not user:
                return False
            
            # Get user with password hash
            async with self.db_manager.async_session() as session:
                from .database import User
                from sqlalchemy import select, update
                
                result = await session.execute(
                    select(User).where(User.id == user_id)
                )
                db_user = result.scalar_one_or_none()
                
                if not db_user:
                    return False
                
                # Verify current password
                if not self.verify_password(current_password, db_user.hashed_password):
                    logger.warning("Password change failed - invalid current password", user_id=user_id)
                    return False
                
                # Hash new password and update
                new_hashed_password = self.hash_password(new_password)
                
                await session.execute(
                    update(User)
                    .where(User.id == user_id)
                    .values(hashed_password=new_hashed_password)
                )
                await session.commit()
            
            logger.info("Password changed successfully", user_id=user_id)
            return True
            
        except Exception as e:
            logger.error("Password change failed", error=str(e), user_id=user_id)
            return False
    
    async def reset_password(self, email: str) -> Optional[str]:
        """Generate password reset token"""
        try:
            # Get user
            user = await self.db_manager.get_user_by_email(email)
            if not user:
                # Don't reveal if email exists
                logger.warning("Password reset requested for non-existent email", email=email)
                return None
            
            # Create reset token (expires in 1 hour)
            reset_token = self.create_access_token(
                data={"sub": str(user.id), "type": "password_reset"},
                expires_delta=timedelta(hours=1)
            )
            
            logger.info("Password reset token generated", user_id=user.id, email=email)
            return reset_token
            
        except Exception as e:
            logger.error("Password reset token generation failed", error=str(e), email=email)
            return None
    
    async def verify_reset_token(self, token: str, new_password: str) -> bool:
        """Verify reset token and update password"""
        try:
            # Decode token
            payload = jwt.decode(token, self.secret_key, algorithms=[self.algorithm])
            user_id = payload.get("sub")
            token_type = payload.get("type")
            
            if user_id is None or token_type != "password_reset":
                return False
            
            # Update password
            new_hashed_password = self.hash_password(new_password)
            
            async with self.db_manager.async_session() as session:
                from .database import User
                from sqlalchemy import update
                
                await session.execute(
                    update(User)
                    .where(User.id == int(user_id))
                    .values(hashed_password=new_hashed_password)
                )
                await session.commit()
            
            logger.info("Password reset completed", user_id=user_id)
            return True
            
        except JWTError as e:
            logger.warning("Invalid password reset token", error=str(e))
            return False
        except Exception as e:
            logger.error("Password reset failed", error=str(e))
            return False
    
    def generate_api_key(self) -> str:
        """Generate API key for external integrations"""
        return f"pgx_{secrets.token_urlsafe(32)}"
    
    async def verify_api_key(self, api_key: str) -> Optional[UserModel]:
        """Verify API key and return associated user"""
        try:
            # In production, store API keys hashed in database
            # For now, this is a placeholder
            
            # Hash the provided key
            key_hash = self.hash_password(api_key)
            
            # Query database for matching API key
            async with self.db_manager.async_session() as session:
                from .database import APIKey, User
                from sqlalchemy import select, join
                
                result = await session.execute(
                    select(User)
                    .join(APIKey)
                    .where(
                        APIKey.key_hash == key_hash,
                        APIKey.is_active == True,
                        User.is_active == True
                    )
                )
                user_row = result.scalar_one_or_none()
                
                if user_row:
                    # Update last used timestamp
                    from sqlalchemy import update
                    await session.execute(
                        update(APIKey)
                        .where(APIKey.key_hash == key_hash)
                        .values(last_used=datetime.now())
                    )
                    await session.commit()
                    
                    return UserModel(
                        id=user_row.id,
                        email=user_row.email,
                        full_name=user_row.full_name,
                        is_active=user_row.is_active,
                        created_at=user_row.created_at,
                        last_login=user_row.last_login
                    )
            
            return None
            
        except Exception as e:
            logger.error("API key verification failed", error=str(e))
            return None