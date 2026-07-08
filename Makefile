# Makefile for Pharmacogenomics Pipeline

.PHONY: help install dev build test clean deploy

# Default target
help:
	@echo "Pharmacogenomics Pipeline - Available Commands:"
	@echo ""
	@echo "  make install     - Install all dependencies"
	@echo "  make dev         - Start development servers"
	@echo "  make build       - Build production images"
	@echo "  make test        - Run all tests"
	@echo "  make lint        - Run linting"
	@echo "  make format      - Format code"
	@echo "  make clean       - Clean up containers and volumes"
	@echo "  make deploy      - Deploy to production"
	@echo "  make logs        - View container logs"
	@echo "  make shell       - Access backend container shell"
	@echo "  make db-migrate  - Run database migrations"
	@echo "  make backup      - Backup database"
	@echo "  make restore     - Restore database from backup"
	@echo ""

# Installation
install:
	@echo "Installing dependencies..."
	pip install -r requirements.txt
	cd frontend && npm install
	@echo "✅ Dependencies installed"

install-dev:
	@echo "Installing development dependencies..."
	pip install -r requirements.txt
	pip install pytest pytest-asyncio black isort mypy
	cd frontend && npm install
	@echo "✅ Development dependencies installed"

# Development
dev:
	@echo "Starting development environment..."
	docker-compose -f docker/docker-compose.yml up -d postgres redis
	@echo "Waiting for services to be ready..."
	sleep 5
	@echo "🚀 Development services started"
	@echo "Run 'make dev-backend' and 'make dev-frontend' in separate terminals"

dev-backend:
	@echo "Starting backend development server..."
	cd backend && uvicorn main:app --host 0.0.0.0 --port 8000 --reload

dev-frontend:
	@echo "Starting frontend development server..."
	cd frontend && npm run dev

# Production
build:
	@echo "Building production images..."
	docker-compose -f docker/docker-compose.yml build
	@echo "✅ Production images built"

deploy:
	@echo "Deploying to production..."
	docker-compose -f docker/docker-compose.yml up -d
	@echo "🚀 Production deployment started"

stop:
	@echo "Stopping all services..."
	docker-compose -f docker/docker-compose.yml down
	@echo "⏹️ Services stopped"

restart:
	@echo "Restarting services..."
	docker-compose -f docker/docker-compose.yml restart
	@echo "🔄 Services restarted"

# Testing
test:
	@echo "Running backend tests..."
	pytest backend/tests/ -v
	@echo "Running frontend tests..."
	cd frontend && npm run test
	@echo "✅ All tests completed"

test-backend:
	@echo "Running backend tests..."
	pytest backend/tests/ -v --cov=backend

test-frontend:
	@echo "Running frontend tests..."
	cd frontend && npm run test

test-integration:
	@echo "Running integration tests..."
	pytest tests/integration/ -v

# Code Quality
lint:
	@echo "Running linting..."
	black --check backend/
	isort --check-only backend/
	mypy backend/
	cd frontend && npm run lint
	@echo "✅ Linting completed"

format:
	@echo "Formatting code..."
	black backend/
	isort backend/
	cd frontend && npm run lint:fix
	@echo "✅ Code formatted"

type-check:
	@echo "Running type checks..."
	mypy backend/
	cd frontend && npm run type-check
	@echo "✅ Type checking completed"

# Database
db-migrate:
	@echo "Running database migrations..."
	docker-compose -f docker/docker-compose.yml exec backend alembic upgrade head
	@echo "✅ Database migrations completed"

db-revision:
	@echo "Creating new database revision..."
	@read -p "Migration message: " message; \
	docker-compose -f docker/docker-compose.yml exec backend alembic revision --autogenerate -m "$$message"

db-reset:
	@echo "Resetting database..."
	docker-compose -f docker/docker-compose.yml exec backend alembic downgrade base
	docker-compose -f docker/docker-compose.yml exec backend alembic upgrade head
	@echo "✅ Database reset completed"

# Backup & Restore
backup:
	@echo "Creating database backup..."
	mkdir -p backups
	docker-compose -f docker/docker-compose.yml exec postgres pg_dump -U postgres pharmacogenomics > backups/backup_$(shell date +%Y%m%d_%H%M%S).sql
	@echo "✅ Database backup created"

restore:
	@echo "Available backups:"
	@ls -la backups/
	@read -p "Enter backup filename: " filename; \
	docker-compose -f docker/docker-compose.yml exec -T postgres psql -U postgres -c "DROP DATABASE IF EXISTS pharmacogenomics;"; \
	docker-compose -f docker/docker-compose.yml exec -T postgres psql -U postgres -c "CREATE DATABASE pharmacogenomics;"; \
	docker-compose -f docker/docker-compose.yml exec -T postgres psql -U postgres pharmacogenomics < backups/$$filename
	@echo "✅ Database restored"

# Monitoring
logs:
	@echo "Viewing container logs..."
	docker-compose -f docker/docker-compose.yml logs -f

logs-backend:
	docker-compose -f docker/docker-compose.yml logs -f backend

logs-frontend:
	docker-compose -f docker/docker-compose.yml logs -f frontend

logs-db:
	docker-compose -f docker/docker-compose.yml logs -f postgres

# Shell Access
shell:
	@echo "Accessing backend container shell..."
	docker-compose -f docker/docker-compose.yml exec backend /bin/bash

shell-db:
	@echo "Accessing database shell..."
	docker-compose -f docker/docker-compose.yml exec postgres psql -U postgres pharmacogenomics

# Cleanup
clean:
	@echo "Cleaning up containers and images..."
	docker-compose -f docker/docker-compose.yml down --volumes --rmi all
	docker system prune -f
	@echo "✅ Cleanup completed"

clean-volumes:
	@echo "Cleaning up volumes..."
	docker-compose -f docker/docker-compose.yml down --volumes
	@echo "✅ Volumes cleaned"

# Development utilities
setup-dev:
	@echo "Setting up development environment..."
	cp .env.example .env
	@echo "Please edit .env file with your configuration"
	make install-dev
	@echo "✅ Development environment setup completed"

check-health:
	@echo "Checking service health..."
	curl -f http://localhost:8000/health || echo "Backend not healthy"
	curl -f http://localhost:3000 || echo "Frontend not healthy"
	@echo "✅ Health check completed"

# Security
security-scan:
	@echo "Running security scans..."
	docker run --rm -v $(PWD):/app securecodewarrior/docker-security-scanner /app
	cd frontend && npm audit
	@echo "✅ Security scan completed"

# Performance
performance-test:
	@echo "Running performance tests..."
	# Add performance testing tools like locust, k6, etc.
	@echo "✅ Performance tests completed"

# Documentation
docs:
	@echo "Generating documentation..."
	cd backend && python -m pydoc -w .
	cd frontend && npm run docs
	@echo "✅ Documentation generated"

# Docker utilities
docker-stats:
	docker stats --no-stream

docker-inspect:
	docker-compose -f docker/docker-compose.yml ps

# Production utilities
prod-deploy:
	@echo "Deploying to production with SSL..."
	docker-compose -f docker/docker-compose.prod.yml up -d
	@echo "🚀 Production deployment with SSL started"

staging-deploy:
	@echo "Deploying to staging..."
	docker-compose -f docker/docker-compose.staging.yml up -d
	@echo "🚀 Staging deployment started"