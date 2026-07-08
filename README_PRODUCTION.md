# Pharmacogenomics Pipeline - Production Deployment Guide

## 🚀 Production-Ready Architecture

This is a comprehensive, production-ready web application for pharmacogenomics analysis with the following architecture:

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   React SPA     │    │  FastAPI + API  │    │  PostgreSQL     │
│   (Frontend)    │◄──►│   (Backend)     │◄──►│   (Database)    │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         │              ┌─────────────────┐    ┌─────────────────┐
         │              │  Celery Workers │    │     Redis       │
         │              │ (Bg Processing) │◄──►│   (Cache/Queue) │
         │              └─────────────────┘    └─────────────────┘
         │
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   Nginx Proxy   │    │   Prometheus    │    │    Grafana      │
│ (Load Balancer) │    │   (Metrics)     │    │  (Monitoring)   │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

## 🛠️ Technology Stack

### Backend (FastAPI)
- **FastAPI**: Modern, fast web framework for building APIs
- **SQLAlchemy**: Async ORM with PostgreSQL
- **Celery**: Distributed task queue for background processing
- **Redis**: Caching and message broker
- **Pydantic**: Data validation and serialization
- **JWT**: Authentication and authorization
- **Prometheus**: Metrics collection
- **Structured Logging**: JSON-formatted logs

### Frontend (React)
- **React 18**: Modern React with concurrent features
- **TypeScript**: Type-safe JavaScript
- **Vite**: Fast build tool and dev server
- **TailwindCSS**: Utility-first CSS framework
- **React Query**: Server state management
- **React Router**: Client-side routing
- **React Hook Form**: Form handling
- **Framer Motion**: Animations
- **Recharts**: Data visualization

### Infrastructure
- **Docker**: Containerization
- **Nginx**: Reverse proxy and static file serving
- **PostgreSQL**: Primary database
- **Redis**: Caching and task queue
- **Prometheus + Grafana**: Monitoring and alerting

## 🚀 Quick Start

### 1. Prerequisites

```bash
# Install required software
- Docker & Docker Compose
- Make (optional, for convenience commands)
- Git
```

### 2. Clone and Setup

```bash
# Clone the repository
git clone <your-repo-url>
cd minki

# Copy environment configuration
cp .env.example .env

# Edit configuration (IMPORTANT!)
nano .env
```

### 3. Start Development Environment

```bash
# Using Make (recommended)
make dev

# Or using Docker Compose directly
docker-compose -f docker/docker-compose.yml up -d
```

### 4. Access the Application

- **Frontend**: http://localhost:3000
- **Backend API**: http://localhost:8000
- **API Documentation**: http://localhost:8000/docs
- **Grafana**: http://localhost:3001
- **Prometheus**: http://localhost:9090

## 🔧 Development Workflow

### Backend Development

```bash
# Start backend in development mode
make dev-backend

# Run tests
make test-backend

# Format code
make format

# Type checking
make type-check
```

### Frontend Development

```bash
# Start frontend in development mode
make dev-frontend

# Run tests
make test-frontend

# Build for production
cd frontend && npm run build
```

### Database Operations

```bash
# Run migrations
make db-migrate

# Create new migration
make db-revision

# Reset database
make db-reset

# Backup database
make backup

# Restore database
make restore
```

## 🚀 Production Deployment

### 1. Server Requirements

**Minimum Requirements:**
- CPU: 4 cores
- RAM: 8GB
- Storage: 100GB SSD
- Network: 100 Mbps

**Recommended for Production:**
- CPU: 8+ cores
- RAM: 16GB+
- Storage: 500GB+ SSD
- Network: 1 Gbps
- Load balancer (for high availability)

### 2. Production Setup

```bash
# 1. Clone repository on production server
git clone <your-repo-url>
cd minki

# 2. Configure environment
cp .env.example .env
# Edit .env with production values:
# - Strong passwords
# - Production database URLs
# - SSL certificate paths
# - External API keys

# 3. Build and deploy
make build
make deploy

# 4. Setup SSL certificates (Let's Encrypt recommended)
# Follow SSL setup guide below
```

### 3. SSL Configuration

```bash
# Install certbot
sudo apt install certbot

# Generate SSL certificates
sudo certbot certonly --standalone -d yourdomain.com

# Update nginx configuration with SSL paths
# Uncomment SSL blocks in docker/nginx.conf
```

### 4. Environment Variables (Production)

```env
# Security (CHANGE THESE!)
SECRET_KEY=your-ultra-secure-secret-key-at-least-32-chars
POSTGRES_PASSWORD=ultra-secure-db-password
REDIS_PASSWORD=ultra-secure-redis-password

# Database
DATABASE_URL=postgresql+asyncpg://postgres:password@postgres:5432/pharmacogenomics

# API Keys
DRUGBANK_API_KEY=your-drugbank-api-key
ENSEMBL_VEP_API_KEY=optional-vep-api-key

# Domain configuration
CORS_ORIGINS=["https://yourdomain.com"]
VITE_API_URL=https://yourdomain.com

# Monitoring
GRAFANA_PASSWORD=secure-grafana-password

# Email notifications
SMTP_HOST=smtp.gmail.com
SMTP_USER=your-email@domain.com
SMTP_PASSWORD=your-app-password
```

## 📊 Monitoring and Logging

### Prometheus Metrics
- HTTP request metrics
- Analysis duration
- Queue sizes
- Database connections
- Custom business metrics

### Grafana Dashboards
- System overview
- Application performance
- Analysis pipeline metrics
- Error rates and alerts

### Log Aggregation
```bash
# View logs
make logs

# View specific service logs
make logs-backend
make logs-frontend

# Monitor live logs
docker-compose -f docker/docker-compose.yml logs -f --tail=100
```

## 🔒 Security Features

### Authentication & Authorization
- JWT-based authentication
- Role-based access control
- API key authentication for external integrations
- Password hashing (bcrypt)

### Data Security
- Input validation (Pydantic)
- SQL injection prevention (SQLAlchemy)
- File upload restrictions
- Rate limiting
- CORS configuration

### Infrastructure Security
- Non-root container users
- Network isolation
- SSL/TLS encryption
- Security headers
- Regular security updates

## 📈 Scalability

### Horizontal Scaling
```bash
# Scale backend workers
docker-compose -f docker/docker-compose.yml up -d --scale backend=3

# Scale Celery workers
docker-compose -f docker/docker-compose.yml up -d --scale celery-worker=5
```

### Load Balancing
```nginx
# Add to nginx.conf upstream block
upstream backend {
    server backend1:8000;
    server backend2:8000;
    server backend3:8000;
}
```

### Database Optimization
- Connection pooling
- Read replicas
- Query optimization
- Proper indexing

## 🔧 Maintenance

### Regular Tasks
```bash
# Update dependencies
make install

# Security updates
make security-scan

# Database maintenance
make db-migrate

# Cleanup old data
docker-compose exec backend python scripts/cleanup.py

# Backup database
make backup
```

### Health Checks
```bash
# Check system health
make check-health

# View system status
curl http://localhost:8000/health

# Monitor resource usage
make docker-stats
```

## 🚨 Troubleshooting

### Common Issues

**1. Database Connection Issues**
```bash
# Check database status
docker-compose ps postgres

# Check database logs
make logs-db

# Reset database connection
docker-compose restart postgres
```

**2. High Memory Usage**
```bash
# Check memory usage
docker stats

# Restart services
make restart

# Scale down if needed
docker-compose scale backend=1
```

**3. SSL Certificate Issues**
```bash
# Renew certificates
sudo certbot renew

# Restart nginx
docker-compose restart nginx
```

## 📚 API Documentation

### Interactive Documentation
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

### Key Endpoints
- `POST /api/auth/login` - User authentication
- `POST /api/upload` - Upload VCF files
- `GET /api/analysis/{id}` - Get analysis status
- `POST /api/pipeline/parse-vcf` - Parse VCF files
- `POST /api/pipeline/match-drugs` - Match drug interactions

## 🔄 CI/CD Pipeline

### GitHub Actions (Example)
```yaml
# .github/workflows/deploy.yml
name: Deploy to Production

on:
  push:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run tests
        run: make test

  deploy:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - name: Deploy to production
        run: make prod-deploy
```

## 📞 Support

### Getting Help
1. Check the troubleshooting section
2. Review application logs
3. Check system metrics in Grafana
4. Create an issue in the repository

### Performance Optimization
- Monitor Grafana dashboards
- Optimize database queries
- Scale services based on load
- Use CDN for static assets

---

## 🎉 Production Checklist

- [ ] Environment variables configured
- [ ] SSL certificates installed
- [ ] Database migrations completed
- [ ] Monitoring setup (Grafana/Prometheus)
- [ ] Backup strategy implemented
- [ ] Security scan completed
- [ ] Load testing performed
- [ ] Documentation updated
- [ ] Team training completed

**Your production-ready pharmacogenomics analysis platform is now ready!** 🧬💊

For additional support or custom deployment needs, please refer to the documentation or contact the development team.