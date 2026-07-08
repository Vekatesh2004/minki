# 🧬 Pharmacogenomics Pipeline - Quick Start Guide

## 🚀 **Get Started in 3 Simple Steps**

### **Step 1: Run the Local Setup**
```bash
# Option A: Automated setup (recommended)
./install_local.sh

# Option B: Manual setup
python3 run_local.py
```

### **Step 2: Start the Server**
```bash
python3 simple_backend.py
```

### **Step 3: Open Your Browser**
Visit: **http://localhost:8000**

---

## 🎯 **What You Get**

### **✨ Web Interface Features:**
- 📁 **File Upload**: Drag & drop VCF files
- 🔍 **Real-time Analysis**: Watch progress in real-time
- 📊 **Interactive Results**: View analysis results instantly
- 🧬 **Pharmacogenomics**: Drug interaction predictions
- 📈 **Quality Control**: Comprehensive QC metrics
- 🎮 **Demo Mode**: Try with sample data

### **🔗 Key URLs:**
- **Main Interface**: http://localhost:8000
- **API Documentation**: http://localhost:8000/docs
- **Demo Analysis**: http://localhost:8000/demo
- **Health Check**: http://localhost:8000/health

---

## 📋 **Quick Test Workflow**

### **1. Upload a VCF File**
- Go to http://localhost:8000
- Click "Choose File" and select a VCF file
- Enter a sample ID (e.g., "TEST_SAMPLE")
- Click "Upload and Analyze"

### **2. Try the Demo**
- Visit http://localhost:8000/demo
- Uses built-in sample VCF with known pharmacogenomic variants
- See instant results for CYP2D6, CYP2C19, CYP2C9 genes

### **3. Explore API**
- Visit http://localhost:8000/docs
- Interactive API documentation
- Test endpoints directly in browser

---

## 🛠️ **What's Running**

### **Backend Components:**
- ✅ **FastAPI Server** (http://localhost:8000)
- ✅ **VCF Parser** (quality control & variant analysis)
- ✅ **Drug Matcher** (pharmacogenomic database)
- ✅ **File Upload** (secure file handling)
- ✅ **Background Processing** (async analysis)

### **Features Working:**
- 📊 VCF file parsing and QC metrics
- 🧬 Pharmacogene variant analysis
- 💊 Drug interaction matching (CYP2D6, CYP2C19, CYP2C9, TPMT, DPYD)
- 📈 Real-time progress tracking
- 📋 JSON and HTML result export

---

## 🔧 **Troubleshooting**

### **If Setup Fails:**
```bash
# Install dependencies manually
pip install fastapi uvicorn python-multipart requests pandas

# Run setup again
python3 run_local.py
```

### **If Server Won't Start:**
```bash
# Check if port 8000 is free
netstat -an | grep 8000

# Try different port
uvicorn simple_backend:app --host 127.0.0.1 --port 8001
```

### **Common Issues:**

**1. "Module not found" errors:**
- Run: `pip install fastapi uvicorn`
- Make sure you're in the correct directory

**2. "Port already in use":**
- Stop other services on port 8000
- Or change port in `simple_backend.py`

**3. "Permission denied":**
- Run: `chmod +x install_local.sh`
- Or use: `python3 run_local.py` directly

---

## 📁 **File Structure**
```
minki/
├── simple_backend.py          # 🚀 Main server (auto-generated)
├── run_local.py               # 🔧 Setup script
├── install_local.sh           # 📦 Installation script
├── config.json                # ⚙️ Configuration
├── modules/                   # 🧬 Core pipeline
├── examples/                  # 📄 Sample VCF files
├── uploads/                   # 📁 User uploads
└── results/                   # 📊 Analysis results
```

---

## 🎮 **Example Usage**

### **1. Test with Sample Data**
```bash
# Demo analysis runs automatically
curl http://localhost:8000/demo
```

### **2. Upload Your Own VCF**
```bash
# Using curl
curl -X POST "http://localhost:8000/api/upload" \
     -F "file=@your_file.vcf" \
     -F "sample_id=YOUR_SAMPLE"
```

### **3. Check Analysis Status**
```bash
# Get results
curl http://localhost:8000/api/analysis/UPLOAD_ID
```

---

## 🎯 **Next Steps**

### **Once Testing is Complete:**
1. **Production Setup**: Follow `README_PRODUCTION.md`
2. **Docker Deployment**: Use `docker-compose.yml`
3. **Add Features**: Extend with React frontend
4. **Scale Up**: Add database and caching

### **Current Limitations (Local Mode):**
- ⚠️ No persistent database (in-memory storage)
- ⚠️ Limited to basic pharmacogenes
- ⚠️ No real-time WebSocket updates
- ⚠️ Single-user mode only

### **Full Production Features Available:**
- 🔐 User authentication & authorization
- 💾 PostgreSQL database with migrations
- 🚀 Redis caching and task queues
- 📊 Prometheus monitoring & Grafana dashboards
- 🔧 Horizontal scaling with Docker
- ⚡ React frontend with real-time updates

---

## 🆘 **Need Help?**

1. **Check the logs**: Look at terminal output
2. **Verify setup**: Run `python3 run_local.py` again
3. **Test connectivity**: Visit http://localhost:8000/health
4. **Review docs**: Check http://localhost:8000/docs

---

## 🎉 **You're Ready!**

Your pharmacogenomics analysis pipeline is now running locally! 

**Start analyzing genetic variants for drug interactions right away!** 🧬💊

---

*For production deployment, advanced features, and scaling options, see `README_PRODUCTION.md`*