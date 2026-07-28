# Minki / HelixAI — active pharmacogenomics application (simple_backend.py)
# Single-image build that runs the local FastAPI product end to end.
FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PORT=8000 \
    HOST=0.0.0.0

# System build deps (cyvcf2/scipy/statsmodels wheels usually cover this, but
# these headers make source builds reliable across architectures).
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ make \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies first for better layer caching.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code and assets.
COPY simple_backend.py config.json ./
COPY modules/ ./modules/
COPY scripts/ ./scripts/
COPY examples/ ./examples/

# Generate the synthetic demonstration datasets so the ClinVar / common-risk
# panels populate for the bundled example. These are clearly-labelled synthetic
# demo data, not real clinical evidence.
RUN python scripts/make_diabetes_demo_data.py

# Runtime working directories.
RUN mkdir -p /app/uploads /app/results

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=10s --retries=3 \
  CMD curl -fsS http://localhost:8000/health || exit 1

CMD ["python", "simple_backend.py"]
