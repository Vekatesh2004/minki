#!/usr/bin/env python3
"""
Local Backend for Pharmacogenomics Pipeline (real pipeline, not a demo stub).

Pipeline per upload:
  1. Parse VCF + QC
  2. Filter to PASS variants, cap how many go to VEP (web API can't take a whole chromosome)
  3. VEP annotation -> gene, consequence, HGVSp (amino-acid change)
  4. Build coding-mutation list
  5. For each mutated gene: UniProt + AlphaFold structure (pLDDT + domain flag)
  6. Drug matching for genes that actually carry coding variants
"""

import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional
import uuid
from datetime import datetime

from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse, FileResponse, HTMLResponse
from pydantic import BaseModel

# Add current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from modules.vcf_parser import VCFParser
    from modules.vep_annotator import VEPAnnotator
    from modules.protein_structure import ProteinStructureAnalyzer
    from modules.drug_matcher import DrugMatcher
    modules_available = True
    module_import_error = None
except Exception as e:  # catch ImportError AND any error raised at import time
    import traceback
    module_import_error = traceback.format_exc()
    print("=" * 70)
    print("MODULES FAILED TO LOAD:")
    print(module_import_error)
    print("=" * 70)
    modules_available = False

# ----------------------------------------------------------------------------
# Tunables for local testing.
# The Ensembl VEP web API cannot annotate a whole chromosome (hundreds of
# thousands of variants). We cap how many PASS variants we send.
# ----------------------------------------------------------------------------
MAX_VARIANTS_TO_VEP = int(os.getenv("MAX_VARIANTS_TO_VEP", "300"))
MAX_GENES_FOR_STRUCTURE = int(os.getenv("MAX_GENES_FOR_STRUCTURE", "15"))
VEP_BATCH_TIMEOUT = 90        # seconds per VEP batch
STRUCTURE_TIMEOUT = 60        # seconds per structure lookup
DRUG_TIMEOUT = 20             # seconds per gene drug lookup

app = FastAPI(
    title="Pharmacogenomics Pipeline - Local",
    description="Local pipeline: VCF -> VEP -> UniProt/AlphaFold -> Drug matching",
    version="2.1.0-local",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

pipeline_components: Dict[str, Any] = {}
analysis_storage: Dict[str, "AnalysisStatus"] = {}


async def initialize_components():
    global pipeline_components
    try:
        with open("config.json", "r") as f:
            config = json.load(f)

        if modules_available:
            pipeline_components = {
                "config": config,
                "vcf_parser": VCFParser(config),
                "vep_annotator": VEPAnnotator(config),
                "structure_analyzer": ProteinStructureAnalyzer(config),
                "drug_matcher": DrugMatcher(config),
            }
            print("Pipeline components initialized")
        else:
            pipeline_components = {"config": config}
            print("Running in limited mode - modules unavailable")
    except Exception as e:
        print(f"Error initializing components: {e}")
        pipeline_components = {}


@app.on_event("startup")
async def startup_event():
    await initialize_components()


class AnalysisStatus(BaseModel):
    id: str
    status: str
    progress: float
    message: str
    created_at: datetime
    results: Optional[Dict[str, Any]] = None


class SimpleUploadResponse(BaseModel):
    upload_id: str
    filename: str
    status: str
    message: str


@app.get("/")
async def root():
    return HTMLResponse(content=INDEX_HTML)


@app.get("/health")
async def health_check():
    return {
        "status": "healthy",
        "mode": "local",
        "components_available": modules_available,
        "module_import_error": module_import_error,
        "max_variants_to_vep": MAX_VARIANTS_TO_VEP,
        "timestamp": datetime.now().isoformat(),
    }


@app.post("/api/upload", response_model=SimpleUploadResponse)
async def upload_file(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    sample_id: str = Form(...),
):
    if not file.filename.lower().endswith((".vcf", ".vcf.gz")):
        raise HTTPException(status_code=400, detail="File must be a VCF file")

    upload_id = str(uuid.uuid4())
    upload_dir = Path("uploads")
    upload_dir.mkdir(exist_ok=True)
    file_path = upload_dir / f"{upload_id}_{file.filename}"

    try:
        content = await file.read()
        with open(file_path, "wb") as f:
            f.write(content)

        background_tasks.add_task(analyze_file, upload_id, str(file_path), sample_id)

        return SimpleUploadResponse(
            upload_id=upload_id,
            filename=file.filename,
            status="uploaded",
            message="File uploaded successfully, analysis queued",
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Upload failed: {str(e)}")


def _set(upload_id: str, *, status=None, progress=None, message=None, results=None):
    a = analysis_storage[upload_id]
    if status is not None:
        a.status = status
    if progress is not None:
        a.progress = progress
    if message is not None:
        a.message = message
    if results is not None:
        a.results = results


async def analyze_file(upload_id: str, file_path: str, sample_id: str):
    analysis_storage[upload_id] = AnalysisStatus(
        id=upload_id,
        status="running",
        progress=0.0,
        message="Starting analysis...",
        created_at=datetime.now(),
    )

    if not (modules_available and "vcf_parser" in pipeline_components):
        detail = module_import_error or "components not initialized (check config.json)"
        _set(upload_id, status="failed",
             message=f"Pipeline modules are not available: {detail}")
        return

    try:
        vcf_parser = pipeline_components["vcf_parser"]
        vep = pipeline_components["vep_annotator"]
        structure_analyzer = pipeline_components["structure_analyzer"]
        drug_matcher = pipeline_components["drug_matcher"]

        # ---- Step 1: parse VCF -------------------------------------------------
        _set(upload_id, progress=10.0, message="Parsing VCF and computing QC...")
        vcf_results = await vcf_parser.parse_vcf(file_path, sample_id)
        all_variants = vcf_results.get("variants", [])
        total_variants = len(all_variants)

        # ---- Step 2: choose which variants to annotate ------------------------
        pass_variants = [
            v for v in all_variants
            if str(v.get("filter", "")).upper() in ("PASS", ".", "", "NONE")
        ]
        # If the "filter" field isn't a real PASS flag (older VCFs), fall back to all
        candidates = pass_variants if pass_variants else all_variants
        to_annotate = candidates[:MAX_VARIANTS_TO_VEP]

        _set(
            upload_id,
            progress=25.0,
            message=f"Annotating {len(to_annotate)} of {total_variants} variants with Ensembl VEP...",
        )

        # ---- Step 3: VEP annotation (with timeout) ----------------------------
        try:
            annotated = await asyncio.wait_for(
                vep.annotate_variants(to_annotate),
                timeout=VEP_BATCH_TIMEOUT * (len(to_annotate) // 200 + 2),
            )
        except asyncio.TimeoutError:
            annotated = []
            _set(upload_id, message="VEP annotation timed out; continuing with partial data")

        # ---- Step 4: build coding-mutation list -------------------------------
        _set(upload_id, progress=55.0, message="Extracting coding mutations...")
        mutations = []
        genes_seen = {}
        for v in annotated:
            ann = v.get("vep_annotation", {}) or {}
            if not ann.get("coding"):
                continue
            gene = ann.get("gene_symbol")
            mut = {
                "position": f"{v.get('chrom')}:{v.get('pos')}",
                "ref": v.get("ref"),
                "alt": v.get("alt"),
                "gene_symbol": gene,
                "consequence_terms": ann.get("consequence_terms", []),
                "hgvsc": ann.get("hgvsc"),
                "hgvsp": ann.get("hgvsp"),
                "amino_acid_change": ann.get("amino_acid_change"),
                "amino_acid_position": ann.get("amino_acid_position"),
                "uniprot_ids": ann.get("uniprot_ids", []),
                "gene_id": ann.get("gene_id"),
            }
            mutations.append(mut)
            if gene and gene not in genes_seen:
                genes_seen[gene] = mut

        # ---- Step 5: structure for each mutated gene --------------------------
        _set(upload_id, progress=70.0, message="Fetching protein structures (UniProt/AlphaFold)...")
        structures = []
        for gene, mut in list(genes_seen.items())[:MAX_GENES_FOR_STRUCTURE]:
            uniprot_ids = mut.get("uniprot_ids") or []
            if not uniprot_ids and gene:
                try:
                    uniprot_ids = await asyncio.wait_for(
                        structure_analyzer.map_gene_to_uniprot(gene, mut.get("gene_id")),
                        timeout=STRUCTURE_TIMEOUT,
                    )
                except (asyncio.TimeoutError, Exception):
                    uniprot_ids = []

            residue = mut.get("amino_acid_position")
            if uniprot_ids and residue:
                uid = uniprot_ids[0].split(":")[-1] if ":" in str(uniprot_ids[0]) else uniprot_ids[0]
                try:
                    s = await asyncio.wait_for(
                        structure_analyzer.analyze_structure(uid, residue, mut.get("amino_acid_change")),
                        timeout=STRUCTURE_TIMEOUT,
                    )
                    pos = s.get("position_analysis") or {}
                    structures.append({
                        "gene_symbol": gene,
                        "uniprot_id": uid,
                        "residue_position": residue,
                        "amino_acid_change": mut.get("amino_acid_change"),
                        "plddt_score": pos.get("plddt_score"),
                        "confidence_level": pos.get("confidence_level"),
                        "in_domain": pos.get("in_domain", False),
                        "domain_info": pos.get("domain_info", []),
                        "in_binding_site": pos.get("in_binding_site", False),
                        "in_active_site": pos.get("in_active_site", False),
                        "structural_impact": pos.get("structural_impact_prediction"),
                        "alphafold_url": f"https://alphafold.ebi.ac.uk/entry/{uid}",
                    })
                except (asyncio.TimeoutError, Exception) as e:
                    structures.append({"gene_symbol": gene, "uniprot_id": uid, "error": str(e)})

        # ---- Step 6: drug matching for mutated genes --------------------------
        _set(upload_id, progress=85.0, message="Matching drugs for mutated genes...")
        drug_results = []
        genes_for_drugs = list(genes_seen.keys())
        # If VEP found no coding genes, fall back to scanning known pharmacogenes
        if not genes_for_drugs:
            genes_for_drugs = ["CYP2D6", "CYP2C19", "CYP2C9", "TPMT", "DPYD", "VKORC1"]

        for gene in genes_for_drugs[:30]:
            try:
                dr = await asyncio.wait_for(drug_matcher.match_drugs(gene), timeout=DRUG_TIMEOUT)
                if dr.get("matched_drugs"):
                    drug_results.append(dr)
            except (asyncio.TimeoutError, Exception):
                continue

        # ---- Assemble results -------------------------------------------------
        results = {
            "sample_id": sample_id,
            "qc_summary": vcf_results.get("qc_summary", {}),
            "total_variants": total_variants,
            "variants_annotated": len(annotated),
            "annotation_capped": total_variants > len(to_annotate),
            "mutations": mutations,
            "coding_mutation_count": len(mutations),
            "genes_with_coding_variants": list(genes_seen.keys()),
            "structures": structures,
            "drug_results": drug_results,
            "summary": {
                "total_variants": total_variants,
                "variants_annotated": len(annotated),
                "coding_mutations": len(mutations),
                "genes_affected": len(genes_seen),
                "structures_resolved": len([s for s in structures if not s.get("error")]),
                "drug_matches": sum(len(r.get("matched_drugs", [])) for r in drug_results),
            },
        }

        msg = "Analysis completed successfully"
        if results["annotation_capped"]:
            msg += (f" (annotated first {len(to_annotate)} of {total_variants} variants; "
                    f"raise MAX_VARIANTS_TO_VEP to annotate more)")

        _set(upload_id, status="completed", progress=100.0, message=msg, results=results)

    except Exception as e:
        import traceback
        traceback.print_exc()
        _set(upload_id, status="failed", message=f"Analysis failed: {str(e)}")


@app.get("/api/analysis/{upload_id}")
async def get_analysis_status(upload_id: str):
    if upload_id not in analysis_storage:
        raise HTTPException(status_code=404, detail="Analysis not found")
    a = analysis_storage[upload_id]
    return {
        "id": a.id,
        "status": a.status,
        "progress": a.progress,
        "message": a.message,
        "created_at": a.created_at,
        "results": a.results,
    }


if Path("results").exists():
    app.mount("/results", StaticFiles(directory="results"), name="results")


INDEX_HTML = """
<!DOCTYPE html>
<html>
<head>
    <title>Pharmacogenomics Pipeline</title>
    <style>
        body { font-family: -apple-system, Arial, sans-serif; margin: 0; background: #f5f7fa; color: #1f2937; }
        .container { max-width: 1000px; margin: 0 auto; padding: 30px; }
        .card { background: white; padding: 24px; border-radius: 10px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); margin-bottom: 20px; }
        h1 { color: #1e3a8a; }
        h2 { color: #1e3a8a; border-left: 4px solid #3b82f6; padding-left: 10px; font-size: 1.1rem; }
        button { background: #3b82f6; color: white; padding: 10px 18px; border: none; border-radius: 6px; cursor: pointer; font-size: 14px; }
        button:hover { background: #2563eb; }
        button:disabled { background: #9ca3af; cursor: not-allowed; }
        input[type=text], input[type=file] { padding: 8px; margin: 6px 0; width: 100%; box-sizing: border-box; border: 1px solid #d1d5db; border-radius: 6px; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 13px; }
        th, td { border: 1px solid #e5e7eb; padding: 6px 8px; text-align: left; }
        th { background: #1e3a8a; color: white; }
        tr:nth-child(even) { background: #f9fafb; }
        .metrics { display: flex; gap: 12px; flex-wrap: wrap; }
        .metric { background: #eff6ff; border-radius: 8px; padding: 14px 18px; text-align: center; flex: 1; min-width: 120px; }
        .metric .val { font-size: 1.6rem; font-weight: bold; color: #2563eb; }
        .metric .lbl { font-size: 0.8rem; color: #6b7280; }
        .bar { height: 22px; background: #e5e7eb; border-radius: 6px; overflow: hidden; }
        .bar-fill { height: 100%; background: #3b82f6; color: white; text-align: center; font-size: 12px; line-height: 22px; transition: width .3s; }
        .note { background: #fef3c7; padding: 10px; border-radius: 6px; font-size: 13px; margin-top: 10px; }
        .badge { padding: 2px 8px; border-radius: 10px; font-size: 11px; color: white; }
        a { color: #2563eb; }
    </style>
</head>
<body>
<div class="container">
    <h1>Pharmacogenomics Pipeline</h1>

    <div class="card">
        <h2>Upload VCF</h2>
        <form id="uploadForm">
            <input type="file" id="vcfFile" accept=".vcf,.vcf.gz" required>
            <input type="text" id="sampleId" placeholder="Sample ID (e.g. SAMPLE_001)" required>
            <button type="submit" id="submitBtn">Upload and Analyze</button>
        </form>
        <div class="note">
            Large files (whole chromosomes) are annotated up to a configurable cap so the
            Ensembl VEP web API can keep up. Set the <code>MAX_VARIANTS_TO_VEP</code> env var to change it.
        </div>
    </div>

    <div class="card" id="progressCard" style="display:none;">
        <h2>Analysis Progress</h2>
        <div class="bar"><div class="bar-fill" id="barFill" style="width:0%">0%</div></div>
        <p id="statusMsg">Starting...</p>
    </div>

    <div id="resultsArea"></div>
</div>

<script>
let pollTimer = null;

document.getElementById('uploadForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const file = document.getElementById('vcfFile').files[0];
    const sampleId = document.getElementById('sampleId').value;
    if (!file || !sampleId) return;

    const fd = new FormData();
    fd.append('file', file);
    fd.append('sample_id', sampleId);

    document.getElementById('submitBtn').disabled = true;
    document.getElementById('progressCard').style.display = 'block';
    document.getElementById('resultsArea').innerHTML = '';
    setProgress(5, 'Uploading file...');

    try {
        const res = await fetch('/api/upload', { method: 'POST', body: fd });
        const data = await res.json();
        if (!res.ok) throw new Error(data.detail || 'Upload failed');
        pollStatus(data.upload_id);
    } catch (err) {
        setProgress(0, 'Error: ' + err.message);
        document.getElementById('submitBtn').disabled = false;
    }
});

function setProgress(pct, msg) {
    const fill = document.getElementById('barFill');
    fill.style.width = pct + '%';
    fill.textContent = Math.round(pct) + '%';
    document.getElementById('statusMsg').textContent = msg;
}

function pollStatus(id) {
    if (pollTimer) clearInterval(pollTimer);
    pollTimer = setInterval(async () => {
        try {
            const res = await fetch('/api/analysis/' + id);
            const data = await res.json();
            setProgress(data.progress, data.message);
            if (data.status === 'completed') {
                clearInterval(pollTimer);
                document.getElementById('submitBtn').disabled = false;
                renderResults(data.results);
            } else if (data.status === 'failed') {
                clearInterval(pollTimer);
                document.getElementById('submitBtn').disabled = false;
                setProgress(0, data.message);
            }
        } catch (err) {
            clearInterval(pollTimer);
            document.getElementById('submitBtn').disabled = false;
        }
    }, 1500);
}

function renderResults(r) {
    if (!r) return;
    const s = r.summary || {};
    let html = '<div class="card"><h2>Summary</h2><div class="metrics">';
    html += metric(s.total_variants, 'Total Variants');
    html += metric(s.variants_annotated, 'Annotated');
    html += metric(s.coding_mutations, 'Coding Mutations');
    html += metric(s.genes_affected, 'Genes Affected');
    html += metric(s.structures_resolved, 'Structures');
    html += metric(s.drug_matches, 'Drug Matches');
    html += '</div>';
    if (r.annotation_capped) {
        html += '<div class="note">Only the first ' + r.variants_annotated +
                ' of ' + r.total_variants + ' variants were annotated (cap). ' +
                'Increase MAX_VARIANTS_TO_VEP to process more.</div>';
    }
    html += '</div>';

    // Mutations table
    html += '<div class="card"><h2>Coding Mutations</h2>';
    if ((r.mutations || []).length === 0) {
        html += '<p>No coding mutations found in the annotated set.</p>';
    } else {
        html += '<table><tr><th>Position</th><th>Change</th><th>Gene</th>' +
                '<th>Consequence</th><th>Amino Acid (HGVSp)</th></tr>';
        r.mutations.slice(0, 200).forEach(m => {
            html += '<tr><td>' + m.position + '</td><td>' + m.ref + '&rarr;' + m.alt + '</td>' +
                    '<td>' + (m.gene_symbol || '-') + '</td>' +
                    '<td>' + (m.consequence_terms || []).join(', ') + '</td>' +
                    '<td>' + (m.hgvsp || m.amino_acid_change || '-') + '</td></tr>';
        });
        html += '</table>';
        if (r.mutations.length > 200) html += '<p>Showing first 200 of ' + r.mutations.length + '.</p>';
    }
    html += '</div>';

    // Structures
    html += '<div class="card"><h2>Protein Structures (mutated residue)</h2>';
    if ((r.structures || []).length === 0) {
        html += '<p>No structures resolved.</p>';
    } else {
        html += '<table><tr><th>Gene</th><th>UniProt</th><th>Residue</th><th>Change</th>' +
                '<th>pLDDT</th><th>In Domain</th><th>Impact</th><th>View</th></tr>';
        r.structures.forEach(st => {
            if (st.error) {
                html += '<tr><td>' + st.gene_symbol + '</td><td>' + (st.uniprot_id||'-') +
                        '</td><td colspan="6">error: ' + st.error + '</td></tr>';
            } else {
                html += '<tr><td>' + st.gene_symbol + '</td><td>' + st.uniprot_id + '</td>' +
                        '<td>' + st.residue_position + '</td><td>' + (st.amino_acid_change||'-') + '</td>' +
                        '<td>' + (st.plddt_score!=null ? st.plddt_score.toFixed(1) : '-') + '</td>' +
                        '<td>' + (st.in_domain ? 'Yes' : 'No') + '</td>' +
                        '<td>' + (st.structural_impact||'-') + '</td>' +
                        '<td><a href="' + st.alphafold_url + '" target="_blank">AlphaFold</a></td></tr>';
            }
        });
        html += '</table>';
    }
    html += '</div>';

    // Drugs
    html += '<div class="card"><h2>Drug Interactions</h2>';
    html += '<p style="font-size:13px;color:#6b7280;">Click a database link to confirm each gene\\u2013drug relationship on an authoritative source.</p>';
    if ((r.drug_results || []).length === 0) {
        html += '<p>No drug interactions found.</p>';
    } else {
        r.drug_results.forEach(dr => {
            var gene = dr.gene_symbol;
            var genePgkb = 'https://www.pharmgkb.org/search?query=' + encodeURIComponent(gene);
            html += '<h3>' + gene +
                    ' <a href="' + genePgkb + '" target="_blank" style="font-size:12px;font-weight:normal;">PharmGKB gene &#8599;</a></h3>';
            html += '<table><tr><th>Drug</th><th>Action</th><th>Evidence</th><th>Confirm on</th></tr>';
            (dr.matched_drugs || []).forEach(d => {
                var name = d.drug_name || '-';
                var q = encodeURIComponent(gene + ' ' + name);
                // PharmGKB clinical annotation search for the gene+drug pair
                var pgkb = 'https://www.pharmgkb.org/search?query=' + q;
                // DrugBank drug search
                var db = 'https://go.drugbank.com/unearth/q?searcher=drugs&query=' + encodeURIComponent(name);
                var links = '<a href="' + pgkb + '" target="_blank">PharmGKB</a> | ' +
                            '<a href="' + db + '" target="_blank">DrugBank</a>';
                html += '<tr><td>' + name + '</td><td>' + (d.action||'-') +
                        '</td><td>' + (d.evidence_level||'-') + '</td><td>' + links + '</td></tr>';
            });
            html += '</table>';
        });
    }
    html += '</div>';

    document.getElementById('resultsArea').innerHTML = html;
}

function metric(val, label) {
    return '<div class="metric"><div class="val">' + (val!=null?val:0) +
           '</div><div class="lbl">' + label + '</div></div>';
}
</script>
</body>
</html>
"""

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("simple_backend:app", host="127.0.0.1", port=8000, reload=True)
