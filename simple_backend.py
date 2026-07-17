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
    from modules import diabetes_kb
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
BASE_DIR = Path(__file__).resolve().parent
EXAMPLE_VCF_PATH = BASE_DIR / "examples" / "sample_pharmacogenomics.vcf"
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


@app.post("/api/example", response_model=SimpleUploadResponse)
async def analyze_example(background_tasks: BackgroundTasks):
    """Queue analysis of the bundled example VCF without requiring an upload."""
    if not EXAMPLE_VCF_PATH.is_file():
        raise HTTPException(status_code=503, detail="Bundled example VCF is unavailable")

    upload_id = str(uuid.uuid4())
    background_tasks.add_task(
        analyze_file,
        upload_id,
        str(EXAMPLE_VCF_PATH),
        "SAMPLE_PHARMACOGENOMICS_EXAMPLE",
    )
    return SimpleUploadResponse(
        upload_id=upload_id,
        filename=EXAMPLE_VCF_PATH.name,
        status="queued",
        message="Example analysis queued",
    )


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

        # ---- Mutation-type breakdown (for charts) -----------------------------
        mutation_type_counts = {}
        for m in mutations:
            for term in (m.get("consequence_terms") or ["unknown"]):
                mutation_type_counts[term] = mutation_type_counts.get(term, 0) + 1

        # ---- Diabetes interpretation ------------------------------------------
        diabetes_findings = []
        for m in mutations:
            gene = m.get("gene_symbol")
            kb = diabetes_kb.lookup(gene)
            if kb:
                diabetes_findings.append({
                    "gene_symbol": gene,
                    "position": m.get("position"),
                    "change": f"{m.get('ref')}>{m.get('alt')}",
                    "consequence_terms": m.get("consequence_terms", []),
                    "amino_acid_change": m.get("hgvsp") or m.get("amino_acid_change"),
                    "category": kb["category"],
                    "role": kb["role"],
                    "significance": kb["significance"],
                    "treatment": kb["treatment"],
                    "pharmgkb_url": kb["pharmgkb_url"],
                    "omim_url": kb["omim_url"],
                    "gwas_url": kb["gwas_url"],
                })

        # ---- Assemble results -------------------------------------------------
        results = {
            "sample_id": sample_id,
            "qc_summary": vcf_results.get("qc_summary", {}),
            "total_variants": total_variants,
            "variants_annotated": len(annotated),
            "annotation_capped": total_variants > len(to_annotate),
            "mutations": mutations,
            "mutation_type_counts": mutation_type_counts,
            "coding_mutation_count": len(mutations),
            "genes_with_coding_variants": list(genes_seen.keys()),
            "structures": structures,
            "drug_results": drug_results,
            "diabetes_findings": diabetes_findings,
            "summary": {
                "total_variants": total_variants,
                "variants_annotated": len(annotated),
                "coding_mutations": len(mutations),
                "genes_affected": len(genes_seen),
                "structures_resolved": len([s for s in structures if not s.get("error")]),
                "drug_matches": sum(len(r.get("matched_drugs", [])) for r in drug_results),
                "diabetes_genes_found": len(diabetes_findings),
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


@app.get("/api/analysis/{upload_id}/mutations.csv")
async def download_mutations_csv(upload_id: str):
    """Download the identified mutations as a CSV file."""
    if upload_id not in analysis_storage:
        raise HTTPException(status_code=404, detail="Analysis not found")
    a = analysis_storage[upload_id]
    if not a.results:
        raise HTTPException(status_code=400, detail="Analysis not complete")

    import csv, io
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([
        "position", "ref", "alt", "gene_symbol", "consequence_terms",
        "hgvsc", "hgvsp", "amino_acid_change", "diabetes_gene",
    ])
    for m in a.results.get("mutations", []):
        writer.writerow([
            m.get("position"), m.get("ref"), m.get("alt"), m.get("gene_symbol"),
            "|".join(m.get("consequence_terms") or []),
            m.get("hgvsc"), m.get("hgvsp"), m.get("amino_acid_change"),
            "yes" if diabetes_kb.is_diabetes_gene(m.get("gene_symbol")) else "no",
        ])

    from fastapi.responses import Response
    fname = f"{a.results.get('sample_id', 'sample')}_mutations.csv"
    return Response(
        content=buf.getvalue(),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={fname}"},
    )


if Path("results").exists():
    app.mount("/results", StaticFiles(directory="results"), name="results")


INDEX_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="theme-color" content="#f4f8ff">
    <meta name="description" content="AI-assisted pharmacogenomics, molecular structure, and precision medicine analysis.">
    <title>HelixAI | Pharmacogenomics Intelligence</title>
    <style>
        :root {
            --ink: #10213b; --muted: #61718a; --blue: #2764eb; --blue-2: #4f8cff;
            --teal: #0d9488; --violet: #7357d9; --line: rgba(93, 125, 176, .19);
            --glass: rgba(255, 255, 255, .74); --glass-strong: rgba(255, 255, 255, .9);
            --shadow: 0 22px 60px rgba(37, 73, 132, .11); --radius: 22px;
            --parallax-x: 0px; --parallax-y: 0px;
        }
        * { box-sizing: border-box; }
        html { scroll-behavior: smooth; }
        body {
            margin: 0; min-width: 320px; color: var(--ink); overflow-x: hidden;
            font-family: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            background:
                radial-gradient(circle at 8% 8%, rgba(90, 140, 255, .14), transparent 30rem),
                radial-gradient(circle at 90% 22%, rgba(22, 186, 173, .1), transparent 26rem),
                linear-gradient(145deg, #f9fbff 0%, #f1f6fd 48%, #f8fbff 100%);
            line-height: 1.55;
        }
        body::before { content: ""; position: fixed; inset: 0; z-index: -3; opacity: .32; background-image: linear-gradient(rgba(59, 100, 168, .05) 1px, transparent 1px), linear-gradient(90deg, rgba(59, 100, 168, .05) 1px, transparent 1px); background-size: 44px 44px; mask-image: linear-gradient(to bottom, black, transparent 78%); }
        button, input { font: inherit; }
        a { color: var(--blue); text-decoration: none; font-weight: 650; transition: color .2s ease, opacity .2s ease; }
        a:hover { color: #1749b8; }
        :focus-visible { outline: 3px solid rgba(39, 100, 235, .35); outline-offset: 3px; }
        .skip-link { position: fixed; z-index: 100; left: 16px; top: -60px; padding: 10px 15px; border-radius: 10px; background: #fff; box-shadow: var(--shadow); }
        .skip-link:focus { top: 12px; }
        .science-bg { position: fixed; inset: 0; z-index: -2; pointer-events: none; overflow: hidden; color: #316ed9; }
        .science-bg svg { width: 100%; height: 100%; opacity: .095; transform: translate3d(var(--parallax-x), var(--parallax-y), 0); transition: transform .25s ease-out; }
        .float-a { transform-origin: 18% 25%; animation: floatA 18s ease-in-out infinite; }
        .float-b { transform-origin: 82% 22%; animation: floatB 24s ease-in-out infinite; }
        .float-c { transform-origin: 78% 78%; animation: floatA 21s ease-in-out -7s infinite reverse; }
        .particle { animation: particle 8s ease-in-out infinite; }
        .particle:nth-child(2n) { animation-delay: -3s; }
        @keyframes floatA { 0%,100% { transform: translate3d(0,0,0) rotate(0); } 50% { transform: translate3d(18px,-18px,0) rotate(5deg); } }
        @keyframes floatB { 0%,100% { transform: translate3d(0,0,0) rotate(0) scale(1); } 50% { transform: translate3d(-16px,22px,0) rotate(-7deg) scale(1.04); } }
        @keyframes particle { 0%,100% { opacity: .25; transform: translateY(0); } 50% { opacity: .85; transform: translateY(-16px); } }
        .topbar { position: relative; z-index: 10; border-bottom: 1px solid rgba(113, 143, 190, .16); background: rgba(250, 252, 255, .7); backdrop-filter: blur(18px); -webkit-backdrop-filter: blur(18px); }
        .topbar-inner { width: min(1180px, calc(100% - 40px)); min-height: 72px; margin: auto; display: flex; align-items: center; justify-content: space-between; gap: 20px; }
        .brand { display: inline-flex; align-items: center; gap: 11px; color: var(--ink); }
        .brand-mark { width: 38px; height: 38px; display: grid; place-items: center; border: 1px solid rgba(75, 126, 216, .22); border-radius: 12px; color: white; background: linear-gradient(145deg, var(--blue), #22a6a0); box-shadow: 0 9px 24px rgba(39, 100, 235, .24); }
        .brand-mark svg { width: 22px; height: 22px; }
        .brand-copy { display: grid; line-height: 1.12; }
        .brand-copy strong { font-size: .98rem; letter-spacing: -.01em; }
        .brand-copy span { margin-top: 3px; color: var(--muted); font-size: .7rem; font-weight: 650; letter-spacing: .07em; text-transform: uppercase; }
        .platform-status { display: inline-flex; align-items: center; gap: 8px; padding: 8px 12px; color: #31605e; border: 1px solid rgba(13, 148, 136, .18); border-radius: 999px; background: rgba(236, 253, 250, .72); font-size: .76rem; font-weight: 700; }
        .status-dot { width: 7px; height: 7px; border-radius: 50%; background: #13a89e; box-shadow: 0 0 0 5px rgba(19, 168, 158, .1); animation: statusPulse 2.4s ease-out infinite; }
        @keyframes statusPulse { 70%,100% { box-shadow: 0 0 0 8px rgba(19,168,158,0); } }
        .container { position: relative; width: min(1180px, calc(100% - 40px)); margin: 0 auto; padding: 70px 0 80px; }
        .hero { display: grid; grid-template-columns: minmax(0, 1.35fr) minmax(300px, .65fr); gap: 54px; align-items: center; margin-bottom: 42px; animation: enter .75s both; }
        .eyebrow { display: inline-flex; align-items: center; gap: 9px; margin-bottom: 18px; color: #2761ca; font-size: .75rem; font-weight: 800; letter-spacing: .12em; text-transform: uppercase; }
        .eyebrow::before { content: ""; width: 28px; height: 1px; background: linear-gradient(90deg, var(--blue), var(--teal)); }
        h1 { max-width: 780px; margin: 0; color: #10284d; font-size: clamp(2.35rem, 5vw, 4.55rem); line-height: .99; letter-spacing: -.055em; }
        .gradient-text { color: transparent; background: linear-gradient(105deg, #245ed5 10%, #198f9c 53%, #7357d9 94%); background-clip: text; -webkit-background-clip: text; }
        .hero-copy > p { max-width: 720px; margin: 24px 0 26px; color: var(--muted); font-size: clamp(1rem, 1.8vw, 1.17rem); }
        .capabilities { display: flex; flex-wrap: wrap; gap: 9px; }
        .capability { display: inline-flex; align-items: center; gap: 7px; padding: 8px 11px; border: 1px solid rgba(83, 119, 177, .16); border-radius: 999px; color: #425875; background: rgba(255,255,255,.56); font-size: .76rem; font-weight: 700; backdrop-filter: blur(8px); }
        .capability svg { width: 14px; height: 14px; color: var(--blue); }
        .hero-visual { position: relative; min-height: 300px; display: grid; place-items: center; }
        .molecule-orbit { position: absolute; width: 290px; aspect-ratio: 1; border: 1px solid rgba(57, 111, 202, .13); border-radius: 50%; animation: orbit 26s linear infinite; }
        .molecule-orbit::before, .molecule-orbit::after { content: ""; position: absolute; width: 10px; height: 10px; border-radius: 50%; background: var(--blue); box-shadow: 0 0 20px rgba(39,100,235,.55); }
        .molecule-orbit::before { left: 18px; top: 62px; } .molecule-orbit::after { right: 24px; bottom: 48px; background: var(--teal); }
        @keyframes orbit { to { transform: rotate(360deg); } }
        .helix-panel { width: 230px; height: 250px; padding: 22px; border: 1px solid rgba(91, 127, 184, .18); border-radius: 34px; background: linear-gradient(150deg, rgba(255,255,255,.78), rgba(239,247,255,.5)); box-shadow: var(--shadow); backdrop-filter: blur(18px); transform: rotate(4deg); }
        .helix-panel svg { width: 100%; height: 100%; filter: drop-shadow(0 12px 15px rgba(45,94,172,.13)); }
        .card { position: relative; margin-bottom: 24px; padding: 30px; border: 1px solid var(--line); border-radius: var(--radius); background: linear-gradient(145deg, var(--glass-strong), var(--glass)); box-shadow: var(--shadow); backdrop-filter: blur(20px); -webkit-backdrop-filter: blur(20px); transition: transform .25s ease, box-shadow .25s ease, border-color .25s ease; }
        .card::after { content: ""; position: absolute; inset: 0; z-index: -1; border-radius: inherit; opacity: 0; background: radial-gradient(circle at 15% 0%, rgba(61,128,246,.08), transparent 38%); transition: opacity .25s ease; }
        .card:hover { border-color: rgba(72, 123, 209, .28); box-shadow: 0 27px 70px rgba(37, 73, 132, .14); transform: translateY(-2px); }
        .card:hover::after { opacity: 1; }
        .upload-card { display: grid; grid-template-columns: minmax(0, 1fr) minmax(240px, .38fr); gap: 34px; overflow: hidden; }
        .upload-card::before { content: ""; position: absolute; width: 250px; height: 250px; right: -90px; top: -110px; border-radius: 50%; background: radial-gradient(circle, rgba(72,132,238,.13), transparent 68%); }
        .upload-main, .upload-aside { position: relative; z-index: 1; }
        .upload-aside { padding: 24px; border: 1px solid rgba(91,127,184,.14); border-radius: 17px; background: rgba(244,248,255,.62); }
        .section-kicker { margin: 0 0 6px; color: var(--teal); font-size: .71rem; font-weight: 800; letter-spacing: .1em; text-transform: uppercase; }
        h2 { display: flex; align-items: center; gap: 11px; margin: 0 0 22px; color: #15305a; font-size: 1.12rem; line-height: 1.3; letter-spacing: -.015em; }
        h2::before { content: ""; width: 9px; height: 9px; flex: 0 0 auto; border: 4px solid rgba(39,100,235,.18); border-radius: 50%; background: var(--blue); box-shadow: 0 0 0 5px rgba(39,100,235,.055); }
        h3 { color: #18355e; }
        .section-intro { margin: -12px 0 24px; color: var(--muted); font-size: .92rem; }
        .field { margin-bottom: 17px; }
        .field label { display: block; margin-bottom: 7px; color: #304765; font-size: .78rem; font-weight: 750; }
        input[type=text], input[type=file] { width: 100%; min-height: 50px; padding: 12px 14px; color: var(--ink); border: 1px solid rgba(92,123,173,.24); border-radius: 13px; background: rgba(255,255,255,.78); box-shadow: inset 0 1px 2px rgba(20,45,80,.025); transition: border-color .2s, box-shadow .2s, background .2s; }
        input[type=file] { padding: 8px; cursor: pointer; color: var(--muted); }
        input[type=file]::file-selector-button { height: 33px; margin-right: 12px; padding: 0 13px; border: 0; border-radius: 9px; color: #2454b6; background: #eaf1ff; font-weight: 750; cursor: pointer; }
        input:focus { outline: 0; border-color: rgba(39,100,235,.7); background: white; box-shadow: 0 0 0 4px rgba(39,100,235,.1), 0 8px 24px rgba(39,100,235,.07); }
        .action-row { display: flex; flex-wrap: wrap; gap: 11px; margin-top: 20px; }
        button { position: relative; min-height: 45px; display: inline-flex; align-items: center; justify-content: center; gap: 8px; padding: 10px 18px; overflow: hidden; border: 0; border-radius: 12px; color: white; background: linear-gradient(120deg, #2865e8, #4c87f1); box-shadow: 0 10px 24px rgba(39,100,235,.2); font-size: .84rem; font-weight: 760; cursor: pointer; transition: transform .18s ease, box-shadow .18s ease, filter .18s ease; }
        button::after { content: ""; position: absolute; inset: -80% -30%; background: linear-gradient(110deg, transparent 42%, rgba(255,255,255,.26), transparent 58%); transform: translateX(-80%); transition: transform .55s ease; }
        button:hover { filter: saturate(1.08); box-shadow: 0 14px 30px rgba(39,100,235,.28); transform: translateY(-2px); }
        button:hover::after { transform: translateX(80%); }
        button:active { transform: translateY(0) scale(.98); }
        button.secondary { color: #08786f; border: 1px solid rgba(13,148,136,.2); background: rgba(236,253,250,.86); box-shadow: 0 8px 20px rgba(13,148,136,.1); }
        button.secondary:hover { box-shadow: 0 12px 25px rgba(13,148,136,.16); }
        button:disabled, button.secondary:disabled { color: white; background: #aeb9c9; box-shadow: none; cursor: not-allowed; transform: none; }
        button svg { width: 17px; height: 17px; }
        .example-help { margin: 10px 0 0; color: var(--muted); font-size: .76rem; }
        code { padding: 2px 6px; border: 1px solid rgba(84,118,172,.12); border-radius: 7px; color: #385477; background: rgba(232,240,252,.75); font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .86em; }
        .aside-list { display: grid; gap: 14px; margin: 18px 0 0; padding: 0; list-style: none; }
        .aside-list li { display: flex; gap: 10px; color: #52677f; font-size: .79rem; }
        .aside-list svg { width: 17px; height: 17px; flex: 0 0 auto; color: var(--teal); }
        .note { margin-top: 19px; padding: 13px 15px; border: 1px solid rgba(217,154,38,.17); border-radius: 12px; color: #745826; background: rgba(255,249,231,.68); font-size: .79rem; }
        #progressCard { overflow: hidden; }
        .progress-head { display: flex; align-items: center; justify-content: space-between; gap: 18px; }
        .analysis-pulse { width: 38px; height: 38px; border: 1px solid rgba(39,100,235,.16); border-radius: 50%; background: radial-gradient(circle, rgba(39,100,235,.18) 0 22%, transparent 24%); animation: analyzePulse 1.8s ease-in-out infinite; }
        @keyframes analyzePulse { 50% { box-shadow: 0 0 0 12px rgba(39,100,235,0); transform: scale(1.06); } }
        .bar { height: 11px; margin-top: 18px; overflow: hidden; border: 1px solid rgba(71,105,159,.1); border-radius: 999px; background: rgba(218,228,242,.68); }
        .bar-fill { height: 100%; min-width: 0; border-radius: inherit; color: transparent; background: linear-gradient(90deg, var(--blue), #18a8a0, var(--violet)); background-size: 200% 100%; box-shadow: 0 0 18px rgba(39,100,235,.24); transition: width .45s cubic-bezier(.22,.8,.25,1); animation: progressFlow 2s linear infinite; }
        @keyframes progressFlow { to { background-position: -200% 0; } }
        #statusMsg { margin: 13px 0 0; color: var(--muted); font-size: .84rem; }
        #resultsArea { display: grid; gap: 0; }
        #resultsArea > .card { animation: enter .55s both; }
        #resultsArea > .card:nth-child(2) { animation-delay: .06s; } #resultsArea > .card:nth-child(3) { animation-delay: .12s; } #resultsArea > .card:nth-child(4) { animation-delay: .18s; } #resultsArea > .card:nth-child(5) { animation-delay: .24s; } #resultsArea > .card:nth-child(6) { animation-delay: .3s; }
        @keyframes enter { from { opacity: 0; transform: translateY(18px) scale(.99); } to { opacity: 1; transform: none; } }
        .metrics { display: grid; grid-template-columns: repeat(6, minmax(115px, 1fr)); gap: 11px; }
        .metric { position: relative; min-height: 112px; display: flex; flex-direction: column; justify-content: center; padding: 17px 14px; overflow: hidden; border: 1px solid rgba(83,119,177,.13); border-radius: 15px; background: linear-gradient(150deg, rgba(244,249,255,.92), rgba(255,255,255,.7)); text-align: left; transition: transform .2s ease, border-color .2s ease; }
        .metric::after { content: ""; position: absolute; width: 54px; height: 54px; right: -18px; top: -20px; border-radius: 50%; background: linear-gradient(145deg, rgba(39,100,235,.12), rgba(13,148,136,.08)); }
        .metric:hover { transform: translateY(-3px); border-color: rgba(39,100,235,.27); }
        .metric .val { color: #245bc9; font-size: 1.75rem; font-weight: 800; line-height: 1; letter-spacing: -.04em; }
        .metric .lbl { margin-top: 9px; color: var(--muted); font-size: .69rem; font-weight: 750; letter-spacing: .045em; text-transform: uppercase; }
        table { width: 100%; min-width: 680px; display: block; overflow-x: auto; border-spacing: 0; border-collapse: separate; border: 1px solid rgba(84,117,169,.14); border-radius: 14px; color: #314965; font-size: .79rem; }
        tbody { width: 100%; }
        th, td { min-width: 110px; padding: 13px 14px; border: 0; border-bottom: 1px solid rgba(84,117,169,.1); text-align: left; vertical-align: top; }
        th { position: sticky; top: 0; z-index: 1; color: #365271; background: rgba(235,242,252,.96); font-size: .68rem; font-weight: 800; letter-spacing: .055em; text-transform: uppercase; backdrop-filter: blur(10px); }
        tr:last-child td { border-bottom: 0; }
        tbody tr { background: rgba(255,255,255,.42); transition: background .18s ease; }
        tbody tr:nth-child(even) { background: rgba(243,247,253,.62); }
        tbody tr:hover { background: rgba(228,239,255,.72); }
        .badge { display: inline-flex; align-items: center; padding: 4px 8px; border-radius: 999px; color: white; box-shadow: 0 4px 12px rgba(40,65,100,.11); font-size: .64rem; font-weight: 800; line-height: 1; letter-spacing: .02em; }
        .section-description { margin: -11px 0 22px; color: var(--muted) !important; font-size: .82rem !important; }
        .finding-card { margin-bottom: 12px; padding: 17px 18px; border: 1px solid rgba(88,120,171,.14); border-radius: 14px; background: rgba(249,251,255,.7); transition: transform .2s, border-color .2s, box-shadow .2s; }
        .finding-card:hover { transform: translateX(3px); border-color: rgba(39,100,235,.24); box-shadow: 0 10px 24px rgba(46,78,128,.07); }
        .finding-card h3 { margin: 0 0 10px !important; display: flex; align-items: center; flex-wrap: wrap; gap: 8px; }
        .finding-meta { margin: 6px 0 !important; color: #465d78; font-size: .8rem !important; }
        .evidence-links { margin: 12px 0 0 !important; padding-top: 10px; border-top: 1px solid rgba(87,120,170,.1); font-size: .75rem !important; }
        .gene-title { display: flex; align-items: center; gap: 10px; margin: 24px 0 10px; }
        .gene-title::before { content: "GENE"; padding: 4px 7px; border-radius: 6px; color: #2764c9; background: #eaf1ff; font-size: .58rem; letter-spacing: .08em; }
        .structure-gene-button { min-height: 0; padding: 2px 0; overflow: visible; border-radius: 5px; color: #245bc9; background: none; box-shadow: none; font-size: inherit; font-weight: 800; text-decoration: underline; text-decoration-color: rgba(36,91,201,.28); text-underline-offset: 3px; }
        .structure-gene-button::after { display: none; }
        .structure-gene-button:hover { color: #123f9c; background: none; box-shadow: none; filter: none; transform: none; text-decoration-color: currentColor; }
        body.viewer-open { overflow: hidden; }
        .structure-modal[hidden] { display: none; }
        .structure-modal { position: fixed; inset: 0; z-index: 80; display: grid; place-items: center; padding: 22px; }
        .structure-modal-backdrop { position: absolute; inset: 0; background: rgba(11,25,48,.58); backdrop-filter: blur(12px); -webkit-backdrop-filter: blur(12px); animation: modalFade .2s ease both; }
        .structure-dialog { position: relative; width: min(1040px, 100%); max-height: calc(100vh - 44px); display: flex; flex-direction: column; overflow: hidden; border: 1px solid rgba(255,255,255,.52); border-radius: 24px; background: rgba(249,252,255,.97); box-shadow: 0 35px 100px rgba(7,24,53,.3); animation: modalEnter .26s cubic-bezier(.2,.8,.2,1) both; }
        @keyframes modalFade { from { opacity: 0; } }
        @keyframes modalEnter { from { opacity: 0; transform: translateY(16px) scale(.98); } }
        .structure-dialog-header { display: flex; align-items: flex-start; justify-content: space-between; gap: 20px; padding: 20px 22px 15px; border-bottom: 1px solid rgba(86,119,169,.15); }
        .structure-dialog-header h2 { margin-bottom: 5px; }
        .structure-subtitle { margin: 0; color: var(--muted); font-size: .78rem; }
        .structure-close { width: 38px; min-height: 38px; flex: 0 0 auto; padding: 0; border: 1px solid rgba(83,116,166,.17); border-radius: 11px; color: #3c526e; background: rgba(238,244,252,.86); box-shadow: none; font-size: 1.25rem; }
        .structure-close::after, .viewer-control::after { display: none; }
        .structure-close:hover { color: #17345e; background: #e4edfa; box-shadow: none; transform: none; }
        .structure-meta { display: flex; flex-wrap: wrap; gap: 8px; padding: 12px 22px; border-bottom: 1px solid rgba(86,119,169,.12); }
        .structure-meta-item { padding: 6px 9px; border: 1px solid rgba(78,115,176,.14); border-radius: 8px; color: #465e7c; background: rgba(238,245,255,.75); font-size: .69rem; font-weight: 700; }
        .viewer-shell { position: relative; min-height: 360px; height: min(59vh, 610px); overflow: hidden; background: radial-gradient(circle at 50% 45%, #ffffff, #edf3fb); }
        #nglViewport { width: 100%; height: 100%; }
        .viewer-loading { position: absolute; inset: 0; z-index: 3; display: grid; place-items: center; padding: 24px; color: #536984; background: rgba(244,248,253,.88); font-size: .82rem; text-align: center; backdrop-filter: blur(6px); }
        .viewer-loading[hidden] { display: none; }
        .viewer-loading::before { content: ""; width: 30px; height: 30px; margin-bottom: 48px; position: absolute; border: 3px solid #d5e2f5; border-top-color: var(--blue); border-radius: 50%; animation: viewerSpin .8s linear infinite; }
        .viewer-loading.is-error { color: #9f2838; background: rgba(255,245,246,.94); }
        .viewer-loading.is-error::before { display: none; }
        @keyframes viewerSpin { to { transform: rotate(360deg); } }
        .variant-label { position: absolute; z-index: 2; top: 14px; left: 14px; max-width: calc(100% - 28px); padding: 8px 11px; border: 1px solid rgba(239,68,68,.2); border-radius: 9px; color: #a41f2d; background: rgba(255,247,247,.9); box-shadow: 0 8px 24px rgba(97,31,42,.1); font-size: .72rem; font-weight: 800; pointer-events: none; backdrop-filter: blur(8px); }
        .plddt-legend { position: absolute; z-index: 2; right: 14px; bottom: 14px; display: flex; flex-wrap: wrap; justify-content: flex-end; gap: 7px 10px; max-width: calc(100% - 28px); padding: 8px 10px; border: 1px solid rgba(78,109,157,.14); border-radius: 9px; color: #536984; background: rgba(255,255,255,.86); box-shadow: 0 8px 24px rgba(35,62,103,.08); font-size: .62rem; backdrop-filter: blur(8px); pointer-events: none; }
        .plddt-legend span { display: inline-flex; align-items: center; gap: 4px; white-space: nowrap; }
        .plddt-dot { width: 7px; height: 7px; border-radius: 50%; }
        .structure-dialog-footer { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 10px; padding: 13px 22px; border-top: 1px solid rgba(86,119,169,.14); }
        .viewer-actions { display: flex; flex-wrap: wrap; gap: 8px; }
        .viewer-control { min-height: 35px; padding: 7px 11px; border: 1px solid rgba(69,108,171,.18); border-radius: 9px; color: #315681; background: rgba(235,243,254,.9); box-shadow: none; font-size: .7rem; }
        .viewer-control:hover { color: #173d72; background: #e1ecfc; box-shadow: none; transform: none; }
        .alphafold-secondary-link { font-size: .73rem; }
        .mutation-chart-host { min-height: 320px; }
        .mutation-chart-layout { display: grid; grid-template-columns: minmax(220px, 1fr) minmax(280px, 1.25fr); gap: 34px; align-items: center; }
        .mutation-legend { display: flex; flex-direction: column; gap: 7px; }
        .mutation-legend-item { width: 100%; min-height: 44px; display: grid; grid-template-columns: 14px minmax(0, 1fr) auto; gap: 10px; align-items: center; padding: 9px 11px; border: 1px solid transparent; border-radius: 11px; color: var(--ink); background: transparent; box-shadow: none; text-align: left; }
        .mutation-legend-item::after { display: none; }
        .mutation-legend-item:hover, .mutation-legend-item.is-active { color: var(--ink); border-color: #c8d9f5; background: #eff5ff; box-shadow: none; transform: translateX(3px); }
        .mutation-swatch { width: 11px; height: 11px; border-radius: 3px; box-shadow: 0 3px 8px rgba(40,65,100,.16); }
        .mutation-name { overflow-wrap: anywhere; font-size: .77rem; font-weight: 650; }
        .mutation-value { color: var(--muted); font-size: .72rem; font-variant-numeric: tabular-nums; white-space: nowrap; }
        .donut-figure { position: relative; display: grid; place-items: center; min-width: 0; }
        .donut-chart { width: min(100%, 340px); height: auto; overflow: visible; filter: drop-shadow(0 16px 22px rgba(41,78,137,.1)); }
        .donut-segment { cursor: pointer; stroke-width: 34; transition: stroke-width .15s ease, filter .15s ease, opacity .15s ease; }
        .donut-segment.is-active, .donut-segment:focus { stroke-width: 42; filter: drop-shadow(0 2px 3px rgba(31,41,55,.2)); }
        .donut-segment:focus { outline: none; }
        .donut-center-total { fill: #17366b; font-size: 30px; font-weight: 800; text-anchor: middle; }
        .donut-center-label { fill: var(--muted); font-size: 10px; text-anchor: middle; }
        .chart-tooltip { position: absolute; z-index: 5; max-width: 230px; padding: 9px 11px; border: 1px solid rgba(255,255,255,.12); border-radius: 10px; color: white; background: rgba(16,33,59,.94); box-shadow: 0 12px 30px rgba(15,31,55,.22); font-size: .72rem; line-height: 1.45; pointer-events: none; opacity: 0; transform: translateY(3px); transition: opacity .12s, transform .12s; backdrop-filter: blur(10px); }
        .chart-tooltip.visible { opacity: 1; transform: translateY(0); }
        .chart-tooltip strong { display: block; overflow-wrap: anywhere; }
        .chart-caption { margin: 12px 0 0; color: var(--muted); font-size: .72rem; text-align: center; }
        .footer { padding: 4px 0 35px; color: #73839a; text-align: center; font-size: .72rem; }
        .footer strong { color: #435a77; }
        @media (max-width: 980px) { .hero { grid-template-columns: 1fr; gap: 20px; } .hero-visual { display: none; } .upload-card { grid-template-columns: 1fr; } .metrics { grid-template-columns: repeat(3, 1fr); } }
        @media (max-width: 680px) { .topbar-inner, .container { width: min(100% - 28px, 1180px); } .topbar-inner { min-height: 64px; } .platform-status { font-size: 0; padding: 9px; } .container { padding-top: 45px; } h1 { font-size: clamp(2.15rem, 12vw, 3.2rem); } .hero { margin-bottom: 30px; } .card { padding: 21px; border-radius: 18px; } .upload-aside { padding: 18px; } .metrics { grid-template-columns: repeat(2, 1fr); } .mutation-chart-layout { grid-template-columns: 1fr; } .mutation-legend { order: 2; } .donut-figure { order: 1; } .action-row button { width: 100%; } }
        @media (max-width: 410px) { .brand-copy span { display: none; } .metrics { grid-template-columns: 1fr 1fr; } .metric { min-width: 0; padding: 14px 10px; } }
        @media (prefers-reduced-motion: reduce) { *, *::before, *::after { scroll-behavior: auto !important; animation-duration: .01ms !important; animation-iteration-count: 1 !important; transition-duration: .01ms !important; } .science-bg svg { transform: none; } }
        @media print { .science-bg, .topbar, .hero, .upload-card, #progressCard, .footer { display: none !important; } body { background: white; } .container { width: 100%; padding: 0; } .card { box-shadow: none; break-inside: avoid; } }
    </style>
</head>
<body>
<a class="skip-link" href="#mainContent">Skip to analysis workspace</a>
<div class="science-bg" aria-hidden="true">
    <svg viewBox="0 0 1440 1000" preserveAspectRatio="xMidYMid slice">
        <g class="float-a" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M92 115c80 55 80 140 0 195s-80 140 0 195"/><path d="M188 115c-80 55-80 140 0 195s80 140 0 195"/>
            <path d="M112 145h56M94 202h94M94 260h94M110 318h60M94 376h94M112 434h56"/>
        </g>
        <g class="float-b" transform="translate(1120 120)" fill="none" stroke="currentColor" stroke-width="2">
            <polygon points="70,0 140,40 140,120 70,160 0,120 0,40"/><polygon points="188,68 235,95 235,149 188,176 141,149 141,95"/>
            <line x1="140" y1="82" x2="159" y2="93"/><circle cx="70" cy="0" r="7" fill="currentColor"/><circle cx="235" cy="149" r="7" fill="currentColor"/>
        </g>
        <g class="float-c" transform="translate(1045 690)" fill="none" stroke="currentColor" stroke-width="1.5">
            <circle cx="100" cy="100" r="82"/><circle cx="100" cy="100" r="45"/><ellipse cx="100" cy="100" rx="118" ry="36" transform="rotate(35 100 100)"/><ellipse cx="100" cy="100" rx="118" ry="36" transform="rotate(-35 100 100)"/>
            <circle cx="37" cy="56" r="7" fill="currentColor"/><circle cx="168" cy="45" r="7" fill="currentColor"/><circle cx="174" cy="159" r="7" fill="currentColor"/>
        </g>
        <g fill="currentColor">
            <circle class="particle" cx="360" cy="170" r="3"/><circle class="particle" cx="480" cy="90" r="2"/><circle class="particle" cx="830" cy="180" r="3"/><circle class="particle" cx="930" cy="430" r="2"/><circle class="particle" cx="270" cy="740" r="3"/><circle class="particle" cx="720" cy="850" r="2"/>
        </g>
        <g opacity=".55" fill="none" stroke="currentColor" stroke-width="1">
            <path d="M360 170L480 90L620 215L830 180L930 430"/><path d="M270 740L460 620L720 850L930 720"/>
            <circle cx="620" cy="215" r="5"/><circle cx="460" cy="620" r="5"/><circle cx="930" cy="720" r="5"/>
        </g>
    </svg>
</div>
<header class="topbar">
    <div class="topbar-inner">
        <a class="brand" href="/" aria-label="HelixAI home">
            <span class="brand-mark" aria-hidden="true"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8"><path d="M7 3c7 4 3 14 10 18M17 3C10 7 14 17 7 21M8.5 6h7M7 12h10M8.5 18h7"/></svg></span>
            <span class="brand-copy"><strong>HelixAI</strong><span>Precision Medicine Intelligence</span></span>
        </a>
        <span class="platform-status"><span class="status-dot"></span> Research pipeline online</span>
    </div>
</header>
<main class="container" id="mainContent">
    <section class="hero" aria-labelledby="pageTitle">
        <div class="hero-copy">
            <div class="eyebrow">AI-assisted pharmacogenomics</div>
            <h1 id="pageTitle">Decode variants.<br><span class="gradient-text">Discover response.</span></h1>
            <p>Transform genomic variants into interpretable molecular, structural, and pharmacological evidence with an integrated precision medicine workflow.</p>
            <div class="capabilities" aria-label="Platform capabilities">
                <span class="capability"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M7 3c7 4 3 14 10 18M17 3C10 7 14 17 7 21"/></svg> Ensembl VEP</span>
                <span class="capability"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="3"/><path d="M4.9 4.9l4.9 4.9M14.2 14.2l4.9 4.9M19.1 4.9l-4.9 4.9M9.8 14.2l-4.9 4.9"/></svg> AlphaFold</span>
                <span class="capability"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M9 3h6l4 7-7 11-7-11z"/><path d="M5 10h14"/></svg> Drug Evidence</span>
                <span class="capability"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="2"/><circle cx="5" cy="6" r="2"/><circle cx="19" cy="6" r="2"/><path d="M6.5 7.5l4 3M17.5 7.5l-4 3M12 14v5"/></svg> Molecular Intelligence</span>
            </div>
        </div>
        <div class="hero-visual" aria-hidden="true">
            <div class="molecule-orbit"></div>
            <div class="helix-panel"><svg viewBox="0 0 180 210" fill="none" stroke-linecap="round"><defs><linearGradient id="helix" x1="0" y1="0" x2="180" y2="210"><stop stop-color="#2d68df"/><stop offset=".52" stop-color="#14a19a"/><stop offset="1" stop-color="#765bd9"/></linearGradient></defs><path d="M47 8c93 45-5 147 87 194M133 8C40 53 138 155 46 202" stroke="url(#helix)" stroke-width="7"/><g stroke="#7e9cc7" stroke-width="2"><path d="M61 23h58M43 54h94M42 88h96M44 122h92M43 157h94M60 189h60"/></g><g fill="#fff" stroke="url(#helix)" stroke-width="3"><circle cx="61" cy="23" r="5"/><circle cx="119" cy="23" r="5"/><circle cx="43" cy="54" r="5"/><circle cx="137" cy="54" r="5"/><circle cx="42" cy="88" r="5"/><circle cx="138" cy="88" r="5"/><circle cx="44" cy="122" r="5"/><circle cx="136" cy="122" r="5"/><circle cx="43" cy="157" r="5"/><circle cx="137" cy="157" r="5"/><circle cx="60" cy="189" r="5"/><circle cx="120" cy="189" r="5"/></g></svg></div>
        </div>
    </section>

    <section class="card upload-card" aria-labelledby="uploadHeading">
        <div class="upload-main">
            <p class="section-kicker">Analysis workspace</p>
            <h2 id="uploadHeading">Upload genomic variants</h2>
            <p class="section-intro">Submit a VCF dataset or launch the validated demonstration workflow.</p>
            <form id="uploadForm">
                <div class="field"><label for="sampleId">Sample identifier</label><input type="text" id="sampleId" placeholder="e.g. SAMPLE_001" autocomplete="off" required></div>
                <div class="field"><label for="vcfFile">Variant Call Format file</label><input type="file" id="vcfFile" accept=".vcf,.vcf.gz" required></div>
                <div class="action-row">
                    <button type="submit" id="submitBtn"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 16V4M7 9l5-5 5 5M5 20h14"/></svg> Upload and Analyze</button>
                    <button type="button" id="exampleBtn" class="secondary"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M8 5v14l11-7z"/></svg> Run Example Analysis</button>
                </div>
                <p class="example-help">Example source: <code>sample_pharmacogenomics.vcf</code></p>
            </form>
        </div>
        <aside class="upload-aside" aria-label="Analysis capabilities">
            <p class="section-kicker">Integrated workflow</p>
            <strong>From sequence to evidence</strong>
            <ul class="aside-list">
                <li><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 6L9 17l-5-5"/></svg><span>Variant annotation and coding consequence prioritization</span></li>
                <li><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 6L9 17l-5-5"/></svg><span>Protein mapping with structural confidence analysis</span></li>
                <li><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M20 6L9 17l-5-5"/></svg><span>PharmGKB evidence and diabetes-focused interpretation</span></li>
            </ul>
            <div class="note">Large datasets are annotated to the configured VEP cap to maintain reliable API throughput.</div>
        </aside>
    </section>

    <section class="card" id="progressCard" style="display:none;" aria-live="polite">
        <div class="progress-head"><div><p class="section-kicker">Computational pipeline</p><h2>Analysis in progress</h2></div><div class="analysis-pulse" aria-hidden="true"></div></div>
        <div class="bar" role="progressbar" aria-label="Analysis progress"><div class="bar-fill" id="barFill" style="width:0%">0%</div></div>
        <p id="statusMsg">Starting...</p>
    </section>

    <div id="resultsArea" aria-live="polite"></div>
</main>
<footer class="footer"><strong>HelixAI Pharmacogenomics Intelligence</strong> &middot; Research and educational use only &middot; Confirm findings with primary evidence sources</footer>

<div id="structureModal" class="structure-modal" role="dialog" aria-modal="true" aria-labelledby="structureViewerTitle" hidden>
    <div class="structure-modal-backdrop" data-close-structure-viewer></div>
    <section class="structure-dialog" tabindex="-1">
        <header class="structure-dialog-header">
            <div><p class="section-kicker">Interactive molecular model</p><h2 id="structureViewerTitle">Protein structure</h2><p id="structureViewerSubtitle" class="structure-subtitle"></p></div>
            <button type="button" id="structureViewerClose" class="structure-close" aria-label="Close structure viewer">&times;</button>
        </header>
        <div id="structureViewerMeta" class="structure-meta"></div>
        <div class="viewer-shell" id="viewerShell" aria-busy="true">
            <div id="nglViewport" aria-label="Interactive three-dimensional protein structure"></div>
            <div id="structureViewerLoading" class="viewer-loading">Loading the AlphaFold structure...</div>
            <div id="structureVariantLabel" class="variant-label"></div>
            <div class="plddt-legend" aria-label="AlphaFold confidence colors"><span><i class="plddt-dot" style="background:#0053D6"></i>Very high</span><span><i class="plddt-dot" style="background:#65CBF3"></i>Confident</span><span><i class="plddt-dot" style="background:#FFDB13"></i>Low</span><span><i class="plddt-dot" style="background:#FF7D45"></i>Very low</span></div>
        </div>
        <footer class="structure-dialog-footer">
            <div class="viewer-actions"><button type="button" id="focusResidueBtn" class="viewer-control">Focus mutated residue</button><button type="button" id="fullStructureBtn" class="viewer-control">View full protein</button></div>
            <a id="alphafoldDbLink" class="alphafold-secondary-link" href="#" target="_blank" rel="noopener noreferrer">Open in AlphaFold DB &#8599;</a>
        </footer>
    </section>
</div>

<script>
if (!window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
    window.addEventListener('pointermove', function(event) {
        const x = (event.clientX / window.innerWidth - 0.5) * 8;
        const y = (event.clientY / window.innerHeight - 0.5) * 8;
        document.documentElement.style.setProperty('--parallax-x', x.toFixed(2) + 'px');
        document.documentElement.style.setProperty('--parallax-y', y.toFixed(2) + 'px');
    }, { passive: true });
}
let pollTimer = null;
let currentAnalysisId = null;

function setAnalysisControlsDisabled(disabled) {
    document.getElementById('submitBtn').disabled = disabled;
    document.getElementById('exampleBtn').disabled = disabled;
}

function beginAnalysis(uploadId, message) {
    currentAnalysisId = uploadId;
    document.getElementById('progressCard').style.display = 'block';
    document.getElementById('resultsArea').innerHTML = '';
    setProgress(5, message);
    pollStatus(uploadId);
}

document.getElementById('uploadForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const file = document.getElementById('vcfFile').files[0];
    const sampleId = document.getElementById('sampleId').value;
    if (!file || !sampleId) return;

    const fd = new FormData();
    fd.append('file', file);
    fd.append('sample_id', sampleId);

    setAnalysisControlsDisabled(true);
    document.getElementById('progressCard').style.display = 'block';
    document.getElementById('resultsArea').innerHTML = '';
    setProgress(5, 'Uploading file...');

    try {
        const res = await fetch('/api/upload', { method: 'POST', body: fd });
        const data = await res.json();
        if (!res.ok) throw new Error(data.detail || 'Upload failed');
        beginAnalysis(data.upload_id, 'Upload complete. Starting analysis...');
    } catch (err) {
        setProgress(0, 'Error: ' + err.message);
        setAnalysisControlsDisabled(false);
    }
});

document.getElementById('exampleBtn').addEventListener('click', async () => {
    setAnalysisControlsDisabled(true);
    document.getElementById('progressCard').style.display = 'block';
    document.getElementById('resultsArea').innerHTML = '';
    setProgress(5, 'Loading example VCF...');

    try {
        const res = await fetch('/api/example', { method: 'POST' });
        const data = await res.json();
        if (!res.ok) throw new Error(data.detail || 'Could not start example analysis');
        beginAnalysis(data.upload_id, 'Example loaded. Starting analysis...');
    } catch (err) {
        setProgress(0, 'Error: ' + err.message);
        setAnalysisControlsDisabled(false);
    }
});

function setProgress(pct, msg) {
    const fill = document.getElementById('barFill');
    const progressBar = fill.parentElement;
    fill.style.width = pct + '%';
    fill.textContent = Math.round(pct) + '%';
    progressBar.setAttribute('aria-valuemin', '0');
    progressBar.setAttribute('aria-valuemax', '100');
    progressBar.setAttribute('aria-valuenow', String(Math.round(pct)));
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
                setAnalysisControlsDisabled(false);
                renderResults(data.results);
            } else if (data.status === 'failed') {
                clearInterval(pollTimer);
                setAnalysisControlsDisabled(false);
                setProgress(0, data.message);
            }
        } catch (err) {
            clearInterval(pollTimer);
            setAnalysisControlsDisabled(false);
            setProgress(0, 'Error checking analysis status: ' + err.message);
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

    // Interactive mutation-type donut chart (rendered after the HTML is inserted).
    html += '<div class="card"><h2>Mutation Types</h2>' +
            '<div id="mutationDonutChart" class="mutation-chart-host" ' +
            'aria-label="Mutation consequence distribution"></div></div>';

    // Mutations table + download
    html += '<div class="card"><h2>Coding Mutations</h2>';
    if (currentAnalysisId) {
        html += '<p><a href="/api/analysis/' + currentAnalysisId + '/mutations.csv" ' +
                'target="_blank"><button type="button">&#8681; Download mutations (CSV)</button></a></p>';
    }
    if ((r.mutations || []).length === 0) {
        html += '<p>No coding mutations found in the annotated set.</p>';
    } else {
        html += '<table><tr><th>Position</th><th>Change</th><th>Gene</th>' +
                '<th>Consequence</th><th>Amino Acid (HGVSp)</th></tr>';
        r.mutations.slice(0, 200).forEach(m => {
            var isDiab = (r.diabetes_findings || []).some(function(f){return f.gene_symbol === m.gene_symbol;});
            var geneCell = (m.gene_symbol || '-') + (isDiab ? ' <span class="badge">diabetes</span>' : '');
            html += '<tr><td>' + m.position + '</td><td>' + m.ref + '&rarr;' + m.alt + '</td>' +
                    '<td>' + geneCell + '</td>' +
                    '<td>' + (m.consequence_terms || []).join(', ') + '</td>' +
                    '<td>' + (m.hgvsp || m.amino_acid_change || '-') + '</td></tr>';
        });
        html += '</table>';
        if (r.mutations.length > 200) html += '<p>Showing first 200 of ' + r.mutations.length + '.</p>';
    }
    html += '</div>';

    // Diabetes interpretation
    html += '<div class="card"><h2>Diabetes-Associated Findings</h2>';
    html += '<p class="section-description">Interpretation from curated diabetes gene knowledge ' +
            '(OMIM, PharmGKB, GWAS Catalog). For research/education only, not diagnostic. ' +
            'Use the links to confirm each finding.</p>';
    var df = r.diabetes_findings || [];
    if (df.length === 0) {
        html += '<p>No mutations in known diabetes-associated genes were found in the annotated set.</p>';
    } else {
        df.forEach(function(f){
            html += '<div class="finding-card">';
            html += '<h3>' + f.gene_symbol +
                    ' <span class="badge">' + f.category + '</span></h3>';
            html += '<p class="finding-meta"><b>Variant:</b> ' + f.position + ' ' + f.change +
                    ' | ' + (f.consequence_terms||[]).join(', ') +
                    (f.amino_acid_change ? ' | ' + f.amino_acid_change : '') + '</p>';
            html += '<p class="finding-meta"><b>Gene role:</b> ' + f.role + '</p>';
            html += '<p class="finding-meta"><b>Diabetes significance:</b> ' + f.significance + '</p>';
            html += '<p class="finding-meta"><b>Treatment relevance:</b> ' + f.treatment + '</p>';
            html += '<p class="evidence-links">Confirm on: ' +
                    '<a href="' + f.pharmgkb_url + '" target="_blank" rel="noopener noreferrer">PharmGKB</a> &middot; ' +
                    '<a href="' + f.omim_url + '" target="_blank" rel="noopener noreferrer">OMIM</a> &middot; ' +
                    '<a href="' + f.gwas_url + '" target="_blank" rel="noopener noreferrer">GWAS Catalog</a></p>';
            html += '</div>';
        });
    }
    html += '</div>';

    // Structures
    html += '<div class="card"><h2>Protein Structures (mutated residue)</h2>';
    if ((r.structures || []).length === 0) {
        html += '<p>No structures resolved.</p>';
    } else {
        html += '<table><tr><th>Gene</th><th>UniProt</th><th>Residue</th><th>Change</th>' +
                '<th>pLDDT</th><th>In Domain</th><th>Impact</th><th>View</th></tr>';
        r.structures.forEach(function(st, structureIndex) {
            if (st.error) {
                html += '<tr><td>' + st.gene_symbol + '</td><td>' + (st.uniprot_id||'-') +
                        '</td><td colspan="6">error: ' + st.error + '</td></tr>';
            } else {
                html += '<tr><td><button type="button" class="structure-gene-button" data-structure-index="' + structureIndex + '" ' +
                        'aria-label="View ' + st.gene_symbol + ' structure at residue ' + st.residue_position + '">' + st.gene_symbol + '</button></td>' +
                        '<td>' + st.uniprot_id + '</td>' +
                        '<td>' + st.residue_position + '</td><td>' + (st.amino_acid_change||'-') + '</td>' +
                        '<td>' + (st.plddt_score!=null ? st.plddt_score.toFixed(1) : '-') + '</td>' +
                        '<td>' + (st.in_domain ? 'Yes' : 'No') + '</td>' +
                        '<td>' + (st.structural_impact||'-') + '</td>' +
                        '<td><button type="button" class="viewer-control" data-structure-index="' + structureIndex + '">3D view</button></td></tr>';
            }
        });
        html += '</table>';
    }
    html += '</div>';

    // Drugs
    html += '<div class="card"><h2>Drug Interactions</h2>';
    html += '<p class="section-description">Evidence comes from PharmGKB clinical annotations. ' +
            'Open the exact annotation to confirm each gene\\u2013drug association. ' +
            'Levels 1A/1B are strongest, 2A/2B moderate, and 3/4 lower. ' +
            '<i>no_pgx_evidence</i> indicates no PharmGKB clinical annotation was found for that pair.</p>';
    if ((r.drug_results || []).length === 0) {
        html += '<p>No drug interactions found.</p>';
    } else {
        r.drug_results.forEach(dr => {
            var gene = dr.gene_symbol;
            var genePgkb = 'https://www.pharmgkb.org/search?query=' + encodeURIComponent(gene);
            html += '<h3 class="gene-title">' + gene +
                    ' <a href="' + genePgkb + '" target="_blank" rel="noopener noreferrer">PharmGKB gene &#8599;</a></h3>';
            html += '<table><tr><th>Drug</th><th>Action</th><th>Evidence</th><th>Confirm on</th></tr>';
            (dr.matched_drugs || []).forEach(d => {
                var name = d.drug_name || '-';
                var q = encodeURIComponent(gene + ' ' + name);
                // Prefer the exact PharmGKB clinical annotation page if we have one,
                // otherwise fall back to a gene+drug search.
                var pgkb = d.pharmgkb_url || ('https://www.pharmgkb.org/search?query=' + q);
                // DrugBank drug search
                var db = 'https://go.drugbank.com/unearth/q?searcher=drugs&query=' + encodeURIComponent(name);
                var links = '<a href="' + pgkb + '" target="_blank">PharmGKB</a> | ' +
                            '<a href="' + db + '" target="_blank">DrugBank</a>';
                // Colour-code the evidence level
                var ev = d.evidence_level || 'no_pgx_evidence';
                var evColor = ev.indexOf('1') > -1 ? '#059669' :
                              ev.indexOf('2') > -1 ? '#2563eb' :
                              ev.indexOf('3') > -1 || ev.indexOf('4') > -1 ? '#d97706' : '#6b7280';
                var evCell = '<span style="color:' + evColor + ';font-weight:600;">' + ev + '</span>';
                html += '<tr><td>' + name + '</td><td>' + (d.action||'-') +
                        '</td><td>' + evCell + '</td><td>' + links + '</td></tr>';
            });
            html += '</table>';
        });
    }
    html += '</div>';

    document.getElementById('resultsArea').innerHTML = html;
    document.getElementById('resultsArea').querySelectorAll('table').forEach(function(table) {
        table.setAttribute('tabindex', '0');
        table.setAttribute('aria-label', 'Scrollable scientific results table');
    });
    document.getElementById('resultsArea').querySelectorAll('[data-structure-index]').forEach(function(control) {
        control.addEventListener('click', function() {
            openStructureViewer((r.structures || [])[Number(control.dataset.structureIndex)], control);
        });
    });
    renderMutationDonut(r.mutation_type_counts || {});
}

function renderMutationDonut(rawCounts) {
    const host = document.getElementById('mutationDonutChart');
    if (!host) return;

    const entries = Object.keys(rawCounts || {}).map(function(name) {
        return { name: name, count: Number(rawCounts[name]) };
    }).filter(function(item) {
        return item.name && Number.isFinite(item.count) && item.count > 0;
    }).sort(function(a, b) {
        return b.count - a.count || a.name.localeCompare(b.name);
    });

    host.replaceChildren();
    if (entries.length === 0) {
        const empty = document.createElement('p');
        empty.textContent = 'No coding mutation types to chart.';
        host.appendChild(empty);
        host.style.minHeight = '0';
        return;
    }

    host.style.minHeight = '';
    const total = entries.reduce(function(sum, item) { return sum + item.count; }, 0);
    const colors = ['#2563eb', '#7c3aed', '#0891b2', '#059669', '#d97706', '#dc2626', '#db2777', '#4f46e5', '#65a30d', '#ea580c'];
    const svgNS = 'http://www.w3.org/2000/svg';
    const layout = document.createElement('div');
    layout.className = 'mutation-chart-layout';
    const legend = document.createElement('div');
    legend.className = 'mutation-legend';
    legend.setAttribute('aria-label', 'Mutation type legend');
    const figure = document.createElement('div');
    figure.className = 'donut-figure';
    const svg = document.createElementNS(svgNS, 'svg');
    svg.classList.add('donut-chart');
    svg.setAttribute('viewBox', '0 0 260 260');
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', 'Donut chart of mutation consequence calls');
    const title = document.createElementNS(svgNS, 'title');
    title.textContent = 'Mutation consequence distribution: ' + total + ' total calls';
    svg.appendChild(title);

    const background = document.createElementNS(svgNS, 'circle');
    background.setAttribute('cx', '130');
    background.setAttribute('cy', '130');
    background.setAttribute('r', '86');
    background.setAttribute('fill', 'none');
    background.setAttribute('stroke', '#e5e7eb');
    background.setAttribute('stroke-width', '34');
    svg.appendChild(background);

    const circumference = 2 * Math.PI * 86;
    let offset = 0;
    const segments = [];
    const legendItems = [];
    const tooltip = document.createElement('div');
    tooltip.className = 'chart-tooltip';
    tooltip.setAttribute('role', 'tooltip');

    function setActive(index, active) {
        segments[index].classList.toggle('is-active', active);
        legendItems[index].classList.toggle('is-active', active);
        segments.forEach(function(segment, i) {
            segment.style.opacity = active && i !== index ? '0.45' : '1';
        });
    }

    function showTooltip(index, clientX, clientY) {
        const item = entries[index];
        const pct = item.count / total * 100;
        tooltip.replaceChildren();
        const strong = document.createElement('strong');
        strong.textContent = item.name;
        const detail = document.createElement('span');
        detail.textContent = item.count + ' call' + (item.count === 1 ? '' : 's') + ' \u2022 ' + pct.toFixed(1) + '%';
        tooltip.append(strong, detail);
        tooltip.classList.add('visible');
        const rect = figure.getBoundingClientRect();
        let left = clientX == null ? rect.width / 2 : clientX - rect.left;
        let top = clientY == null ? rect.height / 2 : clientY - rect.top;
        left = Math.max(8, Math.min(left + 12, rect.width - tooltip.offsetWidth - 8));
        top = Math.max(8, Math.min(top + 12, rect.height - tooltip.offsetHeight - 8));
        tooltip.style.left = left + 'px';
        tooltip.style.top = top + 'px';
    }

    function hideTooltip(index) {
        setActive(index, false);
        tooltip.classList.remove('visible');
    }

    entries.forEach(function(item, index) {
        const fraction = item.count / total;
        const segment = document.createElementNS(svgNS, 'circle');
        segment.classList.add('donut-segment');
        segment.setAttribute('cx', '130');
        segment.setAttribute('cy', '130');
        segment.setAttribute('r', '86');
        segment.setAttribute('fill', 'none');
        segment.setAttribute('stroke', colors[index % colors.length]);
        segment.setAttribute('stroke-dasharray', (fraction * circumference) + ' ' + circumference);
        segment.setAttribute('stroke-dashoffset', String(-offset * circumference));
        segment.setAttribute('transform', 'rotate(-90 130 130)');
        segment.setAttribute('tabindex', '0');
        segment.setAttribute('role', 'img');
        segment.setAttribute('aria-label', item.name + ': ' + item.count + ', ' + (fraction * 100).toFixed(1) + ' percent');
        offset += fraction;
        svg.appendChild(segment);
        segments.push(segment);

        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'mutation-legend-item';
        button.setAttribute('aria-label', item.name + ': ' + item.count + ', ' + (fraction * 100).toFixed(1) + ' percent');
        const swatch = document.createElement('span');
        swatch.className = 'mutation-swatch';
        swatch.style.backgroundColor = colors[index % colors.length];
        const name = document.createElement('span');
        name.className = 'mutation-name';
        name.textContent = item.name;
        const value = document.createElement('span');
        value.className = 'mutation-value';
        value.textContent = item.count + ' (' + (fraction * 100).toFixed(1) + '%)';
        button.append(swatch, name, value);
        legend.appendChild(button);
        legendItems.push(button);

        segment.addEventListener('pointerenter', function(event) { setActive(index, true); showTooltip(index, event.clientX, event.clientY); });
        segment.addEventListener('pointermove', function(event) { showTooltip(index, event.clientX, event.clientY); });
        segment.addEventListener('pointerleave', function() { hideTooltip(index); });
        segment.addEventListener('focus', function() { setActive(index, true); showTooltip(index); });
        segment.addEventListener('blur', function() { hideTooltip(index); });
        button.addEventListener('pointerenter', function() { setActive(index, true); showTooltip(index); });
        button.addEventListener('pointerleave', function() { hideTooltip(index); });
        button.addEventListener('focus', function() { setActive(index, true); showTooltip(index); });
        button.addEventListener('blur', function() { hideTooltip(index); });
    });

    const totalText = document.createElementNS(svgNS, 'text');
    totalText.classList.add('donut-center-total');
    totalText.setAttribute('x', '130');
    totalText.setAttribute('y', '126');
    totalText.textContent = String(total);
    svg.appendChild(totalText);
    const labelText = document.createElementNS(svgNS, 'text');
    labelText.classList.add('donut-center-label');
    labelText.setAttribute('x', '130');
    labelText.setAttribute('y', '146');
    labelText.textContent = 'consequence calls';
    svg.appendChild(labelText);

    figure.append(svg, tooltip);
    const caption = document.createElement('p');
    caption.className = 'chart-caption';
    caption.textContent = 'Hover over, tap, or focus a segment to inspect its count and percentage.';
    figure.appendChild(caption);
    layout.append(legend, figure);
    host.appendChild(layout);
}

function metric(val, label) {
    return '<div class="metric"><div class="val">' + (val!=null?val:0) +
           '</div><div class="lbl">' + label + '</div></div>';
}

let nglScriptPromise = null;
let structureStage = null;
let structureComponent = null;
let activeResidueSelection = null;
let structureViewerTrigger = null;

function loadNglViewer() {
    if (window.NGL) return Promise.resolve(window.NGL);
    if (nglScriptPromise) return nglScriptPromise;
    nglScriptPromise = new Promise(function(resolve, reject) {
        const script = document.createElement('script');
        script.src = 'https://cdnjs.cloudflare.com/ajax/libs/ngl/2.0.0-dev.39/ngl.min.js';
        script.async = true;
        script.onload = function() { window.NGL ? resolve(window.NGL) : reject(new Error('NGL did not initialize')); };
        script.onerror = function() { reject(new Error('The 3D viewer library could not be loaded')); };
        document.head.appendChild(script);
    });
    return nglScriptPromise;
}

function alphaFoldConfidenceScheme() {
    return window.NGL.ColormakerRegistry.addScheme(function() {
        this.atomColor = function(atom) {
            if (atom.bfactor >= 90) return 0x0053D6;
            if (atom.bfactor >= 70) return 0x65CBF3;
            if (atom.bfactor >= 50) return 0xFFDB13;
            return 0xFF7D45;
        };
    });
}

async function resolveAlphaFoldModel(uniprotId) {
    const fallback = 'https://alphafold.ebi.ac.uk/files/AF-' + encodeURIComponent(uniprotId) + '-F1-model_v4.pdb';
    try {
        const response = await fetch('https://alphafold.ebi.ac.uk/api/prediction/' + encodeURIComponent(uniprotId));
        if (!response.ok) return fallback;
        const records = await response.json();
        return records && records[0] && records[0].pdbUrl ? records[0].pdbUrl : fallback;
    } catch (error) {
        return fallback;
    }
}

async function openStructureViewer(structure, trigger) {
    if (!structure || structure.error || !structure.uniprot_id || !structure.residue_position) return;
    structureViewerTrigger = trigger || document.activeElement;
    const modal = document.getElementById('structureModal');
    const loading = document.getElementById('structureViewerLoading');
    const shell = document.getElementById('viewerShell');
    const residue = Number(structure.residue_position);
    activeResidueSelection = String(residue);
    document.getElementById('structureViewerTitle').textContent = (structure.gene_symbol || 'Protein') + ' structure';
    document.getElementById('structureViewerSubtitle').textContent = 'UniProt ' + structure.uniprot_id + ' \u2022 mutated residue ' + residue;
    document.getElementById('structureVariantLabel').textContent = (structure.amino_acid_change || 'Variant') + ' \u2022 residue ' + residue;
    document.getElementById('structureViewerMeta').innerHTML =
        '<span class="structure-meta-item">Gene: ' + (structure.gene_symbol || '-') + '</span>' +
        '<span class="structure-meta-item">UniProt: ' + structure.uniprot_id + '</span>' +
        '<span class="structure-meta-item">Residue: ' + residue + '</span>' +
        '<span class="structure-meta-item">Change: ' + (structure.amino_acid_change || '-') + '</span>' +
        '<span class="structure-meta-item">pLDDT: ' + (structure.plddt_score != null ? Number(structure.plddt_score).toFixed(1) : '-') + '</span>';
    const externalUrl = structure.alphafold_url || ('https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(structure.uniprot_id));
    document.getElementById('alphafoldDbLink').href = externalUrl;
    loading.hidden = false;
    loading.classList.remove('is-error');
    loading.textContent = 'Loading the AlphaFold structure...';
    shell.setAttribute('aria-busy', 'true');
    modal.hidden = false;
    document.body.classList.add('viewer-open');
    document.getElementById('structureViewerClose').focus();

    try {
        await loadNglViewer();
        await new Promise(function(resolve) { requestAnimationFrame(resolve); });
        if (structureStage) structureStage.dispose();
        structureStage = new window.NGL.Stage('nglViewport', { backgroundColor: '#f7faff', quality: 'medium', tooltip: true });
        const modelUrl = await resolveAlphaFoldModel(structure.uniprot_id);
        structureComponent = await structureStage.loadFile(modelUrl, { ext: 'pdb' });
        const confidenceScheme = alphaFoldConfidenceScheme();
        structureComponent.addRepresentation('cartoon', { color: confidenceScheme, smoothSheet: true, quality: 'medium' });
        const mutationSelection = new window.NGL.Selection(activeResidueSelection);
        let residueCount = 0;
        structureComponent.structure.eachAtom(function() { residueCount += 1; }, mutationSelection);
        if (residueCount > 0) {
            structureComponent.addRepresentation('ball+stick', { sele: activeResidueSelection, color: '#e53935', scale: 2.1, aspectRatio: 1.4 });
            structureComponent.addRepresentation('spacefill', { sele: activeResidueSelection, color: '#ff5252', radiusScale: .45, opacity: .85 });
            structureComponent.addRepresentation('label', { sele: activeResidueSelection + ' and .CA', labelType: 'res', color: '#a61120', backgroundColor: '#ffffff', backgroundOpacity: .86, showBackground: true, labelGrouping: 'residue', fixedSize: true, zOffset: 2 });
            structureComponent.autoView(activeResidueSelection, 650);
        } else {
            structureComponent.autoView(undefined, 650);
            document.getElementById('structureVariantLabel').textContent += ' (residue not found in model numbering)';
        }
        loading.hidden = true;
        shell.setAttribute('aria-busy', 'false');
        structureStage.handleResize();
    } catch (error) {
        loading.hidden = false;
        loading.classList.add('is-error');
        loading.textContent = 'Unable to load the interactive model. Use the AlphaFold DB link below or check network/WebGL access.';
        shell.setAttribute('aria-busy', 'false');
        console.error('Structure viewer error:', error);
    }
}

function closeStructureViewer() {
    const modal = document.getElementById('structureModal');
    if (modal.hidden) return;
    modal.hidden = true;
    document.body.classList.remove('viewer-open');
    if (structureStage) { structureStage.dispose(); structureStage = null; structureComponent = null; }
    document.getElementById('nglViewport').replaceChildren();
    if (structureViewerTrigger && document.contains(structureViewerTrigger)) structureViewerTrigger.focus();
}

document.getElementById('structureViewerClose').addEventListener('click', closeStructureViewer);
document.querySelector('[data-close-structure-viewer]').addEventListener('click', closeStructureViewer);
document.getElementById('focusResidueBtn').addEventListener('click', function() {
    if (structureComponent && activeResidueSelection) structureComponent.autoView(activeResidueSelection, 600);
});
document.getElementById('fullStructureBtn').addEventListener('click', function() {
    if (structureComponent) structureComponent.autoView(undefined, 600);
});
window.addEventListener('resize', function() { if (structureStage) structureStage.handleResize(); });
document.addEventListener('keydown', function(event) {
    if (event.key === 'Escape' && !document.getElementById('structureModal').hidden) closeStructureViewer();
    if (event.key === 'Tab' && !document.getElementById('structureModal').hidden) {
        const modal = document.getElementById('structureModal');
        const focusable = Array.from(modal.querySelectorAll('button, a[href], [tabindex]:not([tabindex="-1"])')).filter(function(el) { return !el.disabled; });
        if (!focusable.length) return;
        const first = focusable[0], last = focusable[focusable.length - 1];
        if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus(); }
        else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus(); }
    }
});
</script>
</body>
</html>
"""

if __name__ == "__main__":
    import uvicorn
    # Configurable via env vars for production deployment.
    # HOST=0.0.0.0 to accept external traffic (behind nginx), RELOAD=0 in prod.
    host = os.getenv("HOST", "127.0.0.1")
    port = int(os.getenv("PORT", "8000"))
    reload = os.getenv("RELOAD", "1") == "1"
    uvicorn.run("simple_backend:app", host=host, port=port, reload=reload)
