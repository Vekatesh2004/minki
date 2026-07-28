    # HelixAI / Minki — Pharmacogenomics Analysis Platform
### Client Presentation Report

*Research and educational use only. Outputs are not a medical diagnosis, prescribing advice, or proof of clinical suitability. All findings should be confirmed against primary evidence sources.*

---

## 1. Executive summary

HelixAI (internal name: Minki) is a precision-medicine analysis platform that turns a genomic variant file (VCF) into interpretable molecular, structural, pharmacological, and statistical evidence. A user uploads a VCF through the web application; the pipeline parses and quality-controls the variants, annotates them with gene and protein consequences, maps them onto 3D protein structures, matches them to drug-interaction evidence, adds curated diabetes context, and — when a multi-sample cohort is supplied — computes a genuine genome-wide association (GWAS) Manhattan plot.

The platform is built on open, authoritative bioinformatics data sources (Ensembl, UniProt, AlphaFold, ClinPGx/PharmGKB, ClinVar) and performs real statistical computation rather than illustrative mock-ups.

---

## 2. What the platform delivers

| Capability | Description |
|---|---|
| Variant parsing + QC | Parses VCF/VCF.gz, extracts variants, genotypes, genome build, and computes QC metrics (PASS rate, QUAL, depth, Ts/Tv ratio, genotype counts). |
| Functional annotation | Annotates variants via the Ensembl VEP REST API: gene, consequence, HGVSc/HGVSp, protein residue. |
| 3D protein structure | Resolves UniProt IDs and retrieves AlphaFold predicted structures; renders an interactive, pLDDT-colored 3D viewer with the mutated residue highlighted. |
| Drug-interaction evidence | Matches affected genes to ClinPGx/PharmGKB clinical annotations with evidence-level grading; optional DrugBank enrichment. |
| Diabetes context | Curated knowledge-base gene context, plus optional local ClinVar and common-risk evidence layers. |
| GWAS Manhattan plot | For multi-sample cohorts: real per-variant association statistics (linear/logistic regression) or a genotype-only Hardy-Weinberg plot when no phenotype is available. |
| Interactive dashboard | Summary metrics, mutation table, mutation-type donut chart, 3D structure viewer, Manhattan plot, and CSV export. |

---

## 3. How it works (pipeline flow)

See the diagrams in this folder:
- `MINKI_CLIENT_FLOW.png` — high-level flow for presentation.
- `MINKI_PROJECT_DAG.png` — full audited architecture (all entry points and components).

Step-by-step (active application, `simple_backend.py`):

1. **Upload** — User submits a VCF plus a sample identifier and optional genome-build selection. A bundled synthetic example is available for demonstration.
2. **Parse + QC** — The VCF is parsed; QC metrics and genotype/allele structure are extracted.
3. **Variant selection** — PASS variants are selected (capped for reliable API throughput; default 300).
4. **Annotation (Ensembl VEP)** — Selected variants are annotated in batches; coding mutations are extracted.
5. **Structure (UniProt + AlphaFold)** — Up to 15 unique coding genes are mapped to UniProt IDs and AlphaFold structures.
6. **Drug matching (ClinPGx/PharmGKB)** — Affected genes are matched to clinical drug annotations with evidence grading.
7. **Diabetes context** — Curated gene-level context is added; optional ClinVar/common-risk layers report inspected alleles.
8. **GWAS / Manhattan (optional)** — A multi-sample VCF (with or without a phenotype table) produces a true Manhattan plot.
9. **Results** — Assembled into an interactive dashboard; the browser polls until completion and renders all evidence.

---

## 4. GWAS / Manhattan plot — scientific note

A single-sample VCF contains only genotypes and coordinates; it has **no** association p-values. A phenotype is an external observed trait and cannot be derived from a VCF — inventing one from the genotypes and then testing against it would be statistically invalid. The platform therefore offers two scientifically valid modes, both requiring a **multi-sample** VCF:

- **Trait association (with a phenotype table):** additive genotype dosage per sample is regressed against the trait. Binary traits use logistic regression; quantitative traits use linear regression (via statsmodels, with a scipy fallback). Real p-values are computed with per-variant QC (call rate, minor-allele count, monomorphic filtering) and a genomic-inflation factor (λGC).
- **Genotype-only (no phenotype):** a Hardy-Weinberg equilibrium exact test computed from genotype counts alone. This is genotyping/population-structure QC, clearly labelled and **not** trait or disease significance.

Both modes plot genuine −log10(p) with the 5×10⁻⁸ genome-wide and 1×10⁻⁵ suggestive thresholds. Demonstration datasets can be generated with `scripts/make_gwas_example.py` and `scripts/make_hwe_example.py`.

---

## 5. External APIs and data sources

All external calls are to public, authoritative bioinformatics services.

| Service | Purpose | Endpoint (base) |
|---|---|---|
| Ensembl VEP REST | Variant effect prediction (gene, consequence, HGVS, residue) | `https://rest.ensembl.org/vep/human/region` |
| UniProt | Protein identifier resolution and gene mapping | `https://rest.uniprot.org` (+ `/idmapping`) |
| AlphaFold (EBI) | Predicted 3D protein structures and PDB files | `https://alphafold.ebi.ac.uk/files`, `/api/prediction/{id}` |
| ClinPGx / PharmGKB | Clinical drug–gene annotations and evidence levels | `https://api.clinpgx.org/v1` |
| DrugBank *(optional)* | Additional drug data (requires a valid API key) | `https://go.drugbank.com/api/v1` |
| NCBI ClinVar *(optional)* | Local build-specific clinical-significance index | Built from ClinVar VCF; served locally |
| 3Dmol.js / NGL (CDN) | In-browser 3D structure rendering | `https://3Dmol.csb.pitt.edu`, cdnjs |

Reference links surfaced to users (not queried programmatically): PharmGKB, OMIM, and the GWAS Catalog (EBI).

---

## 6. Active application endpoints

The running product (`simple_backend.py`) exposes:

| Method | Path | Purpose |
|---|---|---|
| GET | `/` | Web application (single-page UI) |
| GET | `/health` | Health check and runtime configuration |
| POST | `/api/upload` | Upload a VCF (+ sample ID, genome build) and queue analysis |
| POST | `/api/example` | Run the bundled synthetic example |
| POST | `/api/gwas` | Compute a Manhattan plot from a multi-sample VCF (optional phenotype) |
| GET | `/api/analysis/{id}` | Poll analysis status/progress and retrieve results |
| GET | `/api/analysis/{id}/mutations.csv` | Download identified mutations as CSV |

---

## 7. Technology stack

- **Backend:** Python, FastAPI, Uvicorn; in-process background tasks for job execution.
- **Scientific computing:** NumPy, SciPy, statsmodels (GWAS association and HWE testing); optional cyvcf2 for accelerated VCF parsing (basic parser fallback otherwise).
- **Frontend (active):** Self-contained HTML/CSS/JavaScript served by the backend; interactive charts drawn on HTML canvas; 3D structures via 3Dmol.js/NGL.
- **Data/storage:** Local SQLite drug cache; per-job uploads directory; optional local ClinVar SQLite index.
- **Deployment:** Local (`127.0.0.1:8000`); optional EC2 with nginx and systemd supervision.

A separate production-oriented FastAPI service (`backend/`) with JWT auth, PostgreSQL, and metrics, an MCP stdio entry point (`pharmacogenomics_mcp.py`), a legacy Flask app, and a React/Docker stack exist in the repository as additional or intended architecture; the audited DAG marks their exact status.

---

## 8. Security, privacy, and compliance posture

- Uploaded VCFs, results, databases, secrets, and reference data are Git-ignored and kept local.
- Patient coordinates are **not** sent to ClinVar; local matching is used for that layer.
- The active local deployment has no built-in authentication — network exposure (e.g. the optional EC2 path) should be placed behind access controls before handling any sensitive data.
- All outputs carry a research/education-only disclaimer.

---

## 9. Current status and roadmap

**Completed and active:** VCF parsing + QC, Ensembl VEP annotation, AlphaFold structure viewer, ClinPGx/PharmGKB drug matching, curated diabetes context, GWAS/HWE Manhattan plotting, interactive dashboard with CSV export, genome-build selection, and inspected-allele reporting for optional evidence layers.

**Implemented foundations (not yet wired into the active dashboard):** local ClinVar annotator and index builder, common-risk annotator, build-aware normalization, and versioned definition validation.

**Planned:** deterministic star-allele calling, activity-score/phenotype mapping, versioned CPIC/DPWG guidance, CNV/SV evidence, and completion of the React/Docker production stack.

---

## 10. Suggested demonstration script

1. Open the web app and click **Run Example Analysis** to show the full annotation → structure → drug → diabetes flow end-to-end.
2. Open the interactive 3D structure viewer and highlight a mutated residue.
3. In the GWAS panel, upload `examples/hwe_example_multisample.vcf` with no phenotype to show a genotype-only Manhattan plot.
4. Upload `examples/gwas_example_genotypes.vcf` with `examples/gwas_example_phenotype.csv` to show a true trait-association Manhattan plot with genome-wide peaks.
5. Download the mutations CSV to show export.
