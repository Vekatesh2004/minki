#!/usr/bin/env python3
"""Generate SYNTHETIC demonstration ClinVar + common-risk datasets.

These datasets exist purely to demonstrate the "ClinVar / Monogenic Diabetes"
and "Common Diabetes Risk Evidence" panels in the UI. Their alleles are matched
to the bundled example VCF (examples/sample_pharmacogenomics.vcf, GRCh38).

IMPORTANT: This is NOT real ClinVar data and NOT a validated risk catalog. The
clinical significances and effect sizes are illustrative placeholders. For real
use, build the ClinVar index with scripts/build_clinvar_index.py and supply an
independently verified common-risk catalog. Every record is labelled synthetic.

Outputs (paths match config.json):
    data/clinvar/clinvar_grch38.sqlite
    data/diabetes_common_risk.json
"""

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BUILD = "GRCh38"

# Alleles chosen to match examples/sample_pharmacogenomics.vcf exactly.
# (chrom, pos, ref, alt, rsid, gene, significance, conditions, review_status)
CLINVAR_RECORDS = [
    ("11", 17388025, "T", "C", "rs5219", "KCNJ11",
     "risk factor",
     "Type 2 diabetes mellitus; Permanent neonatal diabetes mellitus",
     "criteria provided, multiple submitters, no conflicts"),
    ("8", 117172544, "C", "T", "rs13266634", "SLC30A8",
     "risk factor",
     "Type 2 diabetes mellitus",
     "criteria provided, single submitter"),
    ("3", 12351626, "C", "G", "rs1801282", "PPARG",
     "association",
     "Type 2 diabetes mellitus; response to thiazolidinedione",
     "criteria provided, multiple submitters, no conflicts"),
    # A deliberately uncertain/benign one to show non-reportable handling.
    ("12", 21178615, "T", "C", "rs4149056", "SLCO1B1",
     "Uncertain significance",
     "Statin response; not diabetes-specific",
     "criteria provided, single submitter"),
]

# Common-risk association records (synthetic effect sizes / ancestry).
COMMON_RISK_RECORDS = [
    {
        "build": BUILD, "chrom": "11", "pos": 17388025, "ref": "T", "alt": "C",
        "rsid": "rs5219", "gene": "KCNJ11", "trait": "Type 2 diabetes",
        "effect_allele": "C", "odds_ratio": "1.14", "ancestry": "European",
        "study": "SYNTHETIC-DEMO (illustrative; not a real study)",
    },
    {
        "build": BUILD, "chrom": "8", "pos": 117172544, "ref": "C", "alt": "T",
        "rsid": "rs13266634", "gene": "SLC30A8", "trait": "Type 2 diabetes",
        "effect_allele": "C", "odds_ratio": "1.12", "ancestry": "Multi-ancestry",
        "study": "SYNTHETIC-DEMO (illustrative; not a real study)",
    },
    {
        "build": BUILD, "chrom": "3", "pos": 12351626, "ref": "C", "alt": "G",
        "rsid": "rs1801282", "gene": "PPARG", "trait": "Type 2 diabetes (protective)",
        "effect_allele": "G", "odds_ratio": "0.86", "ancestry": "European",
        "study": "SYNTHETIC-DEMO (illustrative; not a real study)",
    },
]


def build_clinvar_index(output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()
    with sqlite3.connect(output) as conn:
        conn.executescript(
            """
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE clinvar_alleles (
                build TEXT NOT NULL, chrom TEXT NOT NULL, pos INTEGER NOT NULL,
                ref TEXT NOT NULL, alt TEXT NOT NULL, variation_id TEXT,
                accession TEXT, rsid TEXT, gene TEXT, significance TEXT,
                conditions TEXT, review_status TEXT, conflict_state INTEGER NOT NULL DEFAULT 0,
                last_evaluated TEXT, source_url TEXT
            );
            CREATE INDEX clinvar_exact_allele
                ON clinvar_alleles(build, chrom, pos, ref, alt);
            """
        )
        insert = "INSERT INTO clinvar_alleles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
        for chrom, pos, ref, alt, rsid, gene, sig, conditions, review in CLINVAR_RECORDS:
            conn.execute(insert, (
                BUILD, chrom, pos, ref, alt,
                None, None, rsid, gene, sig, conditions, review, 0,
                None, "https://www.ncbi.nlm.nih.gov/clinvar/",
            ))
        metadata = {
            "source_name": "SYNTHETIC DEMO — not real ClinVar",
            "source_url": "local-demo",
            "release": "SYNTHETIC-DEMO-v1",
            "build": BUILD,
            "indexed_at": datetime.now(timezone.utc).isoformat(),
            "scope": "diabetes_related_demo",
            "record_count": str(len(CLINVAR_RECORDS)),
            "schema_version": "1.0",
        }
        conn.executemany("INSERT INTO metadata(key,value) VALUES (?,?)", metadata.items())
        conn.commit()
    return len(CLINVAR_RECORDS)


def build_common_risk_catalog(output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source": {
            "name": "SYNTHETIC DEMO common-risk catalog",
            "version": "SYNTHETIC-DEMO-v1",
            "notice": "Illustrative demonstration data only; not validated or for clinical use",
        },
        "records": COMMON_RISK_RECORDS,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return len(COMMON_RISK_RECORDS)


def main() -> None:
    clinvar_path = ROOT / "data" / "clinvar" / "clinvar_grch38.sqlite"
    risk_path = ROOT / "data" / "diabetes_common_risk.json"
    n_clinvar = build_clinvar_index(clinvar_path)
    n_risk = build_common_risk_catalog(risk_path)
    print(f"Wrote {clinvar_path} ({n_clinvar} synthetic ClinVar demo records)")
    print(f"Wrote {risk_path} ({n_risk} synthetic common-risk demo records)")
    print("These are SYNTHETIC demonstration datasets, not real clinical data.")
    print("Run the example analysis to see the ClinVar and Common Risk panels populate.")


if __name__ == "__main__":
    main()
