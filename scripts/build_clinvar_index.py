#!/usr/bin/env python3
"""Build a local exact-allele ClinVar SQLite index from an official VCF.

The input is streamed from a local .vcf/.vcf.gz file or HTTPS URL. By default,
only diabetes-related genes/conditions are retained to keep phase-one indexes
small. No patient data is read or transmitted by this command.
"""

import argparse
import gzip
import hashlib
import io
import json
import sqlite3
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Optional, TextIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from modules.diabetes_kb import DIABETES_GENES
from modules.pgx.models import GenomeBuild
from modules.pgx.normalization import normalize_variant

DEFAULT_URLS = {
    "GRCh37": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    "GRCh38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}
DIABETES_TERMS = (
    "diabetes", "mody", "hyperglyc", "insulin resistance", "wolfram",
    "pancreatic agenesis", "neonatal diabetes",
)


def parse_info(raw: str) -> Dict[str, str]:
    fields = {}
    for item in raw.split(";"):
        key, sep, value = item.partition("=")
        fields[key] = value if sep else "true"
    return fields


def allele_value(value: Optional[str], index: int, alt_count: int) -> Optional[str]:
    if value is None:
        return None
    values = value.split(",")
    return values[index] if len(values) == alt_count else value


def diabetes_related(info: Dict[str, str]) -> bool:
    genes = info.get("GENEINFO", "")
    symbols = {entry.split(":", 1)[0].upper() for entry in genes.split("|") if entry}
    if symbols.intersection(DIABETES_GENES):
        return True
    condition_text = " ".join((info.get("CLNDN", ""), info.get("CLNDISDB", ""))).lower()
    return any(term in condition_text for term in DIABETES_TERMS)

def create_schema(conn: sqlite3.Connection):
    conn.executescript("""
        PRAGMA journal_mode=DELETE;
        PRAGMA synchronous=NORMAL;
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
    """)


def index_stream(handle: TextIO, conn: sqlite3.Connection, build: GenomeBuild, diabetes_only: bool) -> int:
    count = 0
    insert = """INSERT INTO clinvar_alleles VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
    for line_number, line in enumerate(handle, 1):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 8:
            raise ValueError(f"Malformed VCF line {line_number}")
        chrom, pos_raw, variation_id, ref, alts_raw, _, _, info_raw = fields[:8]
        info = parse_info(info_raw)
        if diabetes_only and not diabetes_related(info):
            continue
        alts = alts_raw.split(",")
        for index, alt in enumerate(alts):
            if alt.startswith("<") or "[" in alt or "]" in alt or alt == "*":
                continue
            normalized = normalize_variant(build, chrom, int(pos_raw), ref, alt)
            significance = allele_value(info.get("CLNSIG"), index, len(alts))
            conflict = int("conflict" in str(significance or "").lower())
            variation_id = None if variation_id == "." else variation_id
            accession = f"VCV{int(variation_id):09d}" if variation_id and variation_id.isdigit() else None
            rs_number = allele_value(info.get("RS"), index, len(alts))
            rsid = f"rs{rs_number}" if rs_number and not str(rs_number).lower().startswith("rs") else rs_number
            gene = info.get("GENEINFO")
            source_url = (
                f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"
                if variation_id else "https://www.ncbi.nlm.nih.gov/clinvar/"
            )
            conn.execute(insert, (
                build.value, normalized.chrom, normalized.pos, normalized.ref, normalized.alt,
                variation_id, accession, rsid, gene,
                significance, allele_value(info.get("CLNDN"), index, len(alts)),
                allele_value(info.get("CLNREVSTAT"), index, len(alts)), conflict,
                None, source_url,
            ))
            count += 1
        if count and count % 10000 == 0:
            conn.commit()
    conn.commit()
    return count


def open_source(source: str):
    if source.startswith("https://"):
        response = urllib.request.urlopen(source, timeout=120)
        stream = gzip.GzipFile(fileobj=response) if source.endswith(".gz") else response
        return io.TextIOWrapper(stream, encoding="utf-8"), response
    path = Path(source).expanduser()
    binary = path.open("rb")
    stream = gzip.GzipFile(fileobj=binary) if path.name.endswith(".gz") else binary
    return io.TextIOWrapper(stream, encoding="utf-8"), binary

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", required=True, choices=("GRCh37", "GRCh38"))
    parser.add_argument("--source", help="Local ClinVar VCF(.gz) or official HTTPS URL")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--all-conditions", action="store_true", help="Index all exact sequence alleles")
    parser.add_argument("--release", default="unknown")
    args = parser.parse_args()
    build = GenomeBuild.parse(args.build)
    source = args.source or DEFAULT_URLS[args.build]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(dir=args.output.parent, suffix=".sqlite", delete=False) as tmp:
        temp_path = Path(tmp.name)
    try:
        with sqlite3.connect(temp_path) as conn:
            create_schema(conn)
            handle, owner = open_source(source)
            try:
                count = index_stream(handle, conn, build, not args.all_conditions)
            finally:
                handle.close()
                try:
                    owner.close()
                except Exception:
                    pass
            metadata = {
                "source_name": "NCBI ClinVar",
                "source_url": source,
                "release": args.release,
                "build": build.value,
                "indexed_at": datetime.now(timezone.utc).isoformat(),
                "scope": "all_conditions" if args.all_conditions else "diabetes_related",
                "record_count": str(count),
                "schema_version": "1.0",
            }
            conn.executemany("INSERT INTO metadata(key, value) VALUES (?, ?)", metadata.items())
            conn.commit()
        temp_path.replace(args.output)
        print(json.dumps({"output": str(args.output), "records": count, "build": build.value}, indent=2))
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()
