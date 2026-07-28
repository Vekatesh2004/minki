"""Local, build-aware ClinVar annotation for carried VCF alleles.

Patient alleles are matched only against a local SQLite index. This module never
sends coordinates or genotypes to an external service.
"""

import os
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .pgx.models import GenomeBuild
from .pgx.normalization import normalize_variant


VALID_STATUSES = {
    "complete", "disabled", "missing_index", "skipped_unknown_build", "error",
}


class ClinVarAnnotator:
    def __init__(self, config: Dict[str, Any], base_dir: Optional[Path] = None):
        settings = config.get("clinvar", {})
        env_enabled = os.getenv("CLINVAR_ENABLED")
        self.enabled = (
            env_enabled.lower() in {"1", "true", "yes", "on"}
            if env_enabled is not None else bool(settings.get("enabled", True))
        )
        self.base_dir = Path(base_dir or Path.cwd()).resolve()
        self.indexes = dict(settings.get("indexes", {}))
        self.diabetes_only = bool(settings.get("diabetes_only", True))

    def _index_path(self, build: GenomeBuild) -> Optional[Path]:
        env_key = f"CLINVAR_{build.value.upper()}_INDEX"
        configured = os.getenv(env_key) or self.indexes.get(build.value)
        if not configured:
            return None
        path = Path(configured).expanduser()
        return path if path.is_absolute() else self.base_dir / path

    @staticmethod
    def _base_result(status: str, build: GenomeBuild, warning: Optional[str] = None):
        result = {
            "status": status,
            "build": build.value,
            "evaluated_alleles": 0,
            "matched_records": 0,
            "reportable_count": 0,
            "findings": [],
            "source": {},
            "warnings": [],
        }
        if warning:
            result["warnings"].append(warning)
        return result

    def annotate(self, build: str, carried_alleles: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
        genome_build = GenomeBuild.parse(build)
        allele_list = list(carried_alleles)
        inspected = len({
            (
                str(allele.get("chrom", "")).removeprefix("chr"),
                int(allele.get("pos") or 0),
                str(allele.get("ref", "")).upper(),
                str(allele.get("alt", "")).upper(),
            )
            for allele in allele_list
        })

        def unavailable(status: str, warning: str) -> Dict[str, Any]:
            result = self._base_result(status, genome_build, warning)
            result["inspected_alleles"] = inspected
            return result

        if not self.enabled:
            return unavailable("disabled", "Local ClinVar annotation is disabled")
        if genome_build is GenomeBuild.UNKNOWN:
            return unavailable(
                "skipped_unknown_build",
                "ClinVar could not perform build-specific matching because the VCF genome build is unknown",
            )

        index_path = self._index_path(genome_build)
        if not index_path or not index_path.is_file():
            return unavailable(
                "missing_index",
                f"Inspected {inspected} sample allele(s), but no local {genome_build.value} ClinVar index is installed",
            )

        alleles = self._normalized_unique(genome_build, allele_list)
        result = self._base_result("complete", genome_build)
        result["inspected_alleles"] = inspected
        result["evaluated_alleles"] = len(alleles)
        try:
            uri = f"file:{index_path}?mode=ro"
            with sqlite3.connect(uri, uri=True) as conn:
                conn.row_factory = sqlite3.Row
                result["source"] = self._metadata(conn)
                findings: List[Dict[str, Any]] = []
                for allele in alleles:
                    findings.extend(self._lookup(conn, genome_build, allele))
            result["findings"] = findings
            result["matched_records"] = len(findings)
            result["reportable_count"] = sum(bool(item.get("reportable")) for item in findings)
            return result
        except (sqlite3.Error, ValueError, KeyError) as exc:
            failed = self._base_result("error", genome_build, f"ClinVar index error: {exc}")
            failed["inspected_alleles"] = inspected
            failed["evaluated_alleles"] = len(alleles)
            return failed

    @staticmethod
    def _normalized_unique(build: GenomeBuild, carried_alleles: Iterable[Dict[str, Any]]):
        unique = {}
        for allele in carried_alleles:
            normalized = normalize_variant(
                build, allele["chrom"], int(allele["pos"]), allele["ref"], allele["alt"],
            )
            key = (normalized.chrom, normalized.pos, normalized.ref, normalized.alt)
            unique[key] = {
                **allele,
                "chrom": normalized.chrom,
                "pos": normalized.pos,
                "ref": normalized.ref,
                "alt": normalized.alt,
                "normalization_warnings": list(normalized.warnings),
            }
        return list(unique.values())

    @staticmethod
    def _metadata(conn: sqlite3.Connection) -> Dict[str, Any]:
        try:
            return {row["key"]: row["value"] for row in conn.execute("SELECT key, value FROM metadata")}
        except sqlite3.Error:
            return {}

    @staticmethod
    def _lookup(conn: sqlite3.Connection, build: GenomeBuild, allele: Dict[str, Any]):
        rows = conn.execute(
            """SELECT variation_id, accession, rsid, gene, significance, conditions,
                      review_status, conflict_state, last_evaluated, source_url
                 FROM clinvar_alleles
                WHERE build = ? AND chrom = ? AND pos = ? AND ref = ? AND alt = ?""",
            (build.value, allele["chrom"], allele["pos"], allele["ref"], allele["alt"]),
        ).fetchall()
        findings = []
        for row in rows:
            significance = row["significance"] or "Uncertain significance"
            lower = significance.lower()
            conflicting = bool(row["conflict_state"]) or "conflict" in lower
            reportable = not conflicting and any(
                term in lower for term in ("pathogenic", "risk factor", "association")
            ) and "uncertain" not in lower
            variation_id = row["variation_id"]
            findings.append({
                "build": build.value,
                "chrom": allele["chrom"],
                "pos": allele["pos"],
                "ref": allele["ref"],
                "alt": allele["alt"],
                "genotype": allele.get("raw_gt"),
                "sample_id": allele.get("sample_id"),
                "gene": row["gene"],
                "significance": significance,
                "conditions": row["conditions"],
                "review_status": row["review_status"],
                "conflict_state": "conflicting" if conflicting else "not_conflicting",
                "reportable": reportable,
                "variation_id": variation_id,
                "accession": row["accession"],
                "rsid": row["rsid"],
                "last_evaluated": row["last_evaluated"],
                "url": row["source_url"] or (
                    f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/"
                    if variation_id else "https://www.ncbi.nlm.nih.gov/clinvar/"
                ),
                "normalization_warnings": allele.get("normalization_warnings", []),
            })
        return findings
