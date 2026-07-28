"""Versioned exact-allele matching for curated common diabetes-risk evidence.

The catalog is intentionally empty until independently verified, licensed
records are supplied. Gene-level prose is never converted into allele evidence,
and this module does not calculate a polygenic risk score.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from .pgx.models import GenomeBuild
from .pgx.normalization import normalize_variant


class DiabetesCommonRiskAnnotator:
    def __init__(self, config: Dict[str, Any], base_dir: Optional[Path] = None):
        settings = config.get("diabetes_common_risk", {})
        env_enabled = os.getenv("DIABETES_COMMON_RISK_ENABLED")
        self.enabled = (
            env_enabled.lower() in {"1", "true", "yes", "on"}
            if env_enabled is not None else bool(settings.get("enabled", True))
        )
        self.base_dir = Path(base_dir or Path.cwd()).resolve()
        configured = os.getenv("DIABETES_COMMON_RISK_CATALOG") or settings.get("catalog_path")
        self.catalog_path = None
        if configured:
            path = Path(configured).expanduser()
            self.catalog_path = path if path.is_absolute() else self.base_dir / path

    @staticmethod
    def _result(status: str, build: GenomeBuild, warning: Optional[str] = None):
        result = {
            "status": status,
            "build": build.value,
            "evaluated_alleles": 0,
            "matched_records": 0,
            "findings": [],
            "source": {},
            "warnings": [],
            "score_calculated": False,
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
            result = self._result(status, genome_build, warning)
            result["inspected_alleles"] = inspected
            return result

        if not self.enabled:
            return unavailable("disabled", "Common-risk annotation is disabled")
        if genome_build is GenomeBuild.UNKNOWN:
            return unavailable(
                "skipped_unknown_build",
                "Common-risk evidence could not perform build-specific matching because the VCF genome build is unknown",
            )
        if not self.catalog_path or not self.catalog_path.is_file():
            return unavailable(
                "missing_catalog",
                f"Inspected {inspected} sample allele(s), but no independently verified local common-risk catalog is installed",
            )

        try:
            payload = json.loads(self.catalog_path.read_text(encoding="utf-8"))
            records = payload.get("records", [])
            index = {}
            for record in records:
                record_build = GenomeBuild.parse(record.get("build"))
                normalized = normalize_variant(
                    record_build, record["chrom"], int(record["pos"]),
                    record["ref"], record["alt"],
                )
                key = (record_build.value, normalized.chrom, normalized.pos, normalized.ref, normalized.alt)
                index.setdefault(key, []).append(record)

            unique = {}
            for allele in allele_list:
                normalized = normalize_variant(
                    genome_build, allele["chrom"], int(allele["pos"]), allele["ref"], allele["alt"],
                )
                key = (genome_build.value, normalized.chrom, normalized.pos, normalized.ref, normalized.alt)
                unique[key] = allele

            findings = []
            for key, allele in unique.items():
                for record in index.get(key, []):
                    findings.append({
                        **record,
                        "sample_id": allele.get("sample_id"),
                        "genotype": allele.get("raw_gt"),
                    })
            result = self._result("complete", genome_build)
            result.update({
                "inspected_alleles": inspected,
                "evaluated_alleles": len(unique),
                "matched_records": len(findings),
                "findings": findings,
                "source": payload.get("source", {}),
            })
            return result
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            return unavailable("error", f"Common-risk catalog error: {exc}")
