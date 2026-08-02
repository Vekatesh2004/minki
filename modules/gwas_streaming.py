"""Memory-bounded GWAS/HWE processing directly from a multi-sample VCF.

The general PGx parser intentionally preserves every record and call. That is
valuable for clinical evidence, but too memory-heavy for cohort GWAS on small
servers. These functions stream one VCF row at a time and retain only bounded
Manhattan output points.
"""

from __future__ import annotations

import gzip
import math
import re
from array import array
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

from .gwas import (
    GENOME_WIDE_THRESHOLD,
    SUGGESTIVE_THRESHOLD,
    GWASInputError,
    _fit_variant,
    _hwe_exact_p,
)


def _open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open("r", encoding="utf-8")


def _detected_build(headers: List[str], requested: str) -> str:
    if requested in {"GRCh37", "GRCh38"}:
        return requested
    text = " ".join(headers)
    aliases = {"grch37": "GRCh37", "hg19": "GRCh37", "b37": "GRCh37",
               "grch38": "GRCh38", "hg38": "GRCh38", "b38": "GRCh38"}
    for token in re.findall(r"GRCh3[78]|hg(?:19|38)|b3[78]", text, re.IGNORECASE):
        resolved = aliases.get(token.lower())
        if resolved:
            return resolved
    return "unknown"


def _gt_alleles(raw_sample: str, gt_index: int) -> Optional[Tuple[int, int]]:
    values = raw_sample.split(":")
    if gt_index >= len(values):
        return None
    parts = values[gt_index].replace("|", "/").split("/")
    if len(parts) != 2 or any(value in {"", "."} for value in parts):
        return None
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return None


class _PointCollector:
    """Retain bounded, representative output plus the strongest signals."""

    def __init__(self, max_points: int):
        if max_points < 1:
            raise GWASInputError("max_points must be at least 1")
        self.max_points = max_points
        self.significant: List[Dict[str, Any]] = []
        self.ordinary: List[Dict[str, Any]] = []
        self.total = 0
        self.ordinary_seen = 0
        self.genome_wide_hits = 0
        self.suggestive_hits = 0
        # A bounded reservoir is sufficient for lambda GC and avoids retaining
        # one float per tested variant on very large VCFs.
        self.p_values = array("d")
        self._ordinary_rng = np.random.default_rng(20260731)
        self._p_value_rng = np.random.default_rng(20260732)

    def add(self, point: Dict[str, Any]) -> None:
        self.total += 1
        p_value = point["p_value"]
        if p_value <= GENOME_WIDE_THRESHOLD:
            self.genome_wide_hits += 1
        if p_value <= SUGGESTIVE_THRESHOLD:
            self.suggestive_hits += 1

        if len(self.p_values) < self.max_points:
            self.p_values.append(p_value)
        else:
            slot = int(self._p_value_rng.integers(0, self.total))
            if slot < self.max_points:
                self.p_values[slot] = p_value

        if p_value <= SUGGESTIVE_THRESHOLD:
            self.significant.append(point)
            # Keep memory bounded even when a pathological dataset has many
            # significant results. Trimming in batches avoids sorting each row.
            if len(self.significant) > 2 * self.max_points:
                self.significant = sorted(
                    self.significant, key=lambda item: item["p_value"]
                )[:self.max_points]
            return

        self.ordinary_seen += 1
        if len(self.ordinary) < self.max_points:
            self.ordinary.append(point)
            return
        # Algorithm R reservoir sampling gives each ordinary point an equal
        # chance of appearing, avoiding chromosome/file-order bias.
        slot = int(self._ordinary_rng.integers(0, self.ordinary_seen))
        if slot < self.max_points:
            self.ordinary[slot] = point

    def finish(self) -> Tuple[List[Dict[str, Any]], bool]:
        strongest = sorted(self.significant, key=lambda item: item["p_value"])
        if len(strongest) >= self.max_points:
            return strongest[:self.max_points], self.total > self.max_points
        remaining = self.max_points - len(strongest)
        return strongest + self.ordinary[:remaining], self.total > self.max_points


def _lambda_gc(p_values: array) -> Optional[float]:
    if not p_values:
        return None
    values = np.frombuffer(p_values, dtype=np.float64)
    chi2 = stats.chi2.isf(values, df=1)
    chi2 = chi2[np.isfinite(chi2)]
    if not chi2.size:
        return None
    return float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))


def _base_result(
    *, analysis: str, model: str, trait_type: str, genome_build: str,
    samples: int, tested: int, skipped_qc: int, skipped_fit: int,
    collector: _PointCollector, phenotype_column: Optional[str] = None,
    covariate_names: Optional[List[str]] = None,
) -> Dict[str, Any]:
    points, downsampled = collector.finish()
    return {
        "status": "available" if points else "not_available",
        "analysis": analysis,
        "model": model,
        "trait_type": trait_type,
        "phenotype_column": phenotype_column,
        "covariate_names": covariate_names or [],
        "genome_build": genome_build,
        "samples_analyzed": samples,
        "variants_tested": tested,
        "skipped_qc": skipped_qc,
        "skipped_fit": skipped_fit,
        "point_count": len(points),
        "original_point_count": collector.total,
        "downsampled": downsampled,
        "max_points": collector.max_points,
        "genome_wide_threshold": GENOME_WIDE_THRESHOLD,
        "suggestive_threshold": SUGGESTIVE_THRESHOLD,
        "genome_wide_hits": collector.genome_wide_hits,
        "suggestive_hits": collector.suggestive_hits,
        "lambda_gc": _lambda_gc(collector.p_values) if analysis == "gwas" else None,
        "lambda_gc_approximate": analysis == "gwas" and collector.total > len(collector.p_values),
        "p_value_fields": ["logistic" if trait_type == "binary" else "linear"] if analysis == "gwas" else ["hwe"],
        "points": points,
        "message": None if points else "No variants passed association/QC testing.",
    }


def _read_header(path: Path, requested_build: str):
    headers: List[str] = []
    samples: List[str] = []
    with _open_vcf(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#CHROM"):
                columns = line.split("\t")
                samples = columns[9:]
                break
            elif line and not line.startswith("#"):
                raise GWASInputError("The VCF is missing the #CHROM header")
    if not samples:
        raise GWASInputError("The VCF has no genotype sample columns")
    if len(samples) != len(set(samples)):
        raise GWASInputError("The VCF contains duplicate sample IDs")
    return samples, _detected_build(headers, requested_build)


def run_gwas_vcf_streaming(
    file_path: str,
    phenotype_spec: Dict[str, Any],
    requested_build: str = "auto",
    min_call_rate: float = 0.5,
    min_minor_allele_count: int = 3,
    max_points: int = 12000,
):
    """Stream a VCF and compute true phenotype association statistics."""
    path = Path(file_path)
    vcf_samples, genome_build = _read_header(path, requested_build)
    vcf_index = {sample: index for index, sample in enumerate(vcf_samples)}
    pheno_samples = phenotype_spec["samples"]
    analysis_samples = [sample for sample in pheno_samples if sample in vcf_index]
    if len(analysis_samples) < 10:
        raise GWASInputError(
            f"Only {len(analysis_samples)} sample(s) overlap between the VCF and phenotype file; "
            "a meaningful association test needs many samples (>=10 minimum, ideally hundreds+)"
        )

    trait_type = phenotype_spec["trait_type"]
    values = phenotype_spec["phenotype"]
    if trait_type == "binary":
        low, _ = phenotype_spec["distinct_values"]
        pheno_lookup = {sample: (0.0 if value == low else 1.0)
                        for sample, value in zip(pheno_samples, values)}
    else:
        pheno_lookup = dict(zip(pheno_samples, values))
    phenotype = np.array([pheno_lookup[sample] for sample in analysis_samples], dtype=float)

    covariate_names = phenotype_spec.get("covariate_names", [])
    covariates = None
    if covariate_names:
        lookup = {sample: values for sample, values in zip(
            pheno_samples, phenotype_spec.get("covariates", []))}
        covariates = np.array([lookup[sample] for sample in analysis_samples], dtype=float)

    selected_columns = [(vcf_index[sample], position)
                        for position, sample in enumerate(analysis_samples)]
    collector = _PointCollector(max_points)
    tested = skipped_qc = skipped_fit = 0

    with _open_vcf(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                skipped_qc += 1
                continue
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise GWASInputError(f"Invalid VCF position at line {line_number}") from exc
            alts = fields[4].upper().split(",")
            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                skipped_qc += len(alts)
                continue
            gt_index = format_keys.index("GT")

            for alt_index, alt in enumerate(alts, 1):
                dosages = np.full(len(analysis_samples), np.nan, dtype=float)
                for vcf_position, analysis_position in selected_columns:
                    field_index = 9 + vcf_position
                    if field_index >= len(fields):
                        continue
                    alleles = _gt_alleles(fields[field_index], gt_index)
                    if alleles is not None:
                        dosages[analysis_position] = float(
                            (alleles[0] == alt_index) + (alleles[1] == alt_index)
                        )
                observed = ~np.isnan(dosages)
                observed_count = int(observed.sum())
                call_rate = float(observed.mean()) if observed.size else 0.0
                if call_rate < min_call_rate or observed_count < 10:
                    skipped_qc += 1
                    continue
                obs_dosage = dosages[observed]
                alt_count = float(obs_dosage.sum())
                minor_count = int(min(alt_count, 2 * observed_count - alt_count))
                if minor_count < min_minor_allele_count or np.ptp(obs_dosage) == 0:
                    skipped_qc += 1
                    continue

                tested += 1
                fit = _fit_variant(
                    obs_dosage, phenotype[observed],
                    covariates[observed] if covariates is not None else None,
                    trait_type,
                )
                if fit is None:
                    skipped_fit += 1
                    continue
                beta, p_value = fit
                if not math.isfinite(p_value) or not 0 < p_value <= 1:
                    skipped_fit += 1
                    continue
                collector.add({
                    "chrom": fields[0], "pos": pos,
                    "id": None if fields[2] == "." else fields[2],
                    "ref": fields[3].upper(), "alt": alt,
                    "p_value": p_value,
                    "minus_log10_p": -math.log10(max(p_value, 1e-300)),
                    "beta": beta,
                    "maf": minor_count / (2 * observed_count),
                    "n": observed_count,
                    "field": "logistic" if trait_type == "binary" else "linear",
                })

    result = _base_result(
        analysis="gwas",
        model="logistic regression" if trait_type == "binary" else "linear regression",
        trait_type=trait_type, genome_build=genome_build,
        samples=len(analysis_samples), tested=tested,
        skipped_qc=skipped_qc, skipped_fit=skipped_fit,
        collector=collector,
        phenotype_column=phenotype_spec.get("phenotype_column"),
        covariate_names=covariate_names,
    )
    return result, {"sample_count": len(vcf_samples), "genome_build": genome_build}


def run_hwe_vcf_streaming(
    file_path: str,
    requested_build: str = "auto",
    min_call_rate: float = 0.5,
    min_samples: int = 10,
    max_points: int = 12000,
):
    """Stream a VCF and compute genotype-only HWE exact tests."""
    path = Path(file_path)
    samples, genome_build = _read_header(path, requested_build)
    n_samples = len(samples)
    if n_samples < min_samples:
        raise GWASInputError(
            f"The VCF has only {n_samples} genotyped sample(s). A Hardy-Weinberg "
            "Manhattan plot needs many samples in one multi-sample VCF (>=10 minimum). "
            "A single genome cannot yield per-variant p-values of any kind."
        )

    collector = _PointCollector(max_points)
    tested = skipped_qc = 0
    with _open_vcf(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                skipped_qc += 1
                continue
            alts = fields[4].upper().split(",")
            if len(alts) != 1:
                skipped_qc += 1
                continue
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise GWASInputError(f"Invalid VCF position at line {line_number}") from exc
            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                skipped_qc += 1
                continue
            gt_index = format_keys.index("GT")
            hom_ref = het = hom_alt = 0
            for raw_sample in fields[9:9 + n_samples]:
                alleles = _gt_alleles(raw_sample, gt_index)
                if alleles is None or any(value not in {0, 1} for value in alleles):
                    continue
                alt_copies = (alleles[0] == 1) + (alleles[1] == 1)
                if alt_copies == 0:
                    hom_ref += 1
                elif alt_copies == 2:
                    hom_alt += 1
                else:
                    het += 1
            called = hom_ref + het + hom_alt
            if called < min_samples or called / n_samples < min_call_rate:
                skipped_qc += 1
                continue
            if het + min(hom_ref, hom_alt) == 0:
                skipped_qc += 1
                continue
            tested += 1
            p_value = _hwe_exact_p(het, hom_ref, hom_alt)
            if p_value is None or not 0 < p_value <= 1:
                skipped_qc += 1
                continue
            minor = min(2 * hom_ref + het, 2 * hom_alt + het)
            collector.add({
                "chrom": fields[0], "pos": pos,
                "id": None if fields[2] == "." else fields[2],
                "ref": fields[3].upper(), "alt": alts[0],
                "p_value": p_value,
                "minus_log10_p": -math.log10(max(p_value, 1e-300)),
                "maf": minor / (2 * called), "n": called, "field": "hwe",
            })

    result = _base_result(
        analysis="hwe", model="Hardy-Weinberg equilibrium exact test",
        trait_type="none", genome_build=genome_build,
        samples=n_samples, tested=tested, skipped_qc=skipped_qc,
        skipped_fit=0, collector=collector,
    )
    return result, {"sample_count": n_samples, "genome_build": genome_build}
