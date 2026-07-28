"""Genome-wide association testing from multi-sample genotypes + a phenotype.

This computes *real* per-variant association statistics. It never invents
p-values. Given a multi-sample VCF (already parsed into pgx_records) and a
phenotype table keyed by sample id, it fits one regression per tested ALT
allele and returns genome-wide p-values suitable for a true Manhattan plot.

Model selection:
  * Binary phenotype (two distinct values, e.g. control/case) -> logistic
    regression (statsmodels Logit).
  * Quantitative phenotype (>2 distinct numeric values) -> ordinary least
    squares linear regression (statsmodels OLS).

Additive genotype coding is used: dosage = number of copies of the tested ALT
allele (0, 1, 2). Optional numeric covariates are included as fixed effects.
Variants failing QC (low call rate, low minor-allele count, monomorphic) are
skipped and reported, not silently plotted.
"""

from __future__ import annotations

import csv
import io
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:  # statsmodels gives calibrated regression p-values with covariates.
    import statsmodels.api as sm
    _HAVE_SM = True
except Exception:  # pragma: no cover - exercised only without statsmodels
    _HAVE_SM = False

from scipy import stats


GENOME_WIDE_THRESHOLD = 5e-8
SUGGESTIVE_THRESHOLD = 1e-5


class GWASInputError(ValueError):
    """Raised when phenotype/genotype inputs are unusable for association."""


def _hwe_exact_p(obs_hets: int, obs_hom1: int, obs_hom2: int) -> Optional[float]:
    """Wigginton et al. (2005) exact Hardy-Weinberg equilibrium test p-value.

    Uses only genotype counts (no phenotype). Returns None when undefined
    (e.g. no minor alleles). This is a genuine statistical test of whether
    genotype proportions match HWE expectations across the cohort.
    """
    obs_homc = max(obs_hom1, obs_hom2)
    obs_homr = min(obs_hom1, obs_hom2)
    rare = 2 * obs_homr + obs_hets
    n = obs_hets + obs_homc + obs_homr
    if n == 0 or rare == 0:
        return None

    het_probs = [0.0] * (rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if mid % 2 != rare % 2:
        mid += 1
    het_probs[mid] = 1.0
    total = het_probs[mid]

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets > 1:
        het_probs[curr_hets - 2] = (
            het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
            / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        )
        total += het_probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets <= rare - 2:
        het_probs[curr_hets + 2] = (
            het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
            / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        total += het_probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    if total <= 0:
        return None
    target = het_probs[obs_hets]
    p_value = sum(prob for prob in het_probs if prob <= target) / total
    return float(min(max(p_value, 0.0), 1.0))


def run_genotype_only_manhattan(
    pgx_records: Sequence[Dict[str, Any]],
    genome_build: str = "unknown",
    min_call_rate: float = 0.5,
    min_samples: int = 10,
    max_points: int = 200000,
) -> Dict[str, Any]:
    """Compute a Hardy-Weinberg-equilibrium Manhattan plot from genotypes alone.

    This needs NO phenotype. It tests, per variant, whether observed genotype
    proportions across the cohort deviate from HWE expectations. It is a
    genotyping/population-structure QC view, not a trait-association result:
    peaks flag possible genotyping error, batch effects, or population
    structure, and must never be read as disease significance.
    """
    sample_ids = set()
    for record in pgx_records:
        for call in record.get("calls", []) or []:
            if call.get("sample_id"):
                sample_ids.add(call["sample_id"])
    n_samples = len(sample_ids)
    if n_samples < min_samples:
        raise GWASInputError(
            f"The VCF has only {n_samples} genotyped sample(s). A Hardy-Weinberg "
            "Manhattan plot needs many samples in one multi-sample VCF (>=10 minimum). "
            "A single genome cannot yield per-variant p-values of any kind."
        )

    points: List[Dict[str, Any]] = []
    tested = 0
    skipped_qc = 0
    for record in pgx_records:
        alts = record.get("alts", []) or []
        calls = record.get("calls", []) or []
        # Only diallelic sites have a well-defined single-ALT HWE test here.
        if len(alts) != 1:
            skipped_qc += 1
            continue
        hom_ref = het = hom_alt = missing = 0
        for call in calls:
            indices = call.get("allele_indices")
            if not indices or any(value is None for value in indices):
                missing += 1
                continue
            alt_copies = sum(1 for value in indices if value == 1)
            ref_copies = sum(1 for value in indices if value == 0)
            if alt_copies + ref_copies != len(indices):
                missing += 1  # multiallelic call at a diallelic record
                continue
            if alt_copies == 0:
                hom_ref += 1
            elif alt_copies == len(indices):
                hom_alt += 1
            else:
                het += 1
        called = hom_ref + het + hom_alt
        call_rate = called / n_samples if n_samples else 0.0
        if called < min_samples or call_rate < min_call_rate:
            skipped_qc += 1
            continue
        if (het + min(hom_ref, hom_alt)) == 0:
            skipped_qc += 1  # monomorphic
            continue
        p_value = _hwe_exact_p(het, hom_ref, hom_alt)
        tested += 1
        if p_value is None or not (0 < p_value <= 1):
            skipped_qc += 1
            continue
        minor = min(2 * hom_ref + het, 2 * hom_alt + het)
        points.append({
            "chrom": str(record.get("chrom") or ""),
            "pos": int(record.get("pos") or 0),
            "id": record.get("id"),
            "ref": record.get("ref"),
            "alt": alts[0],
            "p_value": p_value,
            "minus_log10_p": -math.log10(max(p_value, 1e-300)),
            "maf": float(minor) / (2 * called) if called else 0.0,
            "n": called,
            "field": "hwe",
        })

    original_count = len(points)
    downsampled = False
    if len(points) > max_points:
        significant = [p for p in points if p["p_value"] <= SUGGESTIVE_THRESHOLD]
        ordinary = [p for p in points if p["p_value"] > SUGGESTIVE_THRESHOLD]
        remaining = max(0, max_points - len(significant))
        if remaining and ordinary:
            step = max(1, math.ceil(len(ordinary) / remaining))
            ordinary = ordinary[::step][:remaining]
        else:
            ordinary = []
        points = significant + ordinary
        downsampled = True

    genome_wide_hits = sum(1 for point in points if point["p_value"] <= GENOME_WIDE_THRESHOLD)
    suggestive_hits = sum(1 for point in points if point["p_value"] <= SUGGESTIVE_THRESHOLD)

    return {
        "status": "available" if points else "not_available",
        "analysis": "hwe",
        "model": "Hardy-Weinberg equilibrium exact test",
        "trait_type": "none",
        "phenotype_column": None,
        "covariate_names": [],
        "genome_build": genome_build,
        "samples_analyzed": n_samples,
        "variants_tested": tested,
        "skipped_qc": skipped_qc,
        "skipped_fit": 0,
        "point_count": len(points),
        "original_point_count": original_count,
        "downsampled": downsampled,
        "max_points": max_points,
        "genome_wide_threshold": GENOME_WIDE_THRESHOLD,
        "suggestive_threshold": SUGGESTIVE_THRESHOLD,
        "genome_wide_hits": genome_wide_hits,
        "suggestive_hits": suggestive_hits,
        "lambda_gc": None,
        "p_value_fields": ["hwe"],
        "points": points,
        "message": None if points else (
            "No diallelic variants passed call-rate and polymorphism QC for a "
            "Hardy-Weinberg test. Check that the VCF has many genotyped samples."
        ),
    }


def parse_phenotype_table(
    text: str,
    phenotype_column: Optional[str] = None,
    sample_column: Optional[str] = None,
) -> Dict[str, Any]:
    """Parse a CSV/TSV phenotype table into samples, phenotype, and covariates.

    The first column (or ``sample_column``) identifies the sample. The
    ``phenotype_column`` (or a column literally named 'phenotype'/'trait', else
    the second column) holds the trait. Remaining numeric columns become
    covariates. Missing phenotype values drop that sample from the analysis.
    """
    sniff = text.lstrip()
    if not sniff:
        raise GWASInputError("The phenotype file is empty")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    rows = [row for row in reader if row and any(cell.strip() for cell in row)]
    if len(rows) < 2:
        raise GWASInputError("The phenotype file needs a header and at least one sample row")

    header = [cell.strip() for cell in rows[0]]
    lower = [name.lower() for name in header]

    sample_idx = 0
    if sample_column:
        if sample_column not in header:
            raise GWASInputError(f"Sample column '{sample_column}' is not present")
        sample_idx = header.index(sample_column)
    else:
        for candidate in ("sample_id", "sample", "iid", "id"):
            if candidate in lower:
                sample_idx = lower.index(candidate)
                break

    if phenotype_column:
        if phenotype_column not in header:
            raise GWASInputError(f"Phenotype column '{phenotype_column}' is not present")
        pheno_idx = header.index(phenotype_column)
    else:
        pheno_idx = None
        for candidate in ("phenotype", "pheno", "trait", "status", "y"):
            if candidate in lower:
                pheno_idx = lower.index(candidate)
                break
        if pheno_idx is None:
            pheno_idx = 1 if len(header) > 1 else None
        if pheno_idx is None or pheno_idx == sample_idx:
            raise GWASInputError("Could not identify a phenotype column")

    covariate_indices = [
        index for index in range(len(header))
        if index not in (sample_idx, pheno_idx)
    ]

    samples: List[str] = []
    phenotype: List[Optional[float]] = []
    covariate_rows: List[List[Optional[float]]] = []
    for row in rows[1:]:
        if len(row) <= max(sample_idx, pheno_idx):
            continue
        sample = row[sample_idx].strip()
        if not sample:
            continue
        raw_pheno = row[pheno_idx].strip()
        if raw_pheno in ("", ".", "NA", "NaN", "nan", "-9"):
            continue
        try:
            pheno_value = float(raw_pheno)
        except ValueError:
            raise GWASInputError(
                f"Phenotype value '{raw_pheno}' for sample '{sample}' is not numeric"
            )
        covariates: List[Optional[float]] = []
        for index in covariate_indices:
            cell = row[index].strip() if index < len(row) else ""
            if cell in ("", ".", "NA", "NaN", "nan", "-9"):
                covariates.append(None)
            else:
                try:
                    covariates.append(float(cell))
                except ValueError:
                    covariates.append(None)
        samples.append(sample)
        phenotype.append(pheno_value)
        covariate_rows.append(covariates)

    if not samples:
        raise GWASInputError("No samples with a usable phenotype value were found")

    # Keep only covariate columns that are fully numeric across kept samples.
    usable_cov: List[str] = []
    usable_cov_indices: List[int] = []
    for local_index, header_index in enumerate(covariate_indices):
        column = [row[local_index] for row in covariate_rows]
        if all(value is not None for value in column) and len({v for v in column}) > 1:
            usable_cov.append(header[header_index])
            usable_cov_indices.append(local_index)

    covariate_matrix = [
        [row[local_index] for local_index in usable_cov_indices]
        for row in covariate_rows
    ]

    distinct = sorted({value for value in phenotype})
    trait_type = "binary" if len(distinct) == 2 else "quantitative"

    return {
        "samples": samples,
        "phenotype": phenotype,
        "phenotype_column": header[pheno_idx],
        "covariate_names": usable_cov,
        "covariates": covariate_matrix,
        "trait_type": trait_type,
        "distinct_values": distinct,
    }


def _dosage_for_call(call: Dict[str, Any], alt_index: int) -> Optional[float]:
    """Additive dosage: number of copies of the given ALT allele index."""
    indices = call.get("allele_indices")
    if not indices:
        return None
    if any(value is None for value in indices):
        return None  # partial/no-call -> missing genotype
    return float(sum(1 for value in indices if value == alt_index))


def _fit_variant(
    dosage: np.ndarray,
    phenotype: np.ndarray,
    covariates: Optional[np.ndarray],
    trait_type: str,
) -> Optional[Tuple[float, float]]:
    """Return (beta, p_value) for the genotype term, or None if the fit fails."""
    columns = [np.ones_like(dosage), dosage]
    if covariates is not None and covariates.size:
        columns.extend(covariates[:, index] for index in range(covariates.shape[1]))
    design = np.column_stack(columns)

    if _HAVE_SM:
        try:
            if trait_type == "binary":
                model = sm.Logit(phenotype, design)
                result = model.fit(disp=0, maxiter=100)
            else:
                model = sm.OLS(phenotype, design)
                result = model.fit()
            beta = float(result.params[1])
            p_value = float(result.pvalues[1])
            if not math.isfinite(p_value):
                return None
            return beta, p_value
        except Exception:
            return None

    # Fallback: linear model via numpy least squares + t-test on the slope.
    try:
        coef, residuals, rank, _ = np.linalg.lstsq(design, phenotype, rcond=None)
        n, k = design.shape
        if n <= k:
            return None
        fitted = design @ coef
        resid = phenotype - fitted
        dof = n - k
        sigma2 = float(resid @ resid) / dof
        xtx_inv = np.linalg.inv(design.T @ design)
        se = math.sqrt(sigma2 * xtx_inv[1, 1])
        if se == 0:
            return None
        t_stat = coef[1] / se
        p_value = float(2 * stats.t.sf(abs(t_stat), dof))
        if not math.isfinite(p_value):
            return None
        return float(coef[1]), p_value
    except Exception:
        return None


def run_gwas(
    pgx_records: Sequence[Dict[str, Any]],
    phenotype_spec: Dict[str, Any],
    genome_build: str = "unknown",
    min_call_rate: float = 0.5,
    min_minor_allele_count: int = 3,
    max_points: int = 200000,
) -> Dict[str, Any]:
    """Compute per-variant association statistics for a true Manhattan plot."""
    pheno_samples = phenotype_spec["samples"]
    trait_type = phenotype_spec["trait_type"]
    phenotype_values = phenotype_spec["phenotype"]

    # Binary traits are recoded to 0/1 in the order of their sorted values.
    if trait_type == "binary":
        low, high = phenotype_spec["distinct_values"]
        pheno_lookup = {sample: (0.0 if value == low else 1.0)
                        for sample, value in zip(pheno_samples, phenotype_values)}
    else:
        pheno_lookup = dict(zip(pheno_samples, phenotype_values))

    covariate_names = phenotype_spec.get("covariate_names", [])
    covariate_lookup = {
        sample: values
        for sample, values in zip(pheno_samples, phenotype_spec.get("covariates", []))
    }

    # Intersect phenotype samples with samples that actually appear in the VCF.
    vcf_samples = set()
    for record in pgx_records:
        for call in record.get("calls", []) or []:
            if call.get("sample_id"):
                vcf_samples.add(call["sample_id"])
    analysis_samples = [s for s in pheno_samples if s in vcf_samples]
    if len(analysis_samples) < 10:
        raise GWASInputError(
            f"Only {len(analysis_samples)} sample(s) overlap between the VCF and phenotype file; "
            "a meaningful association test needs many samples (>=10 minimum, ideally hundreds+)"
        )

    sample_order = {sample: index for index, sample in enumerate(analysis_samples)}
    phenotype = np.array([pheno_lookup[s] for s in analysis_samples], dtype=float)
    covariates = None
    if covariate_names:
        covariates = np.array([covariate_lookup[s] for s in analysis_samples], dtype=float)

    points: List[Dict[str, Any]] = []
    tested = 0
    skipped_qc = 0
    skipped_fit = 0
    fields_seen = {"linear" if trait_type == "quantitative" else "logistic"}

    for record in pgx_records:
        alts = record.get("alts", []) or []
        calls_by_sample = {
            call.get("sample_id"): call for call in record.get("calls", []) or []
        }
        for alt_index, alt in enumerate(alts, start=1):
            dosages = np.full(len(analysis_samples), np.nan, dtype=float)
            for sample, position in sample_order.items():
                call = calls_by_sample.get(sample)
                if call is None:
                    continue
                dosage = _dosage_for_call(call, alt_index)
                if dosage is not None:
                    dosages[position] = dosage

            observed = ~np.isnan(dosages)
            call_rate = float(observed.mean()) if observed.size else 0.0
            if call_rate < min_call_rate or observed.sum() < 10:
                skipped_qc += 1
                continue

            obs_dosage = dosages[observed]
            minor_allele_count = int(min(obs_dosage.sum(), (2 * observed.sum()) - obs_dosage.sum()))
            if minor_allele_count < min_minor_allele_count or len({*obs_dosage.tolist()}) < 2:
                skipped_qc += 1  # monomorphic or too rare
                continue

            y = phenotype[observed]
            cov = covariates[observed] if covariates is not None else None
            fit = _fit_variant(obs_dosage, y, cov, trait_type)
            tested += 1
            if fit is None:
                skipped_fit += 1
                continue
            beta, p_value = fit
            if not (0 < p_value <= 1):
                skipped_fit += 1
                continue

            points.append({
                "chrom": str(record.get("chrom") or ""),
                "pos": int(record.get("pos") or 0),
                "id": record.get("id"),
                "ref": record.get("ref"),
                "alt": alt,
                "p_value": p_value,
                "minus_log10_p": -math.log10(max(p_value, 1e-300)),
                "beta": beta,
                "maf": float(minor_allele_count) / (2 * int(observed.sum())),
                "n": int(observed.sum()),
                "field": "logistic" if trait_type == "binary" else "linear",
            })

    # Genomic inflation factor (lambda_GC): median chi-square vs expected.
    lambda_gc = None
    if points:
        chi2 = np.array([stats.chi2.isf(point["p_value"], df=1) for point in points])
        chi2 = chi2[np.isfinite(chi2)]
        if chi2.size:
            lambda_gc = float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))

    original_count = len(points)
    downsampled = False
    if len(points) > max_points:
        significant = [p for p in points if p["p_value"] <= SUGGESTIVE_THRESHOLD]
        ordinary = [p for p in points if p["p_value"] > SUGGESTIVE_THRESHOLD]
        remaining = max(0, max_points - len(significant))
        if remaining and ordinary:
            step = max(1, math.ceil(len(ordinary) / remaining))
            ordinary = ordinary[::step][:remaining]
        else:
            ordinary = []
        points = significant + ordinary
        downsampled = True

    genome_wide_hits = sum(1 for point in points if point["p_value"] <= GENOME_WIDE_THRESHOLD)
    suggestive_hits = sum(1 for point in points if point["p_value"] <= SUGGESTIVE_THRESHOLD)

    return {
        "status": "available" if points else "not_available",
        "analysis": "gwas",
        "model": "logistic regression" if trait_type == "binary" else "linear regression",
        "trait_type": trait_type,
        "phenotype_column": phenotype_spec.get("phenotype_column"),
        "covariate_names": covariate_names,
        "genome_build": genome_build,
        "samples_analyzed": len(analysis_samples),
        "variants_tested": tested,
        "skipped_qc": skipped_qc,
        "skipped_fit": skipped_fit,
        "point_count": len(points),
        "original_point_count": original_count,
        "downsampled": downsampled,
        "max_points": max_points,
        "genome_wide_threshold": GENOME_WIDE_THRESHOLD,
        "suggestive_threshold": SUGGESTIVE_THRESHOLD,
        "genome_wide_hits": genome_wide_hits,
        "suggestive_hits": suggestive_hits,
        "lambda_gc": lambda_gc,
        "p_value_fields": sorted(fields_seen),
        "points": points,
        "message": None if points else (
            "No variants passed QC and association testing. Check sample overlap, "
            "genotype quality, and that the phenotype has variation."
        ),
    }
