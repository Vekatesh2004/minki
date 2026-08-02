"""Gene expression profiling across human tissues (Objective 4, part 1).

Retrieves real median gene-level expression from the GTEx portal for the
patient's mutated genes, and scores each gene for relevance to diabetes based
on expression in metabolically-relevant tissues (pancreas, liver, skeletal
muscle, adipose).

Honesty constraints:
  * GTEx medians are population reference values from post-mortem donors. They
    are NOT the patient's own expression; this pipeline has no RNA input.
    Everything here is tissue-context annotation, not differential expression
    in this individual.
  * "Tissue enrichment" is computed relative to the gene's own median across
    tissues. It is a descriptive ratio, not a statistical DEG test.
  * When GTEx cannot resolve a gene, it is reported as unresolved rather than
    being assigned a default value.
"""

import asyncio
import logging
import statistics
from typing import Any, Dict, Iterable, List, Optional

import aiohttp

logger = logging.getLogger(__name__)

GTEX_BASE = "https://gtexportal.org/api/v2"

# Tissues where expression is most informative for type 2 diabetes biology.
DIABETES_RELEVANT_TISSUES: Dict[str, str] = {
    "Pancreas": "Insulin-secreting tissue (beta-cell context)",
    "Liver": "Hepatic glucose production and insulin clearance",
    "Muscle_Skeletal": "Primary site of insulin-stimulated glucose disposal",
    "Adipose_Subcutaneous": "Insulin sensitivity and adipokine signalling",
    "Adipose_Visceral_Omentum": "Visceral adiposity, strongly linked to insulin resistance",
    "Kidney_Cortex": "Renal glucose handling; diabetic nephropathy target tissue",
    "Artery_Coronary": "Vascular complications of diabetes",
    "Nerve_Tibial": "Peripheral neuropathy target tissue",
    "Brain_Hypothalamus": "Central appetite and energy-balance regulation",
}

# Expression bands for TPM, used for plain-language description only.
EXPRESSION_BANDS = (
    (100.0, "very high"),
    (20.0, "high"),
    (5.0, "moderate"),
    (1.0, "low"),
)


def _expression_band(tpm: float) -> str:
    for threshold, label in EXPRESSION_BANDS:
        if tpm >= threshold:
            return label
    return "minimal"


class GeneExpressionProfiler:
    """Fetches and summarises GTEx median expression for a gene set."""

    schema_version = "expression-v1"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        settings = (config or {}).get("gene_expression", {})
        self.base_url = settings.get("base_url", GTEX_BASE)
        self.dataset_id = settings.get("dataset_id", "gtex_v8")
        self.max_genes = int(settings.get("max_genes", 25))
        self.timeout = float(settings.get("timeout", 45.0))
        self.concurrency = int(settings.get("concurrency", 4))
        self.enabled = bool(settings.get("enabled", True))

    async def _get_json(
        self, session: aiohttp.ClientSession, path: str, params: Dict[str, Any],
    ) -> Any:
        async with session.get(f"{self.base_url}{path}", params=params) as response:
            if response.status != 200:
                body = await response.text()
                raise RuntimeError(f"GTEx {path} HTTP {response.status}: {body[:180]}")
            return await response.json(content_type=None)

    async def _resolve_gencode_id(
        self, session: aiohttp.ClientSession, gene: str,
    ) -> Optional[Dict[str, Any]]:
        """Map a gene symbol to a versioned GENCODE ID GTEx accepts."""
        try:
            payload = await self._get_json(session, "/reference/gene", {"geneId": gene})
        except Exception as exc:
            logger.debug("GTEx gene lookup failed for %s: %s", gene, exc)
            return None
        for record in payload.get("data", []) or []:
            # Require an exact symbol match; GTEx search is fuzzy.
            if str(record.get("geneSymbolUpper") or "").upper() == gene.upper():
                return {
                    "gencode_id": record.get("gencodeId"),
                    "description": record.get("description"),
                    "gene_type": record.get("geneType"),
                    "chromosome": record.get("chromosome"),
                }
        return None

    async def _fetch_medians(
        self, session: aiohttp.ClientSession, gencode_id: str,
    ) -> List[Dict[str, Any]]:
        payload = await self._get_json(
            session, "/expression/medianGeneExpression",
            {"gencodeId": gencode_id, "datasetId": self.dataset_id},
        )
        return payload.get("data", []) or []

    async def _profile_gene(
        self,
        session: aiohttp.ClientSession,
        gene: str,
        semaphore: asyncio.Semaphore,
    ) -> Dict[str, Any]:
        async with semaphore:
            reference = await self._resolve_gencode_id(session, gene)
            if not reference or not reference.get("gencode_id"):
                return {
                    "gene_symbol": gene,
                    "status": "unresolved",
                    "note": "GTEx could not resolve this gene symbol.",
                }
            try:
                rows = await self._fetch_medians(session, reference["gencode_id"])
            except Exception as exc:
                return {
                    "gene_symbol": gene,
                    "status": "error",
                    "note": f"GTEx expression lookup failed: {exc}",
                }

            if not rows:
                return {
                    "gene_symbol": gene,
                    "status": "no_data",
                    "gencode_id": reference["gencode_id"],
                    "note": "GTEx returned no median expression rows.",
                }

            values: Dict[str, float] = {}
            for row in rows:
                tissue = row.get("tissueSiteDetailId")
                try:
                    median = float(row.get("median"))
                except (TypeError, ValueError):
                    continue
                if tissue:
                    values[tissue] = median

            if not values:
                return {
                    "gene_symbol": gene,
                    "status": "no_data",
                    "gencode_id": reference["gencode_id"],
                    "note": "No numeric medians present in the GTEx response.",
                }

            all_medians = sorted(values.values(), reverse=True)
            overall_median = statistics.median(all_medians)
            peak_tissue = max(values.items(), key=lambda kv: kv[1])

            diabetes_tissues = []
            for tissue, rationale in DIABETES_RELEVANT_TISSUES.items():
                if tissue not in values:
                    continue
                tpm = values[tissue]
                # Fold-change vs the gene's own cross-tissue median describes
                # relative tissue specificity without any statistical claim.
                fold = (tpm / overall_median) if overall_median > 0 else None
                diabetes_tissues.append({
                    "tissue": tissue.replace("_", " "),
                    "tissue_id": tissue,
                    "median_tpm": round(tpm, 3),
                    "expression_band": _expression_band(tpm),
                    "fold_vs_gene_median": round(fold, 2) if fold is not None else None,
                    "tissue_relevance": rationale,
                })
            diabetes_tissues.sort(key=lambda t: -t["median_tpm"])

            top_tissues = [
                {
                    "tissue": tissue.replace("_", " "),
                    "median_tpm": round(tpm, 3),
                }
                for tissue, tpm in sorted(values.items(), key=lambda kv: -kv[1])[:5]
            ]

            relevance = self._relevance_score(diabetes_tissues, overall_median)

            return {
                "gene_symbol": gene,
                "status": "complete",
                "gencode_id": reference["gencode_id"],
                "gene_description": reference.get("description"),
                "gene_type": reference.get("gene_type"),
                "tissue_count": len(values),
                "overall_median_tpm": round(overall_median, 3),
                "peak_tissue": peak_tissue[0].replace("_", " "),
                "peak_tissue_tpm": round(peak_tissue[1], 3),
                "top_tissues": top_tissues,
                "diabetes_relevant_tissues": diabetes_tissues,
                "diabetes_tissue_relevance": relevance,
            }

    @staticmethod
    def _relevance_score(
        diabetes_tissues: List[Dict[str, Any]], overall_median: float,
    ) -> Dict[str, Any]:
        """Descriptive relevance summary for metabolic tissues."""
        if not diabetes_tissues:
            return {
                "level": "unknown",
                "rationale": "No diabetes-relevant tissue medians were returned.",
                "expressed_tissue_count": 0,
            }
        expressed = [t for t in diabetes_tissues if t["median_tpm"] >= 5.0]
        enriched = [
            t for t in diabetes_tissues
            if t["fold_vs_gene_median"] is not None and t["fold_vs_gene_median"] >= 2.0
        ]
        if enriched and expressed:
            level = "high"
            rationale = (
                f"Expressed in {len(expressed)} metabolic tissue(s) and enriched "
                f"(>=2x its own cross-tissue median) in "
                f"{', '.join(t['tissue'] for t in enriched[:3])}."
            )
        elif expressed:
            level = "moderate"
            rationale = (
                f"Expressed (>=5 TPM) in {len(expressed)} diabetes-relevant "
                "tissue(s) without marked tissue enrichment."
            )
        else:
            level = "low"
            rationale = (
                "Low expression (<5 TPM) across diabetes-relevant tissues in the "
                "GTEx reference population."
            )
        return {
            "level": level,
            "rationale": rationale,
            "expressed_tissue_count": len(expressed),
            "enriched_tissues": [t["tissue"] for t in enriched],
        }

    async def profile_genes(self, genes: Iterable[str]) -> Dict[str, Any]:
        gene_list: List[str] = []
        for gene in genes or []:
            if not gene:
                continue
            symbol = str(gene).strip().upper()
            if symbol and symbol not in gene_list:
                gene_list.append(symbol)
        gene_list = gene_list[: self.max_genes]

        if not self.enabled:
            return self._empty("disabled", "Gene expression profiling is disabled in config.")
        if not gene_list:
            return self._empty("no_genes", "No mutated genes were available to profile.")

        try:
            timeout = aiohttp.ClientTimeout(total=self.timeout)
            semaphore = asyncio.Semaphore(max(1, self.concurrency))
            async with aiohttp.ClientSession(timeout=timeout) as session:
                profiles = await asyncio.gather(*[
                    self._profile_gene(session, gene, semaphore) for gene in gene_list
                ], return_exceptions=True)
        except asyncio.TimeoutError:
            return self._empty(
                "timeout", f"GTEx did not respond within {self.timeout:.0f}s.",
            )
        except Exception as exc:
            logger.error("Gene expression profiling failed: %s", exc)
            return self._empty("error", f"Gene expression error: {exc}")

        resolved: List[Dict[str, Any]] = []
        warnings: List[str] = []
        for gene, profile in zip(gene_list, profiles):
            if isinstance(profile, Exception):
                warnings.append(f"{gene}: {profile}")
                continue
            if profile.get("status") == "complete":
                resolved.append(profile)
            else:
                warnings.append(f"{gene}: {profile.get('note', profile.get('status'))}")

        resolved.sort(key=lambda p: (
            {"high": 0, "moderate": 1, "low": 2, "unknown": 3}[
                p["diabetes_tissue_relevance"]["level"]
            ],
            -p["overall_median_tpm"],
        ))

        return {
            "status": "complete" if resolved else "no_data",
            "schema_version": self.schema_version,
            "dataset": self.dataset_id,
            "requested_genes": gene_list,
            "profiled_gene_count": len(resolved),
            "profiles": resolved,
            "warnings": warnings,
            "tissue_panel": [
                {"tissue_id": tid, "tissue": tid.replace("_", " "), "relevance": why}
                for tid, why in DIABETES_RELEVANT_TISSUES.items()
            ],
            "source": {
                "name": "GTEx Portal",
                "url": "https://gtexportal.org",
                "dataset": self.dataset_id,
                "unit": "median TPM",
            },
            "interpretation": (
                "GTEx medians describe expression in a healthy reference "
                "population, not this patient. No differential-expression test "
                "is performed and no per-patient RNA data is used."
            ),
        }

    def _empty(self, status: str, warning: str) -> Dict[str, Any]:
        return {
            "status": status,
            "schema_version": self.schema_version,
            "dataset": self.dataset_id,
            "requested_genes": [],
            "profiled_gene_count": 0,
            "profiles": [],
            "warnings": [warning],
            "tissue_panel": [],
            "source": {"name": "GTEx Portal", "url": "https://gtexportal.org"},
            "interpretation": None,
        }
