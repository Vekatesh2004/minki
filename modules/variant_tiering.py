"""Deterministic tier assignment for T2D-susceptible variants (Objective 1).

Segregates annotated variants into transparent tiers for early-detection
triage. Every tier decision records the exact rules that fired, so a reviewer
can audit why a variant landed where it did.

Scope and honesty constraints:
  * This is research prioritisation, NOT ACMG/AMP clinical classification and
    NOT a diagnosis. No pathogenicity verdict is invented.
  * A tier is only as strong as the evidence layers that were actually
    available. When ClinVar or the common-risk catalog is missing, disabled or
    seeded with demonstration data, the tier is explicitly marked as
    evidence-limited rather than silently upgraded.
  * No polygenic risk score is computed.
"""

from typing import Any, Dict, Iterable, List, Optional, Tuple

from . import diabetes_kb

# ---------------------------------------------------------------------------
# Consequence severity. Ensembl VEP terms, ordered by protein impact.
# ---------------------------------------------------------------------------
LOF_CONSEQUENCES = frozenset({
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification",
})

MODERATE_CONSEQUENCES = frozenset({
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
})

LOW_CONSEQUENCES = frozenset({
    "splice_region_variant", "synonymous_variant", "start_retained_variant",
    "stop_retained_variant", "incomplete_terminal_codon_variant",
})

# Diabetes KB categories that represent single-gene (monogenic) disease, where
# a damaging allele carries far more interpretive weight than at a common locus.
MONOGENIC_CATEGORY_MARKERS = (
    "mody", "neonatal", "wolfram", "severe insulin resistance",
    "pancreatic agenesis",
)

# ClinPGx/PharmGKB evidence levels considered clinically actionable.
STRONG_PGX_LEVELS = ("1a", "1b")
MODERATE_PGX_LEVELS = ("2a", "2b")

TIER_LABELS = {
    1: "Tier 1 — Strong evidence, review first",
    2: "Tier 2 — Moderate evidence, likely relevant",
    3: "Tier 3 — Uncertain significance, supporting only",
    4: "Tier 4 — Low priority / no diabetes evidence",
}

TIER_DESCRIPTIONS = {
    1: ("Damaging allele in an established monogenic diabetes gene, or a "
        "non-conflicting pathogenic/risk ClinVar assertion."),
    2: ("Protein-altering allele in a diabetes gene with structural or "
        "pharmacogenomic support, or a curated common-risk association."),
    3: ("Coding change in a diabetes-associated gene without strong "
        "corroborating evidence. Requires manual review."),
    4: "No established diabetes evidence at this allele.",
}


def _is_monogenic_gene(kb_record: Optional[Dict[str, Any]]) -> bool:
    if not kb_record:
        return False
    category = str(kb_record.get("category", "")).lower()
    return any(marker in category for marker in MONOGENIC_CATEGORY_MARKERS)


def _normalize_terms(terms: Iterable[Any]) -> List[str]:
    return [str(t).strip().lower() for t in (terms or []) if str(t).strip()]


def _severity(terms: List[str]) -> str:
    if any(t in LOF_CONSEQUENCES for t in terms):
        return "high"
    if any(t in MODERATE_CONSEQUENCES for t in terms):
        return "moderate"
    if any(t in LOW_CONSEQUENCES for t in terms):
        return "low"
    return "other"


def _allele_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> Optional[Tuple[str, int, str, str]]:
    try:
        return (
            str(chrom).removeprefix("chr"),
            int(pos),
            str(ref).upper(),
            str(alt).upper(),
        )
    except (TypeError, ValueError):
        return None


def _mutation_allele_key(mutation: Dict[str, Any]) -> Optional[Tuple[str, int, str, str]]:
    position = str(mutation.get("position") or "")
    if ":" not in position:
        return None
    chrom, _, pos = position.partition(":")
    return _allele_key(chrom, pos, mutation.get("ref"), mutation.get("alt"))


def _index_findings(findings: Iterable[Dict[str, Any]]) -> Dict[Tuple[str, int, str, str], List[Dict[str, Any]]]:
    index: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]] = {}
    for finding in findings or []:
        key = _allele_key(
            finding.get("chrom"), finding.get("pos"),
            finding.get("ref"), finding.get("alt"),
        )
        if key:
            index.setdefault(key, []).append(finding)
    return index


def _index_structures(structures: Iterable[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    for structure in structures or []:
        gene = structure.get("gene_symbol")
        if gene and not structure.get("error"):
            index.setdefault(str(gene).upper(), structure)
    return index


def _index_pgx(drug_results: Iterable[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Reduce drug results to the strongest evidence level seen per gene."""
    index: Dict[str, Dict[str, Any]] = {}
    for result in drug_results or []:
        gene = result.get("gene_symbol") or result.get("gene")
        if not gene:
            continue
        levels = []
        for drug in result.get("matched_drugs", []) or []:
            level = str(drug.get("evidence_level") or "").strip().lower()
            if level:
                levels.append(level)
        entry = index.setdefault(str(gene).upper(), {"levels": [], "drug_count": 0})
        entry["levels"].extend(levels)
        entry["drug_count"] += len(result.get("matched_drugs", []) or [])
    return index


def _pgx_strength(entry: Optional[Dict[str, Any]]) -> str:
    if not entry:
        return "none"
    levels = entry.get("levels", [])
    if any(level.startswith(STRONG_PGX_LEVELS) for level in levels):
        return "strong"
    if any(level.startswith(MODERATE_PGX_LEVELS) for level in levels):
        return "moderate"
    return "weak" if entry.get("drug_count") else "none"


def _clinvar_call(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Summarise ClinVar rows for one allele without resolving conflicts."""
    reportable = [r for r in records if r.get("reportable")]
    conflicting = [r for r in records if r.get("conflicting")]
    significances = sorted({
        str(r.get("significance")) for r in records if r.get("significance")
    })
    return {
        "matched": bool(records),
        "reportable": bool(reportable),
        "conflicting": bool(conflicting),
        "significances": significances,
        "synthetic": any(
            "synthetic" in str(r.get("source_name", "")).lower() for r in records
        ),
    }


class VariantTiering:
    """Assigns audit-ready tiers to annotated coding variants."""

    schema_version = "tiering-v1"

    def tier_variants(
        self,
        mutations: Iterable[Dict[str, Any]],
        clinvar_findings: Optional[Iterable[Dict[str, Any]]] = None,
        common_risk_findings: Optional[Iterable[Dict[str, Any]]] = None,
        structures: Optional[Iterable[Dict[str, Any]]] = None,
        drug_results: Optional[Iterable[Dict[str, Any]]] = None,
        clinvar_status: Optional[Dict[str, Any]] = None,
        common_risk_status: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        mutation_list = list(mutations or [])
        clinvar_index = _index_findings(clinvar_findings or [])
        risk_index = _index_findings(common_risk_findings or [])
        structure_index = _index_structures(structures or [])
        pgx_index = _index_pgx(drug_results or [])

        limitations = self._collect_limitations(clinvar_status, common_risk_status)

        tiered: List[Dict[str, Any]] = []
        for mutation in mutation_list:
            tiered.append(self._tier_one(
                mutation, clinvar_index, risk_index, structure_index, pgx_index,
            ))

        # Strongest evidence first, then by gene for stable presentation.
        tiered.sort(key=lambda item: (
            item["tier"],
            -item["evidence_points"],
            str(item.get("gene_symbol") or ""),
        ))

        counts = {tier: 0 for tier in (1, 2, 3, 4)}
        for item in tiered:
            counts[item["tier"]] += 1

        return {
            "status": "complete" if mutation_list else "no_coding_variants",
            "schema_version": self.schema_version,
            "evaluated_variants": len(mutation_list),
            "tier_counts": counts,
            "tier_labels": TIER_LABELS,
            "tier_descriptions": TIER_DESCRIPTIONS,
            "variants": tiered,
            "evidence_limitations": limitations,
            "method": (
                "Deterministic rule-based triage over VEP consequence, curated "
                "monogenic/common-risk diabetes gene context, local ClinVar "
                "assertions, AlphaFold structural context and ClinPGx evidence "
                "levels. Not ACMG classification; not a diagnosis."
            ),
        }

    @staticmethod
    def _collect_limitations(
        clinvar_status: Optional[Dict[str, Any]],
        common_risk_status: Optional[Dict[str, Any]],
    ) -> List[str]:
        limitations: List[str] = []
        for label, status in (
            ("ClinVar", clinvar_status or {}),
            ("Common-risk catalog", common_risk_status or {}),
        ):
            state = str(status.get("status") or "unknown")
            if state != "complete":
                limitations.append(
                    f"{label} evidence was not fully evaluated (status: {state}); "
                    "tiers relying on it may be under-called."
                )
            source_name = str(
                (status.get("source") or {}).get("name")
                if isinstance(status.get("source"), dict) else ""
            ).lower()
            if "synthetic" in source_name or "demo" in source_name:
                limitations.append(
                    f"{label} is backed by synthetic demonstration data, so any "
                    "tier it supports is illustrative only."
                )
        return limitations

    def _tier_one(
        self,
        mutation: Dict[str, Any],
        clinvar_index: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
        risk_index: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
        structure_index: Dict[str, Dict[str, Any]],
        pgx_index: Dict[str, Dict[str, Any]],
    ) -> Dict[str, Any]:
        gene = mutation.get("gene_symbol")
        gene_upper = str(gene).upper() if gene else None
        terms = _normalize_terms(mutation.get("consequence_terms"))
        severity = _severity(terms)

        kb_record = diabetes_kb.lookup(gene)
        monogenic = _is_monogenic_gene(kb_record)

        key = _mutation_allele_key(mutation)
        clinvar = _clinvar_call(clinvar_index.get(key, []) if key else [])
        risk_records = risk_index.get(key, []) if key else []
        structure = structure_index.get(gene_upper) if gene_upper else None
        pgx_strength = _pgx_strength(pgx_index.get(gene_upper) if gene_upper else None)

        structural_hit = bool(
            structure and (
                structure.get("in_domain")
                or structure.get("in_binding_site")
                or structure.get("in_active_site")
            )
        )
        confident_structure = bool(
            structure and str(structure.get("confidence_level", "")).lower()
            in {"very high", "high", "confident"}
        )

        reasons: List[str] = []
        points = 0
        tier = 4

        # --- Tier 1 rules -------------------------------------------------
        if clinvar["reportable"]:
            tier = 1
            points += 6
            reasons.append(
                "ClinVar reports a non-conflicting pathogenic/risk assertion: "
                + ", ".join(clinvar["significances"])
            )
        if monogenic and severity == "high":
            tier = min(tier, 1)
            points += 5
            reasons.append(
                f"Predicted loss-of-function ({', '.join(terms)}) in {gene}, an "
                f"established monogenic diabetes gene ({kb_record['category']})."
            )

        # --- Tier 2 rules -------------------------------------------------
        if risk_records:
            tier = min(tier, 2)
            points += 3
            traits = sorted({
                str(r.get("trait") or r.get("phenotype") or "diabetes-related")
                for r in risk_records
            })
            reasons.append(
                "Exact match in the curated common-risk association catalog: "
                + ", ".join(traits)
            )
        if monogenic and severity == "moderate":
            tier = min(tier, 2)
            points += 3
            reasons.append(
                f"Protein-altering change in monogenic diabetes gene {gene}."
            )
        if kb_record and severity in {"high", "moderate"} and structural_hit:
            tier = min(tier, 2)
            points += 2
            sites = [
                name for name, flag in (
                    ("annotated domain", structure.get("in_domain")),
                    ("binding site", structure.get("in_binding_site")),
                    ("active site", structure.get("in_active_site")),
                ) if flag
            ]
            reasons.append(
                f"Mutated residue falls in a functional region ({', '.join(sites)})."
            )
            if confident_structure:
                points += 1
                reasons.append(
                    f"AlphaFold confidence at this residue is "
                    f"{structure.get('confidence_level')} (pLDDT "
                    f"{structure.get('plddt_score')})."
                )
        if kb_record and pgx_strength == "strong":
            tier = min(tier, 2)
            points += 3
            reasons.append(
                "Gene carries high-confidence ClinPGx/PharmGKB evidence "
                "(level 1A/1B) for drug response."
            )

        # --- Tier 3 rules -------------------------------------------------
        if kb_record and tier == 4:
            tier = 3
            points += 1
            reasons.append(
                f"{gene} is a curated diabetes-associated gene "
                f"({kb_record['category']}), but this allele lacks strong "
                "corroborating evidence."
            )
        if pgx_strength == "moderate" and tier == 4:
            tier = 3
            reasons.append(
                "Gene has moderate pharmacogenomic evidence (level 2A/2B)."
            )
        if clinvar["matched"] and not clinvar["reportable"]:
            tier = min(tier, 3)
            reasons.append(
                "ClinVar contains conflicting or uncertain assertions for this "
                "allele; it is not presented as causal."
            )

        if tier == 4 and not reasons:
            reasons.append(
                "No curated diabetes gene context, ClinVar assertion, or "
                "pharmacogenomic evidence matched this allele."
            )

        return {
            "position": mutation.get("position"),
            "ref": mutation.get("ref"),
            "alt": mutation.get("alt"),
            "gene_symbol": gene,
            "consequence_terms": mutation.get("consequence_terms", []),
            "amino_acid_change": mutation.get("hgvsp") or mutation.get("amino_acid_change"),
            "tier": tier,
            "tier_label": TIER_LABELS[tier],
            "evidence_points": points,
            "reasons": reasons,
            "severity": severity,
            "monogenic_gene": monogenic,
            "diabetes_category": kb_record.get("category") if kb_record else None,
            "clinvar_matched": clinvar["matched"],
            "clinvar_reportable": clinvar["reportable"],
            "clinvar_conflicting": clinvar["conflicting"],
            "common_risk_matched": bool(risk_records),
            "structural_hit": structural_hit,
            "plddt_score": structure.get("plddt_score") if structure else None,
            "pgx_strength": pgx_strength,
            "uniprot_ids": mutation.get("uniprot_ids", []),
        }
