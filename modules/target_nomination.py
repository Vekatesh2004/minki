"""Integrative target nomination (Objectives 2, 3 and 4).

Fuses the independent evidence streams the pipeline produces — variant tiering,
PPI network topology, tissue expression, structural context and pharmacogenomic
annotations — into one ranked, auditable list of candidate therapeutic targets.

This is the integration step Objective 4 calls for: genomic variation plus
expression profiling plus network evidence feeding the screening stage.

Honesty constraints:
  * Ranking is a transparent weighted sum over evidence that was actually
    observed. Missing evidence contributes zero; it is never imputed.
  * Every nominated target carries the list of evidence that fired and the
    streams that were unavailable, so the score can be audited.
  * A high rank means "worth investigating", not "druggable" or "validated".
"""

from typing import Any, Dict, Iterable, List, Optional

from . import diabetes_kb

# Weight per evidence stream. Genomic evidence dominates because it is the only
# patient-specific input; the remaining streams are population/reference context.
WEIGHTS = {
    "variant_tier": 0.30,
    "network_hub": 0.20,
    "pathway": 0.10,
    "expression": 0.20,
    "structure": 0.10,
    "pharmacogenomic": 0.10,
}

TIER_CONTRIBUTION = {1: 1.0, 2: 0.7, 3: 0.35, 4: 0.0}
EXPRESSION_CONTRIBUTION = {"high": 1.0, "moderate": 0.6, "low": 0.2, "unknown": 0.0}
PGX_CONTRIBUTION = {"strong": 1.0, "moderate": 0.6, "weak": 0.3, "none": 0.0}


class TargetNominator:
    """Ranks candidate therapeutic targets from fused pipeline evidence."""

    schema_version = "target-nomination-v1"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        settings = (config or {}).get("target_nomination", {})
        self.max_targets = int(settings.get("max_targets", 25))
        self.min_score = float(settings.get("min_score", 0.05))

    def nominate(
        self,
        tiering: Optional[Dict[str, Any]] = None,
        ppi: Optional[Dict[str, Any]] = None,
        expression: Optional[Dict[str, Any]] = None,
        structures: Optional[Iterable[Dict[str, Any]]] = None,
        drug_results: Optional[Iterable[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        candidates: Dict[str, Dict[str, Any]] = {}

        missing_streams = [
            name for name, payload in (
                ("variant tiering", tiering),
                ("PPI network", ppi),
                ("gene expression", expression),
            )
            if not payload or payload.get("status") != "complete"
        ]

        self._add_tier_evidence(candidates, tiering)
        self._add_network_evidence(candidates, ppi)
        self._add_pathway_evidence(candidates, ppi)
        self._add_expression_evidence(candidates, expression)
        self._add_structure_evidence(candidates, structures)
        self._add_pgx_evidence(candidates, drug_results)

        nominated: List[Dict[str, Any]] = []
        for gene, record in candidates.items():
            score = sum(
                WEIGHTS[stream] * value
                for stream, value in record["contributions"].items()
            )
            if score < self.min_score:
                continue
            kb_record = diabetes_kb.lookup(gene)
            nominated.append({
                "gene_symbol": gene,
                "target_score": round(score, 4),
                "evidence": record["evidence"],
                "contributions": {
                    stream: round(WEIGHTS[stream] * value, 4)
                    for stream, value in record["contributions"].items()
                },
                "evidence_stream_count": len(record["contributions"]),
                "patient_variant": record.get("patient_variant", False),
                "best_tier": record.get("best_tier"),
                "diabetes_category": kb_record.get("category") if kb_record else None,
                "known_drug_target": record.get("known_drug_target", False),
                "rationale": self._rationale(gene, record, kb_record),
            })

        nominated.sort(key=lambda item: (
            -item["target_score"], -item["evidence_stream_count"], item["gene_symbol"],
        ))
        nominated = nominated[: self.max_targets]

        return {
            "status": "complete" if nominated else "no_candidates",
            "schema_version": self.schema_version,
            "candidate_count": len(nominated),
            "candidates": nominated,
            "weights": WEIGHTS,
            "missing_evidence_streams": missing_streams,
            "method": (
                "Transparent weighted fusion of variant tier, PPI hub score, "
                "pathway membership, GTEx tissue relevance, structural context "
                "and pharmacogenomic evidence level. Absent evidence scores zero "
                "and is never imputed."
            ),
            "interpretation": (
                "A high score marks a gene as worth investigating further. It is "
                "not evidence of druggability, efficacy, or clinical validity."
            ),
        }

    # -- Evidence collectors ----------------------------------------------
    @staticmethod
    def _entry(candidates: Dict[str, Dict[str, Any]], gene: str) -> Dict[str, Any]:
        return candidates.setdefault(str(gene).upper(), {
            "contributions": {}, "evidence": [],
        })

    def _add_tier_evidence(
        self, candidates: Dict[str, Dict[str, Any]], tiering: Optional[Dict[str, Any]],
    ) -> None:
        if not tiering or tiering.get("status") != "complete":
            return
        for variant in tiering.get("variants", []) or []:
            gene = variant.get("gene_symbol")
            if not gene:
                continue
            tier = variant.get("tier", 4)
            contribution = TIER_CONTRIBUTION.get(tier, 0.0)
            if contribution <= 0:
                continue
            entry = self._entry(candidates, gene)
            # Keep the strongest tier seen for this gene.
            previous = entry["contributions"].get("variant_tier", 0.0)
            if contribution >= previous:
                entry["contributions"]["variant_tier"] = contribution
                entry["best_tier"] = tier
            entry["patient_variant"] = True
            label = variant.get("tier_label", f"Tier {tier}")
            change = variant.get("amino_acid_change") or variant.get("position")
            note = f"{label} variant in this sample ({change})"
            if note not in entry["evidence"]:
                entry["evidence"].append(note)

    def _add_network_evidence(
        self, candidates: Dict[str, Dict[str, Any]], ppi: Optional[Dict[str, Any]],
    ) -> None:
        if not ppi or ppi.get("status") != "complete":
            return
        hubs = ppi.get("hubs", []) or []
        if not hubs:
            return
        top_score = max((h.get("hub_score") or 0.0) for h in hubs) or 1.0
        for hub in hubs:
            gene = hub.get("gene_symbol")
            if not gene:
                continue
            hub_score = hub.get("hub_score") or 0.0
            if hub_score <= 0:
                continue
            entry = self._entry(candidates, gene)
            entry["contributions"]["network_hub"] = min(1.0, hub_score / top_score)
            entry["evidence"].append(
                f"PPI network hub (degree {hub.get('degree')}, hub score "
                f"{hub_score})"
            )
            if not hub.get("is_seed"):
                entry["evidence"].append(
                    "Nominated from network topology despite carrying no variant "
                    "in this sample"
                )

    def _add_pathway_evidence(
        self, candidates: Dict[str, Dict[str, Any]], ppi: Optional[Dict[str, Any]],
    ) -> None:
        if not ppi or ppi.get("status") != "complete":
            return
        terms = ppi.get("enrichment", []) or []
        pathway_terms = [
            t for t in terms if t.get("category") in {"KEGG", "RCTM", "WikiPathways"}
        ]
        if not pathway_terms:
            return
        # Count how many significant pathways each gene participates in.
        counts: Dict[str, int] = {}
        examples: Dict[str, str] = {}
        for term in pathway_terms[:20]:
            for gene in term.get("genes", []) or []:
                key = str(gene).upper()
                counts[key] = counts.get(key, 0) + 1
                examples.setdefault(key, str(term.get("description") or ""))
        if not counts:
            return
        top = max(counts.values())
        for gene, count in counts.items():
            entry = self._entry(candidates, gene)
            entry["contributions"]["pathway"] = min(1.0, count / top)
            entry["evidence"].append(
                f"Member of {count} enriched pathway(s), e.g. {examples[gene]}"
            )

    def _add_expression_evidence(
        self, candidates: Dict[str, Dict[str, Any]], expression: Optional[Dict[str, Any]],
    ) -> None:
        if not expression or expression.get("status") != "complete":
            return
        for profile in expression.get("profiles", []) or []:
            gene = profile.get("gene_symbol")
            if not gene:
                continue
            relevance = profile.get("diabetes_tissue_relevance", {}) or {}
            level = str(relevance.get("level", "unknown"))
            contribution = EXPRESSION_CONTRIBUTION.get(level, 0.0)
            if contribution <= 0:
                continue
            entry = self._entry(candidates, gene)
            entry["contributions"]["expression"] = contribution
            entry["evidence"].append(
                f"GTEx tissue relevance: {level} — {relevance.get('rationale')}"
            )

    def _add_structure_evidence(
        self,
        candidates: Dict[str, Dict[str, Any]],
        structures: Optional[Iterable[Dict[str, Any]]],
    ) -> None:
        for structure in structures or []:
            gene = structure.get("gene_symbol")
            if not gene or structure.get("error"):
                continue
            functional_site = bool(
                structure.get("in_binding_site") or structure.get("in_active_site")
            )
            in_domain = bool(structure.get("in_domain"))
            if not (functional_site or in_domain):
                continue
            entry = self._entry(candidates, gene)
            entry["contributions"]["structure"] = 1.0 if functional_site else 0.5
            site = "binding/active site" if functional_site else "annotated domain"
            entry["evidence"].append(
                f"Mutated residue lies in a {site}; AlphaFold confidence "
                f"{structure.get('confidence_level') or 'unknown'}"
            )

    def _add_pgx_evidence(
        self,
        candidates: Dict[str, Dict[str, Any]],
        drug_results: Optional[Iterable[Dict[str, Any]]],
    ) -> None:
        for result in drug_results or []:
            gene = result.get("gene_symbol") or result.get("gene")
            matched = result.get("matched_drugs", []) or []
            if not gene or not matched:
                continue
            levels = [
                str(drug.get("evidence_level") or "").lower() for drug in matched
            ]
            if any(level.startswith(("1a", "1b")) for level in levels):
                strength = "strong"
            elif any(level.startswith(("2a", "2b")) for level in levels):
                strength = "moderate"
            else:
                strength = "weak"
            entry = self._entry(candidates, gene)
            entry["contributions"]["pharmacogenomic"] = PGX_CONTRIBUTION[strength]
            entry["known_drug_target"] = True
            entry["evidence"].append(
                f"{len(matched)} pharmacogenomic drug annotation(s), "
                f"{strength} evidence"
            )

    @staticmethod
    def _rationale(
        gene: str, record: Dict[str, Any], kb_record: Optional[Dict[str, Any]],
    ) -> str:
        streams = record["contributions"]
        parts: List[str] = []
        if "variant_tier" in streams:
            parts.append(f"carries a Tier {record.get('best_tier')} variant")
        if "network_hub" in streams:
            parts.append("is a hub in the interaction network")
        if "pathway" in streams:
            parts.append("sits in enriched diabetes-relevant pathways")
        if "expression" in streams:
            parts.append("is expressed in metabolically relevant tissues")
        if "structure" in streams:
            parts.append("has a structurally consequential mutation site")
        if "pharmacogenomic" in streams:
            parts.append("has existing pharmacogenomic drug evidence")
        summary = f"{gene} " + ", ".join(parts) if parts else f"{gene} has limited evidence"
        if kb_record:
            summary += f". Curated context: {kb_record['category']}."
        return summary
