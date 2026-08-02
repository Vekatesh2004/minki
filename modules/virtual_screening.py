"""Ligand-based virtual screening for antidiabetic candidates (Objective 4, part 2).

For each candidate protein target nominated by the pipeline, this module
retrieves measured bioactivity data from ChEMBL, ranks the compounds by
potency and drug-likeness, and reports which are already approved drugs versus
novel chemical starting points.

IMPORTANT — what this is and is not:
  * This is LIGAND-BASED screening over experimentally measured affinities
    (IC50/Ki/EC50) curated in ChEMBL. It is NOT structure-based molecular
    docking: no binding pose is generated and no docking score is computed.
    No AutoDock/Vina engine is invoked, and none is required.
  * Potency values are real experimental measurements from the literature, not
    predictions. Assay conditions vary between records, so cross-compound
    comparison is indicative only.
  * Drug-likeness uses ChEMBL's precomputed physicochemical descriptors
    (Lipinski violations, QED, TPSA). These are heuristics for developability,
    not efficacy or safety predictions.
  * Nothing here constitutes a therapeutic recommendation.
"""

import asyncio
import logging
import math
from typing import Any, Dict, Iterable, List, Optional

import aiohttp

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

# Activity types that represent a direct binding/functional potency measurement.
POTENCY_TYPES = ("IC50", "Ki", "EC50", "Kd")

# max_phase in ChEMBL: 4 = approved, 1-3 = clinical phases, 0/None = preclinical.
PHASE_LABELS = {
    4.0: "Approved drug",
    3.0: "Phase 3 clinical",
    2.0: "Phase 2 clinical",
    1.0: "Phase 1 clinical",
    0.5: "Early clinical",
}


def _phase_label(max_phase: Any) -> str:
    try:
        value = float(max_phase)
    except (TypeError, ValueError):
        return "Preclinical / research compound"
    return PHASE_LABELS.get(value, "Preclinical / research compound")


def _is_approved(max_phase: Any) -> bool:
    """ChEMBL returns max_phase as a string (e.g. '4.0'), so compare numerically.

    Only the molecule record's own max_phase is trusted. The activity
    endpoint's ``molecule_max_phase`` filter is silently ignored by ChEMBL, so
    it must never be used to infer approval status.
    """
    try:
        return float(max_phase) >= 4.0
    except (TypeError, ValueError):
        return False


def _pactivity(value_nm: Optional[float]) -> Optional[float]:
    """Convert nM potency to pActivity (-log10 M). Higher is more potent."""
    if not value_nm or value_nm <= 0:
        return None
    try:
        return round(-math.log10(value_nm * 1e-9), 2)
    except (ValueError, OverflowError):
        return None


class VirtualScreeningEngine:
    """Ranks ChEMBL bioactive compounds against pipeline-nominated targets."""

    schema_version = "screening-v1"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        settings = (config or {}).get("virtual_screening", {})
        self.base_url = settings.get("base_url", CHEMBL_BASE)
        self.max_targets = int(settings.get("max_targets", 6))
        self.max_compounds_per_target = int(settings.get("max_compounds_per_target", 12))
        self.activity_fetch_limit = int(settings.get("activity_fetch_limit", 120))
        self.timeout = float(settings.get("timeout", 60.0))
        self.concurrency = int(settings.get("concurrency", 3))
        self.enabled = bool(settings.get("enabled", True))

    async def _get_json(
        self, session: aiohttp.ClientSession, path: str, params: Dict[str, Any],
    ) -> Any:
        params = {**params, "format": "json"}
        async with session.get(f"{self.base_url}{path}", params=params) as response:
            if response.status != 200:
                body = await response.text()
                raise RuntimeError(f"ChEMBL {path} HTTP {response.status}: {body[:180]}")
            return await response.json(content_type=None)

    async def _resolve_target(
        self, session: aiohttp.ClientSession, gene: str,
    ) -> Optional[Dict[str, Any]]:
        """Find the human single-protein ChEMBL target for a gene symbol."""
        payload = await self._get_json(
            session, "/target/search", {"q": gene, "limit": 10},
        )
        candidates = payload.get("targets", []) or []
        human_single = [
            t for t in candidates
            if t.get("organism") == "Homo sapiens"
            and t.get("target_type") == "SINGLE PROTEIN"
        ]
        pool = human_single or [
            t for t in candidates if t.get("organism") == "Homo sapiens"
        ]
        if not pool:
            return None
        # Prefer a target whose synonyms explicitly list the gene symbol.
        for target in pool:
            synonyms = {
                str(c.get("component_synonym", "")).upper()
                for component in target.get("target_components", []) or []
                for c in component.get("target_component_synonyms", []) or []
            }
            if gene.upper() in synonyms:
                return target
        return pool[0]

    async def _fetch_activities(
        self, session: aiohttp.ClientSession, target_chembl_id: str,
    ) -> List[Dict[str, Any]]:
        """Pull measured potencies, best (lowest nM) first."""
        collected: List[Dict[str, Any]] = []
        for activity_type in POTENCY_TYPES:
            try:
                payload = await self._get_json(session, "/activity", {
                    "target_chembl_id": target_chembl_id,
                    "standard_type": activity_type,
                    "standard_units": "nM",
                    "limit": max(20, self.activity_fetch_limit // len(POTENCY_TYPES)),
                    "order_by": "standard_value",
                })
            except Exception as exc:
                logger.debug("ChEMBL activity fetch failed (%s): %s", activity_type, exc)
                continue
            collected.extend(payload.get("activities", []) or [])
        return collected

    async def _fetch_mechanisms(
        self, session: aiohttp.ClientSession, target_chembl_id: str,
    ) -> List[Dict[str, Any]]:
        """Retrieve established drug mechanisms of action for a target.

        The most potent ligands for a target are usually unapproved research
        compounds, so approved drugs rarely survive a potency-only shortlist.
        This endpoint is the authoritative source for clinically advanced drugs
        acting on the target, which is what repurposing needs.
        """
        try:
            payload = await self._get_json(session, "/mechanism", {
                "target_chembl_id": target_chembl_id, "limit": 50,
            })
        except Exception as exc:
            logger.debug("ChEMBL mechanism fetch failed: %s", exc)
            return []

        mechanisms: List[Dict[str, Any]] = []
        for row in payload.get("mechanisms", []) or []:
            molecule_id = row.get("molecule_chembl_id")
            if not molecule_id:
                continue
            mechanisms.append({
                "molecule_chembl_id": molecule_id,
                "action_type": row.get("action_type"),
                "mechanism_of_action": row.get("mechanism_of_action"),
                "max_phase": row.get("max_phase"),
                "is_approved": _is_approved(row.get("max_phase")),
            })
        return mechanisms

    async def _fetch_molecule(
        self, session: aiohttp.ClientSession, molecule_id: str,
    ) -> Optional[Dict[str, Any]]:
        try:
            return await self._get_json(session, f"/molecule/{molecule_id}", {})
        except Exception:
            return None

    @staticmethod
    def _best_per_molecule(activities: Iterable[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Keep the most potent measurement per compound."""
        best: Dict[str, Dict[str, Any]] = {}
        for activity in activities:
            molecule_id = activity.get("molecule_chembl_id")
            if not molecule_id:
                continue
            try:
                value = float(activity.get("standard_value"))
            except (TypeError, ValueError):
                continue
            if value <= 0:
                continue
            relation = str(activity.get("standard_relation") or "=")
            # Censored values ('>') understate potency; skip them.
            if relation.startswith(">"):
                continue
            current = best.get(molecule_id)
            if current is None or value < current["value_nm"]:
                best[molecule_id] = {
                    "molecule_chembl_id": molecule_id,
                    "value_nm": value,
                    "activity_type": activity.get("standard_type"),
                    "relation": relation,
                    "assay_chembl_id": activity.get("assay_chembl_id"),
                    "assay_description": activity.get("assay_description"),
                    "smiles": activity.get("canonical_smiles"),
                    "pchembl_value": activity.get("pchembl_value"),
                }
        return best

    @staticmethod
    def _drug_likeness(properties: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Summarise ChEMBL's curated physicochemical descriptors."""
        if not properties:
            return {"available": False}

        def _num(key: str) -> Optional[float]:
            try:
                return float(properties.get(key))
            except (TypeError, ValueError):
                return None

        violations = properties.get("num_ro5_violations")
        try:
            violations = int(violations)
        except (TypeError, ValueError):
            violations = None
        qed = _num("qed_weighted")

        if violations is None:
            verdict = "unknown"
        elif violations == 0:
            verdict = "passes Lipinski rule of five"
        elif violations == 1:
            verdict = "1 Lipinski violation (often tolerated)"
        else:
            verdict = f"{violations} Lipinski violations"

        return {
            "available": True,
            "molecular_weight": _num("full_mwt"),
            "alogp": _num("alogp"),
            "tpsa": _num("psa"),
            "hbd": properties.get("hbd"),
            "hba": properties.get("hba"),
            "rotatable_bonds": properties.get("rtb"),
            "aromatic_rings": properties.get("aromatic_rings"),
            "lipinski_violations": violations,
            "qed_weighted": qed,
            "verdict": verdict,
            "oral_druglike": bool(violations is not None and violations <= 1),
        }

    @staticmethod
    def _candidate_score(
        pactivity: Optional[float], drug_likeness: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Composite prioritisation score: potency plus developability.

        Deliberately transparent and reproducible. Potency dominates; QED and
        Lipinski compliance modulate. This is a triage ranking, not a
        predicted binding affinity.
        """
        components = []
        score = 0.0
        if pactivity is not None:
            # pActivity 5 (10 uM) -> 0.0 ; 10 (0.1 nM) -> 1.0
            potency_component = max(0.0, min(1.0, (pactivity - 5.0) / 5.0))
            score += potency_component * 0.65
            components.append({
                "component": "measured potency",
                "value": pactivity,
                "weight": 0.65,
            })
        qed = drug_likeness.get("qed_weighted")
        if qed is not None:
            score += max(0.0, min(1.0, qed)) * 0.25
            components.append({
                "component": "QED drug-likeness", "value": qed, "weight": 0.25,
            })
        violations = drug_likeness.get("lipinski_violations")
        if violations is not None:
            lipinski_component = 1.0 if violations == 0 else (0.5 if violations == 1 else 0.0)
            score += lipinski_component * 0.10
            components.append({
                "component": "Lipinski compliance",
                "value": violations,
                "weight": 0.10,
            })
        return {"score": round(score, 4), "components": components}

    async def _screen_target(
        self,
        session: aiohttp.ClientSession,
        gene: str,
        rationale: Optional[str],
        semaphore: asyncio.Semaphore,
    ) -> Dict[str, Any]:
        async with semaphore:
            try:
                target = await self._resolve_target(session, gene)
            except Exception as exc:
                return {
                    "gene_symbol": gene, "status": "error",
                    "note": f"ChEMBL target lookup failed: {exc}",
                }
            if not target:
                return {
                    "gene_symbol": gene, "status": "no_target",
                    "note": "No human ChEMBL target found for this gene.",
                }

            target_id = target.get("target_chembl_id")
            activities, mechanisms = await asyncio.gather(
                self._fetch_activities(session, target_id),
                self._fetch_mechanisms(session, target_id),
                return_exceptions=True,
            )
            if isinstance(activities, Exception):
                activities = []
            if isinstance(mechanisms, Exception):
                mechanisms = []
            best = self._best_per_molecule(activities)
            mechanism_by_molecule = {
                item["molecule_chembl_id"]: item for item in mechanisms
            }
            if not best and not mechanisms:
                return {
                    "gene_symbol": gene,
                    "status": "no_ligands",
                    "target_chembl_id": target_id,
                    "target_name": target.get("pref_name"),
                    "note": (
                        "Target exists in ChEMBL but has no uncensored "
                        "IC50/Ki/EC50/Kd measurements in nM and no recorded "
                        "drug mechanism of action."
                    ),
                }

            ranked = sorted(best.values(), key=lambda item: item["value_nm"])
            shortlist = list(ranked[: self.max_compounds_per_target])
            shortlist_ids = {item["molecule_chembl_id"] for item in shortlist}
            # Always include clinically advanced drugs for this target, even
            # without a potency record, so repurposing options are never hidden
            # behind more potent preclinical chemistry.
            for molecule_id, mechanism in mechanism_by_molecule.items():
                if not mechanism["is_approved"] or molecule_id in shortlist_ids:
                    continue
                shortlist.append({
                    "molecule_chembl_id": molecule_id,
                    "value_nm": best.get(molecule_id, {}).get("value_nm"),
                    "activity_type": best.get(molecule_id, {}).get("activity_type"),
                    "relation": "=",
                    "assay_chembl_id": None,
                    "assay_description": None,
                    "smiles": None,
                    "pchembl_value": None,
                })
                shortlist_ids.add(molecule_id)

            molecules = await asyncio.gather(*[
                self._fetch_molecule(session, item["molecule_chembl_id"])
                for item in shortlist
            ], return_exceptions=True)

            compounds: List[Dict[str, Any]] = []
            for item, molecule in zip(shortlist, molecules):
                molecule = None if isinstance(molecule, Exception) else molecule
                properties = (molecule or {}).get("molecule_properties")
                drug_likeness = self._drug_likeness(properties)
                pactivity = (
                    float(item["pchembl_value"])
                    if item.get("pchembl_value") not in (None, "")
                    else _pactivity(item.get("value_nm"))
                )
                scoring = self._candidate_score(pactivity, drug_likeness)
                structures = (molecule or {}).get("molecule_structures") or {}
                mechanism = mechanism_by_molecule.get(item["molecule_chembl_id"], {})
                # Prefer the mechanism record's max_phase: it is populated for
                # clinically advanced drugs whose molecule record may omit it.
                max_phase = (molecule or {}).get("max_phase")
                if max_phase in (None, "") and mechanism:
                    max_phase = mechanism.get("max_phase")
                compounds.append({
                    "molecule_chembl_id": item["molecule_chembl_id"],
                    "preferred_name": (molecule or {}).get("pref_name"),
                    "development_stage": _phase_label(max_phase),
                    "is_approved_drug": _is_approved(max_phase),
                    "action_type": mechanism.get("action_type"),
                    "mechanism_of_action": mechanism.get("mechanism_of_action"),
                    "activity_type": item.get("activity_type"),
                    "potency_nm": item.get("value_nm"),
                    "pactivity": pactivity,
                    "assay_description": (item.get("assay_description") or "")[:220],
                    "assay_chembl_id": item.get("assay_chembl_id"),
                    "smiles": structures.get("canonical_smiles") or item.get("smiles"),
                    "drug_likeness": drug_likeness,
                    "priority_score": scoring["score"],
                    "score_components": scoring["components"],
                    "chembl_url": (
                        "https://www.ebi.ac.uk/chembl/compound_report_card/"
                        f"{item['molecule_chembl_id']}/"
                    ),
                })

            # Approved drugs first (repurposing is immediately actionable),
            # then by composite priority score.
            compounds.sort(key=lambda c: (
                not c["is_approved_drug"], -c["priority_score"],
            ))
            approved = [c for c in compounds if c["is_approved_drug"]]
            novel = [
                c for c in compounds
                if not c["is_approved_drug"] and c["drug_likeness"].get("oral_druglike")
            ]

            return {
                "gene_symbol": gene,
                "status": "complete",
                "target_chembl_id": target_id,
                "target_name": target.get("pref_name"),
                "target_type": target.get("target_type"),
                "target_rationale": rationale,
                "measured_compound_count": len(best),
                "compounds": compounds,
                "approved_drug_count": len(approved),
                "novel_druglike_count": len(novel),
                "repurposing_candidates": [
                    {
                        "molecule_chembl_id": c["molecule_chembl_id"],
                        "preferred_name": c["preferred_name"],
                        "potency_nm": c["potency_nm"],
                        "activity_type": c["activity_type"],
                        "action_type": c.get("action_type"),
                        "mechanism_of_action": c.get("mechanism_of_action"),
                    }
                    for c in approved[:5]
                ],
                "novel_starting_points": [
                    {
                        "molecule_chembl_id": c["molecule_chembl_id"],
                        "potency_nm": c["potency_nm"],
                        "qed_weighted": c["drug_likeness"].get("qed_weighted"),
                    }
                    for c in novel[:5]
                ],
                "target_chembl_url": (
                    f"https://www.ebi.ac.uk/chembl/target_report_card/{target_id}/"
                ),
            }

    async def screen_targets(
        self, targets: Iterable[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Screen a prioritised target list.

        Each target is a dict with at least ``gene_symbol``, optionally
        ``rationale`` explaining why the pipeline nominated it.
        """
        deduped: List[Dict[str, Any]] = []
        seen = set()
        for target in targets or []:
            gene = target.get("gene_symbol") if isinstance(target, dict) else target
            if not gene:
                continue
            symbol = str(gene).strip().upper()
            if symbol in seen:
                continue
            seen.add(symbol)
            deduped.append({
                "gene_symbol": symbol,
                "rationale": target.get("rationale") if isinstance(target, dict) else None,
            })
        deduped = deduped[: self.max_targets]

        if not self.enabled:
            return self._empty("disabled", "Virtual screening is disabled in config.")
        if not deduped:
            return self._empty(
                "no_targets",
                "No candidate targets were nominated by the upstream analysis.",
            )

        try:
            timeout = aiohttp.ClientTimeout(total=self.timeout)
            semaphore = asyncio.Semaphore(max(1, self.concurrency))
            async with aiohttp.ClientSession(timeout=timeout) as session:
                results = await asyncio.gather(*[
                    self._screen_target(
                        session, target["gene_symbol"], target["rationale"], semaphore,
                    )
                    for target in deduped
                ], return_exceptions=True)
        except asyncio.TimeoutError:
            return self._empty(
                "timeout", f"ChEMBL did not respond within {self.timeout:.0f}s.",
            )
        except Exception as exc:
            logger.error("Virtual screening failed: %s", exc)
            return self._empty("error", f"Virtual screening error: {exc}")

        screened: List[Dict[str, Any]] = []
        warnings: List[str] = []
        for target, result in zip(deduped, results):
            if isinstance(result, Exception):
                warnings.append(f"{target['gene_symbol']}: {result}")
                continue
            if result.get("status") == "complete":
                screened.append(result)
            else:
                warnings.append(
                    f"{target['gene_symbol']}: {result.get('note', result.get('status'))}"
                )

        total_compounds = sum(len(t["compounds"]) for t in screened)
        return {
            "status": "complete" if screened else "no_data",
            "schema_version": self.schema_version,
            "method": "ligand-based (ChEMBL measured bioactivity)",
            "docking_performed": False,
            "screened_targets": screened,
            "target_count": len(screened),
            "total_compounds": total_compounds,
            "approved_drug_count": sum(t["approved_drug_count"] for t in screened),
            "novel_druglike_count": sum(t["novel_druglike_count"] for t in screened),
            "warnings": warnings,
            "source": {
                "name": "ChEMBL",
                "url": "https://www.ebi.ac.uk/chembl/",
                "note": (
                    "Potency values are experimentally measured and literature-"
                    "curated. Assay conditions differ between records."
                ),
            },
            "interpretation": (
                "Ligand-based prioritisation over measured affinities. No binding "
                "pose or docking score is computed, so these are chemical "
                "starting points for follow-up, not validated drug candidates."
            ),
        }

    def _empty(self, status: str, warning: str) -> Dict[str, Any]:
        return {
            "status": status,
            "schema_version": self.schema_version,
            "method": "ligand-based (ChEMBL measured bioactivity)",
            "docking_performed": False,
            "screened_targets": [],
            "target_count": 0,
            "total_compounds": 0,
            "approved_drug_count": 0,
            "novel_druglike_count": 0,
            "warnings": [warning],
            "source": {"name": "ChEMBL", "url": "https://www.ebi.ac.uk/chembl/"},
            "interpretation": None,
        }
