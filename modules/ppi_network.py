"""Protein-protein interaction network construction and analysis (Objective 3).

Builds a real PPI network for the patient's mutated diabetes-relevant genes
using the public STRING database, then derives graph-theoretic centrality to
nominate candidate functional biomarkers, and retrieves pathway/GO functional
enrichment for the same gene set.

Honesty constraints:
  * Hub status is a network-topology observation, not proof that a protein is a
    clinically valid biomarker. Output is labelled as a research hypothesis.
  * STRING combined scores mix predicted and experimental channels. The
    experimental sub-score is reported separately so predicted-only edges are
    never presented as measured interactions.
  * Enrichment FDR values come from STRING; they are not recomputed here.
  * If STRING is unreachable, the module reports an explicit failure status
    rather than fabricating a network.
"""

import asyncio
import logging
from typing import Any, Dict, Iterable, List, Optional, Set

import aiohttp

logger = logging.getLogger(__name__)

try:
    import networkx as nx
    _HAVE_NETWORKX = True
except ImportError:  # graph metrics degrade gracefully
    _HAVE_NETWORKX = False

STRING_BASE = "https://string-db.org/api"
HUMAN_TAXON = 9606

# STRING combined-score bands (scores are returned 0-1 by the JSON API).
CONFIDENCE_BANDS = (
    (0.900, "highest"),
    (0.700, "high"),
    (0.400, "medium"),
    (0.150, "low"),
)

# Categories STRING returns for functional enrichment, mapped to display names.
ENRICHMENT_CATEGORIES = {
    "KEGG": "KEGG pathway",
    "RCTM": "Reactome pathway",
    "Process": "GO biological process",
    "Function": "GO molecular function",
    "Component": "GO cellular component",
    "WikiPathways": "WikiPathways",
}


def _confidence_band(score: float) -> str:
    for threshold, label in CONFIDENCE_BANDS:
        if score >= threshold:
            return label
    return "very low"


class PPINetworkAnalyzer:
    """Constructs and analyses a STRING-derived PPI network."""

    schema_version = "ppi-v1"

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        settings = (config or {}).get("ppi_network", {})
        self.base_url = settings.get("base_url", STRING_BASE)
        self.required_score = int(settings.get("required_score", 400))
        self.max_seed_genes = int(settings.get("max_seed_genes", 40))
        self.expansion_limit = int(settings.get("expansion_limit", 20))
        self.timeout = float(settings.get("timeout", 45.0))
        self.enabled = bool(settings.get("enabled", True))

    # -- HTTP -------------------------------------------------------------
    async def _get_json(
        self, session: aiohttp.ClientSession, endpoint: str, params: Dict[str, Any],
    ) -> Any:
        url = f"{self.base_url}/json/{endpoint}"
        async with session.get(url, params=params) as response:
            if response.status != 200:
                body = await response.text()
                raise RuntimeError(
                    f"STRING {endpoint} returned HTTP {response.status}: {body[:200]}"
                )
            return await response.json(content_type=None)

    @staticmethod
    def _identifier_payload(genes: Iterable[str]) -> str:
        # STRING expects carriage-return separated identifiers.
        return "\r".join(genes)

    # -- Public API -------------------------------------------------------
    async def analyze(
        self,
        genes: Iterable[str],
        expand: bool = True,
    ) -> Dict[str, Any]:
        seed_genes = self._clean_genes(genes)
        if not self.enabled:
            return self._empty("disabled", seed_genes, "PPI analysis is disabled in config.")
        if len(seed_genes) < 2:
            return self._empty(
                "insufficient_genes", seed_genes,
                "At least two mutated genes are required to build an interaction network.",
            )

        try:
            timeout = aiohttp.ClientTimeout(total=self.timeout)
            async with aiohttp.ClientSession(timeout=timeout) as session:
                params: Dict[str, Any] = {
                    "identifiers": self._identifier_payload(seed_genes),
                    "species": HUMAN_TAXON,
                    "required_score": self.required_score,
                    "caller_identity": "minki-pgx-pipeline",
                }
                if expand and self.expansion_limit > 0:
                    # Pull in first-shell partners so hubs outside the mutated
                    # set (candidate biomarkers) can be discovered.
                    params["add_nodes"] = self.expansion_limit

                edges_raw, enrichment_raw = await asyncio.gather(
                    self._get_json(session, "network", params),
                    self._get_json(session, "enrichment", {
                        "identifiers": self._identifier_payload(seed_genes),
                        "species": HUMAN_TAXON,
                        "caller_identity": "minki-pgx-pipeline",
                    }),
                    return_exceptions=True,
                )

            if isinstance(edges_raw, Exception):
                raise edges_raw

            edges = self._parse_edges(edges_raw)
            if not edges:
                return self._empty(
                    "no_interactions", seed_genes,
                    "STRING returned no interactions above the confidence threshold "
                    f"({self.required_score}/1000) for these genes.",
                )

            enrichment = (
                [] if isinstance(enrichment_raw, Exception)
                else self._parse_enrichment(enrichment_raw)
            )
            enrichment_warning = (
                f"Functional enrichment unavailable: {enrichment_raw}"
                if isinstance(enrichment_raw, Exception) else None
            )

            analysis = self._analyze_graph(edges, set(seed_genes))
            warnings: List[str] = []
            if enrichment_warning:
                warnings.append(enrichment_warning)
            if not _HAVE_NETWORKX:
                warnings.append(
                    "networkx is not installed; centrality is limited to degree "
                    "counts computed directly from the edge list."
                )

            return {
                "status": "complete",
                "schema_version": self.schema_version,
                "seed_genes": seed_genes,
                "required_score": self.required_score,
                "expanded": bool(expand and self.expansion_limit > 0),
                "expansion_limit": self.expansion_limit if expand else 0,
                "edges": edges,
                "enrichment": enrichment,
                "warnings": warnings,
                "source": {
                    "name": "STRING",
                    "url": "https://string-db.org",
                    "species": "Homo sapiens (9606)",
                    "note": (
                        "Combined scores integrate experimental and predicted "
                        "evidence channels. Experimental sub-scores are reported "
                        "separately."
                    ),
                },
                "interpretation": (
                    "Network hubs are research hypotheses for functional "
                    "biomarkers, derived from interaction topology. They are not "
                    "validated biomarkers and require experimental confirmation."
                ),
                **analysis,
            }

        except asyncio.TimeoutError:
            return self._empty(
                "timeout", seed_genes,
                f"STRING did not respond within {self.timeout:.0f}s.",
            )
        except Exception as exc:
            logger.error("PPI network construction failed: %s", exc)
            return self._empty("error", seed_genes, f"PPI network error: {exc}")

    # -- Helpers ----------------------------------------------------------
    def _clean_genes(self, genes: Iterable[str]) -> List[str]:
        seen: List[str] = []
        for gene in genes or []:
            if not gene:
                continue
            symbol = str(gene).strip().upper()
            if symbol and symbol not in seen:
                seen.append(symbol)
        return seen[: self.max_seed_genes]

    def _empty(self, status: str, seed_genes: List[str], warning: str) -> Dict[str, Any]:
        return {
            "status": status,
            "schema_version": self.schema_version,
            "seed_genes": seed_genes,
            "required_score": self.required_score,
            "expanded": False,
            "expansion_limit": 0,
            "edges": [],
            "nodes": [],
            "hubs": [],
            "modules": [],
            "enrichment": [],
            "node_count": 0,
            "edge_count": 0,
            "network_density": 0.0,
            "warnings": [warning],
            "source": {"name": "STRING", "url": "https://string-db.org"},
            "interpretation": None,
        }

    @staticmethod
    def _parse_edges(payload: Any) -> List[Dict[str, Any]]:
        edges: List[Dict[str, Any]] = []
        deduped: Set[frozenset] = set()
        for row in payload or []:
            gene_a = str(row.get("preferredName_A") or "").upper()
            gene_b = str(row.get("preferredName_B") or "").upper()
            if not gene_a or not gene_b or gene_a == gene_b:
                continue
            pair = frozenset((gene_a, gene_b))
            if pair in deduped:
                continue
            deduped.add(pair)
            try:
                score = float(row.get("score") or 0.0)
            except (TypeError, ValueError):
                score = 0.0
            try:
                experimental = float(row.get("escore") or 0.0)
            except (TypeError, ValueError):
                experimental = 0.0
            try:
                database = float(row.get("dscore") or 0.0)
            except (TypeError, ValueError):
                database = 0.0
            edges.append({
                "gene_a": gene_a,
                "gene_b": gene_b,
                "combined_score": round(score, 3),
                "confidence": _confidence_band(score),
                "experimental_score": round(experimental, 3),
                "database_score": round(database, 3),
                "has_experimental_support": experimental > 0 or database > 0,
            })
        edges.sort(key=lambda e: -e["combined_score"])
        return edges

    @staticmethod
    def _parse_enrichment(payload: Any) -> List[Dict[str, Any]]:
        terms: List[Dict[str, Any]] = []
        for row in payload or []:
            category = str(row.get("category") or "")
            if category not in ENRICHMENT_CATEGORIES:
                continue
            try:
                fdr = float(row.get("fdr"))
            except (TypeError, ValueError):
                continue
            if fdr > 0.05:
                continue
            genes = row.get("inputGenes")
            gene_list = (
                [g.strip().upper() for g in genes.split(",") if g.strip()]
                if isinstance(genes, str) else
                [str(g).upper() for g in (genes or [])]
            )
            terms.append({
                "category": category,
                "category_label": ENRICHMENT_CATEGORIES[category],
                "term": row.get("term"),
                "description": row.get("description"),
                "fdr": fdr,
                "observed_gene_count": row.get("number_of_genes"),
                "background_gene_count": row.get("number_of_genes_in_background"),
                "genes": gene_list,
            })
        # Most significant first; cap so the UI stays readable.
        terms.sort(key=lambda t: t["fdr"])
        return terms[:40]

    def _analyze_graph(
        self, edges: List[Dict[str, Any]], seeds: Set[str],
    ) -> Dict[str, Any]:
        """Compute centrality, hub ranking, and community structure."""
        if _HAVE_NETWORKX:
            return self._analyze_with_networkx(edges, seeds)
        return self._analyze_degree_only(edges, seeds)

    def _analyze_with_networkx(
        self, edges: List[Dict[str, Any]], seeds: Set[str],
    ) -> Dict[str, Any]:
        graph = nx.Graph()
        for edge in edges:
            graph.add_edge(
                edge["gene_a"], edge["gene_b"],
                weight=edge["combined_score"],
            )

        degree = dict(graph.degree())
        betweenness = nx.betweenness_centrality(graph, weight=None)
        closeness = nx.closeness_centrality(graph)
        try:
            eigenvector = nx.eigenvector_centrality(graph, max_iter=500, tol=1e-6)
        except (nx.PowerIterationFailedConvergence, nx.NetworkXException):
            eigenvector = {node: 0.0 for node in graph.nodes()}
        clustering = nx.clustering(graph)

        max_degree = max(degree.values()) if degree else 0
        nodes: List[Dict[str, Any]] = []
        for node in graph.nodes():
            nodes.append({
                "gene_symbol": node,
                "is_seed": node in seeds,
                "degree": degree.get(node, 0),
                "degree_centrality": round(
                    degree.get(node, 0) / (graph.number_of_nodes() - 1), 4,
                ) if graph.number_of_nodes() > 1 else 0.0,
                "betweenness_centrality": round(betweenness.get(node, 0.0), 4),
                "closeness_centrality": round(closeness.get(node, 0.0), 4),
                "eigenvector_centrality": round(eigenvector.get(node, 0.0), 4),
                "clustering_coefficient": round(clustering.get(node, 0.0), 4),
            })

        # Composite hub score: degree, brokerage, and influence weighted equally
        # after normalising each to 0-1 within this network.
        def _norm(values: List[float]) -> Dict[str, float]:
            top = max(values) if values else 0.0
            return {"max": top or 1.0}

        deg_norm = _norm([n["degree"] for n in nodes])
        btw_norm = _norm([n["betweenness_centrality"] for n in nodes])
        eig_norm = _norm([n["eigenvector_centrality"] for n in nodes])

        for node in nodes:
            node["hub_score"] = round(
                (node["degree"] / deg_norm["max"]) * 0.4
                + (node["betweenness_centrality"] / btw_norm["max"]) * 0.3
                + (node["eigenvector_centrality"] / eig_norm["max"]) * 0.3,
                4,
            )

        nodes.sort(key=lambda n: -n["hub_score"])
        hubs = [
            {
                **node,
                "rationale": self._hub_rationale(node, max_degree),
            }
            for node in nodes[:15]
        ]

        modules = self._detect_modules(graph, seeds)

        return {
            "nodes": nodes,
            "hubs": hubs,
            "modules": modules,
            "node_count": graph.number_of_nodes(),
            "edge_count": graph.number_of_edges(),
            "network_density": round(nx.density(graph), 4),
            "connected_components": nx.number_connected_components(graph),
            "average_clustering": round(nx.average_clustering(graph), 4),
            "metrics_engine": "networkx",
        }

    @staticmethod
    def _hub_rationale(node: Dict[str, Any], max_degree: int) -> str:
        parts = [f"{node['degree']} interaction partner(s)"]
        if max_degree and node["degree"] == max_degree:
            parts.append("highest-degree node in this network")
        if node["betweenness_centrality"] > 0:
            parts.append(
                f"bridges distinct network regions (betweenness "
                f"{node['betweenness_centrality']})"
            )
        if node["is_seed"]:
            parts.append("carries a coding variant in this sample")
        else:
            parts.append("added from the STRING first shell, not mutated here")
        return "; ".join(parts)

    @staticmethod
    def _detect_modules(graph: "nx.Graph", seeds: Set[str]) -> List[Dict[str, Any]]:
        """Greedy modularity communities = candidate functional modules."""
        try:
            communities = nx.community.greedy_modularity_communities(graph)
        except Exception:
            return []
        modules = []
        for index, community in enumerate(communities, start=1):
            members = sorted(str(m) for m in community)
            if len(members) < 2:
                continue
            modules.append({
                "module_id": index,
                "size": len(members),
                "genes": members,
                "seed_genes": sorted(m for m in members if m in seeds),
            })
        modules.sort(key=lambda m: -m["size"])
        return modules[:10]

    @staticmethod
    def _analyze_degree_only(
        edges: List[Dict[str, Any]], seeds: Set[str],
    ) -> Dict[str, Any]:
        degree: Dict[str, int] = {}
        for edge in edges:
            degree[edge["gene_a"]] = degree.get(edge["gene_a"], 0) + 1
            degree[edge["gene_b"]] = degree.get(edge["gene_b"], 0) + 1

        node_count = len(degree)
        max_degree = max(degree.values()) if degree else 1
        nodes = [
            {
                "gene_symbol": gene,
                "is_seed": gene in seeds,
                "degree": count,
                "degree_centrality": round(count / (node_count - 1), 4) if node_count > 1 else 0.0,
                "betweenness_centrality": None,
                "closeness_centrality": None,
                "eigenvector_centrality": None,
                "clustering_coefficient": None,
                "hub_score": round(count / max_degree, 4),
            }
            for gene, count in degree.items()
        ]
        nodes.sort(key=lambda n: -n["degree"])
        possible_edges = node_count * (node_count - 1) / 2
        return {
            "nodes": nodes,
            "hubs": [
                {**node, "rationale": f"{node['degree']} interaction partner(s)"}
                for node in nodes[:15]
            ],
            "modules": [],
            "node_count": node_count,
            "edge_count": len(edges),
            "network_density": round(len(edges) / possible_edges, 4) if possible_edges else 0.0,
            "connected_components": None,
            "average_clustering": None,
            "metrics_engine": "degree-only-fallback",
        }
