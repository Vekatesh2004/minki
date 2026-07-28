"""Validation and loading for versioned, provenance-aware PGx definitions."""

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Union

from .models import (
    AlleleDefinition, DefinitionSource, GeneDefinition, GenomeBuild,
    LocusDefinition,
)

SUPPORTED_SCHEMA_VERSIONS = {"1.0"}


class DefinitionValidationError(ValueError):
    """Raised when a definition could produce scientifically unsafe calls."""


class DefinitionCatalog:
    def __init__(self, genes: Iterable[GeneDefinition] = ()):
        self._genes = {gene.gene: gene for gene in genes}

    @property
    def genes(self):
        return tuple(sorted(self._genes))

    def get(self, gene: str) -> GeneDefinition:
        try:
            return self._genes[gene.upper()]
        except KeyError as exc:
            raise DefinitionValidationError(f"No validated definition for {gene}") from exc

    @classmethod
    def from_json(cls, source: Union[str, Path, Dict[str, Any]]):
        if isinstance(source, dict):
            payload = source
        else:
            path = Path(source)
            payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(_parse_gene(item) for item in payload.get("genes", []))


def _required(data: Dict[str, Any], key: str, context: str):
    value = data.get(key)
    if value is None or value == "":
        raise DefinitionValidationError(f"Missing {key} in {context}")
    return value


def _parse_gene(data: Dict[str, Any]) -> GeneDefinition:
    gene = str(_required(data, "gene", "gene definition")).upper()
    schema_version = str(_required(data, "schema_version", gene))
    if schema_version not in SUPPORTED_SCHEMA_VERSIONS:
        raise DefinitionValidationError(f"Unsupported schema {schema_version} for {gene}")
    build = GenomeBuild.parse(_required(data, "genome_build", gene))
    if build is GenomeBuild.UNKNOWN:
        raise DefinitionValidationError(f"Unknown genome build for {gene}")

    raw_source = data.get("source") or {}
    source = DefinitionSource(
        name=str(_required(raw_source, "name", f"{gene}.source")),
        version=str(_required(raw_source, "version", f"{gene}.source")),
        url=str(_required(raw_source, "url", f"{gene}.source")),
        retrieved_at=str(_required(raw_source, "retrieved_at", f"{gene}.source")),
        checksum=raw_source.get("checksum"),
        license=raw_source.get("license"),
    )

    loci = tuple(_parse_locus(gene, build, item) for item in data.get("loci", []))
    locus_ids = {locus.id for locus in loci}
    if len(locus_ids) != len(loci):
        raise DefinitionValidationError(f"Duplicate locus IDs for {gene}")
    alleles = tuple(_parse_allele(gene, item, locus_ids) for item in data.get("alleles", []))
    if not loci or not alleles:
        raise DefinitionValidationError(f"{gene} must contain loci and alleles")
    return GeneDefinition(gene, schema_version, build, source, loci, alleles)


def _parse_locus(gene: str, build: GenomeBuild, data: Dict[str, Any]) -> LocusDefinition:
    locus_id = str(_required(data, "id", f"{gene}.locus"))
    ref = str(_required(data, "ref", locus_id)).upper()
    alt = str(_required(data, "alt", locus_id)).upper()
    if ref == alt:
        raise DefinitionValidationError(f"REF and ALT are identical for {locus_id}")
    return LocusDefinition(
        id=locus_id,
        build=build,
        chrom=str(_required(data, "chrom", locus_id)).removeprefix("chr"),
        pos=int(_required(data, "pos", locus_id)),
        ref=ref,
        alt=alt,
        rsid=data.get("rsid"),
        required=bool(data.get("required", True)),
    )


def _parse_allele(gene: str, data: Dict[str, Any], locus_ids) -> AlleleDefinition:
    name = str(_required(data, "name", f"{gene}.allele"))
    variants = tuple(str(value) for value in data.get("variants", []))
    unknown = set(variants) - locus_ids
    if unknown:
        raise DefinitionValidationError(
            f"{gene} {name} references unknown loci: {sorted(unknown)}"
        )
    score = data.get("activity_score")
    return AlleleDefinition(
        name=name,
        variants=variants,
        function=data.get("function"),
        activity_score=float(score) if score is not None else None,
        is_reference=bool(data.get("is_reference", False)),
    )
