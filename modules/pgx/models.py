"""Typed native PGx models that preserve uncertainty and provenance."""

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


class GenomeBuild(str, Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    UNKNOWN = "unknown"

    @classmethod
    def parse(cls, value: Optional[str]) -> "GenomeBuild":
        normalized = str(value or "").strip().lower().replace("_", "")
        aliases = {
            "grch37": cls.GRCH37, "hg19": cls.GRCH37, "b37": cls.GRCH37,
            "grch38": cls.GRCH38, "hg38": cls.GRCH38, "b38": cls.GRCH38,
        }
        return aliases.get(normalized, cls.UNKNOWN)


class CallState(str, Enum):
    CALLED = "called"
    NO_CALL = "no_call"
    PARTIAL_NO_CALL = "partial_no_call"
    NOT_PRESENT = "not_present"
    NOT_ASSESSED = "not_assessed"


@dataclass(frozen=True)
class GenotypeCall:
    sample_id: str
    allele_indices: Tuple[Optional[int], ...]
    phased: bool
    phase_set: Optional[str] = None
    raw_gt: str = "."
    fields: Dict[str, Any] = field(default_factory=dict)

    @property
    def state(self) -> CallState:
        if not self.allele_indices or all(a is None for a in self.allele_indices):
            return CallState.NO_CALL
        if any(a is None for a in self.allele_indices):
            return CallState.PARTIAL_NO_CALL
        return CallState.CALLED

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["allele_indices"] = list(self.allele_indices)
        data["state"] = self.state.value
        return data


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alts: Tuple[str, ...]
    record_id: Optional[str] = None
    qual: Optional[float] = None
    filter_status: Optional[str] = None
    info_raw: str = "."
    format_keys: Tuple[str, ...] = ()
    calls: Tuple[GenotypeCall, ...] = ()


@dataclass(frozen=True)
class VCFMetadata:
    file_format: Optional[str]
    reference: Optional[str]
    genome_build: GenomeBuild
    sample_ids: Tuple[str, ...]
    headers: Tuple[str, ...] = ()
    build_source: str = "unknown"
    warnings: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["genome_build"] = self.genome_build.value
        data["sample_ids"] = list(self.sample_ids)
        data["headers"] = list(self.headers)
        data["warnings"] = list(self.warnings)
        return data


@dataclass(frozen=True)
class VCFDocument:
    metadata: VCFMetadata
    records: Tuple[VariantRecord, ...]


@dataclass(frozen=True)
class DefinitionSource:
    name: str
    version: str
    url: str
    retrieved_at: str
    checksum: Optional[str] = None
    license: Optional[str] = None


@dataclass(frozen=True)
class LocusDefinition:
    id: str
    build: GenomeBuild
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: Optional[str] = None
    required: bool = True


@dataclass(frozen=True)
class AlleleDefinition:
    name: str
    variants: Tuple[str, ...]
    function: Optional[str] = None
    activity_score: Optional[float] = None
    is_reference: bool = False


@dataclass(frozen=True)
class GeneDefinition:
    gene: str
    schema_version: str
    genome_build: GenomeBuild
    source: DefinitionSource
    loci: Tuple[LocusDefinition, ...]
    alleles: Tuple[AlleleDefinition, ...]
