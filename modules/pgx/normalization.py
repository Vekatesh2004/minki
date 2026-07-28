"""Build-aware normalization helpers for native PGx definition matching."""

from dataclasses import dataclass, field
from typing import List, Tuple

from .models import GenomeBuild


@dataclass(frozen=True)
class NormalizationResult:
    build: GenomeBuild
    chrom: str
    pos: int
    ref: str
    alt: str
    fully_normalized: bool
    warnings: Tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self):
        return {
            "build": self.build.value, "chrom": self.chrom, "pos": self.pos,
            "ref": self.ref, "alt": self.alt,
            "fully_normalized": self.fully_normalized,
            "warnings": list(self.warnings),
        }


def canonical_chromosome(chrom: str) -> str:
    value = str(chrom).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    if value == "M":
        value = "MT"
    return value


def normalize_variant(build, chrom: str, pos: int, ref: str, alt: str) -> NormalizationResult:
    """Minimally normalize one allele without pretending to left-align it.

    Common suffixes/prefixes are trimmed while retaining a VCF anchor base.
    Indels require a matching reference FASTA for true left alignment, so they
    are explicitly marked as not fully normalized here.
    """
    genome_build = build if isinstance(build, GenomeBuild) else GenomeBuild.parse(build)
    warnings: List[str] = []
    ref, alt = str(ref).upper(), str(alt).upper()
    if genome_build is GenomeBuild.UNKNOWN:
        warnings.append("Genome build is unknown; definition matching must not run")
    if not ref or not alt or alt == ".":
        warnings.append("Allele is empty or missing")
        return NormalizationResult(genome_build, canonical_chromosome(chrom), pos, ref, alt, False, tuple(warnings))

    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1

    fully_normalized = len(ref) == len(alt)
    if not fully_normalized:
        warnings.append("Indel is minimally trimmed but not reference-left-aligned")
    return NormalizationResult(
        genome_build, canonical_chromosome(chrom), int(pos), ref, alt,
        fully_normalized, tuple(warnings),
    )
