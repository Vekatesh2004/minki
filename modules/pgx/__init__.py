"""Native, provenance-aware pharmacogenomics data contracts.

This package contains Minki-owned implementations. Reference repositories are
review inputs only and are never imported or executed here.
"""

from .definitions import DefinitionCatalog, DefinitionValidationError
from .models import (
    AlleleDefinition,
    CallState,
    DefinitionSource,
    GeneDefinition,
    GenomeBuild,
    GenotypeCall,
    LocusDefinition,
    VariantRecord,
    VCFDocument,
    VCFMetadata,
)
from .normalization import NormalizationResult, normalize_variant
from .vcf import VCFInputError, read_vcf_document

__all__ = [
    "AlleleDefinition", "CallState", "DefinitionCatalog",
    "DefinitionSource", "DefinitionValidationError", "GeneDefinition",
    "GenomeBuild", "GenotypeCall", "LocusDefinition",
    "NormalizationResult", "VCFDocument", "VCFInputError",
    "VCFMetadata", "VariantRecord", "normalize_variant",
    "read_vcf_document",
]
