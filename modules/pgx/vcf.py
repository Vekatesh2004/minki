"""Lossless-enough VCF reader for native PGx evidence contracts.

It preserves every ALT, sample GT, phasing delimiter, phase set, and raw FORMAT
value. It does not infer reference calls for absent loci.
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .models import GenomeBuild, GenotypeCall, VariantRecord, VCFDocument, VCFMetadata


class VCFInputError(ValueError):
    pass


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open("r", encoding="utf-8")


def _detect_build(headers: List[str], explicit_build: Optional[str]):
    if explicit_build:
        build = GenomeBuild.parse(explicit_build)
        if build is GenomeBuild.UNKNOWN:
            raise VCFInputError(f"Unsupported explicit genome build: {explicit_build}")
        return build, "explicit"
    text = " ".join(headers)
    for token in re.findall(r"GRCh3[78]|hg(?:19|38)|b3[78]", text, re.IGNORECASE):
        build = GenomeBuild.parse(token)
        if build is not GenomeBuild.UNKNOWN:
            return build, "vcf_header"
    return GenomeBuild.UNKNOWN, "unknown"


def _parse_gt(sample_id: str, format_keys: Tuple[str, ...], raw_sample: str) -> GenotypeCall:
    values = raw_sample.split(":")
    fields = {key: values[i] if i < len(values) else "." for i, key in enumerate(format_keys)}
    raw_gt = fields.get("GT", ".")
    phased = "|" in raw_gt
    parts = re.split(r"[/|]", raw_gt) if raw_gt else ["."]
    alleles = tuple(None if value in ("", ".") else int(value) for value in parts)
    phase_set = fields.get("PS")
    if phase_set in (None, "", "."):
        phase_set = None
    return GenotypeCall(sample_id, alleles, phased, phase_set, raw_gt, fields)


def read_vcf_document(file_path: str, genome_build: Optional[str] = None) -> VCFDocument:
    path = Path(file_path)
    if not path.is_file():
        raise VCFInputError(f"VCF file not found: {file_path}")

    headers: List[str] = []
    sample_ids: Tuple[str, ...] = ()
    records: List[VariantRecord] = []
    column_header_seen = False
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#CHROM"):
                columns = line.split("\t")
                sample_ids = tuple(columns[9:])
                column_header_seen = True
                continue
            if not line or line.startswith("#"):
                continue
            if not column_header_seen:
                raise VCFInputError("VCF is missing the #CHROM header")
            fields = line.split("\t")
            if len(fields) < 8:
                raise VCFInputError(f"Malformed VCF record at line {line_number}")
            try:
                pos = int(fields[1])
                qual = None if fields[5] == "." else float(fields[5])
            except ValueError as exc:
                raise VCFInputError(f"Invalid POS/QUAL at line {line_number}") from exc
            format_keys = tuple(fields[8].split(":")) if len(fields) > 8 and fields[8] != "." else ()
            calls = tuple(
                _parse_gt(sample_id, format_keys, fields[9 + index] if 9 + index < len(fields) else ".")
                for index, sample_id in enumerate(sample_ids)
            )
            records.append(VariantRecord(
                chrom=fields[0], pos=pos, record_id=None if fields[2] == "." else fields[2],
                ref=fields[3].upper(), alts=tuple(value.upper() for value in fields[4].split(",")),
                qual=qual, filter_status=None if fields[6] == "." else fields[6],
                info_raw=fields[7], format_keys=format_keys, calls=calls,
            ))

    build, build_source = _detect_build(headers, genome_build)
    reference = next((line.split("=", 1)[1] for line in headers if line.lower().startswith("##reference=")), None)
    file_format = next((line.split("=", 1)[1] for line in headers if line.lower().startswith("##fileformat=")), None)
    warnings = () if build is not GenomeBuild.UNKNOWN else (
        "Genome build is not declared; star-allele matching is disabled until supplied explicitly",
    )
    metadata = VCFMetadata(file_format, reference, build, sample_ids, tuple(headers), build_source, warnings)
    return VCFDocument(metadata, tuple(records))
