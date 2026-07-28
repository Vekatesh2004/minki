# Native PGx Engine: Repository Review and Adoption Record

## Scope and safety

This review was performed by static inspection only. No reference-repository scripts, JARs, containers, installers, models, or databases were executed. The repositories remain under `reference_repositories/`, which is Git-ignored. Minki will implement selected concepts natively; it will not wrap these projects or silently redistribute their code or knowledge assets.

Reviewed checkouts:

| Project | Commit | Checkout commit date | Code license |
|---|---|---:|---|
| PharmCAT | `3ec443f7c1de79486fd4d8d47543b0e9ccc10b35` | 2026-07-19 | MPL-2.0 |
| PyPGx | `2ef6e014389ac8a060a81964136961b7240ac856` | 2026-06-12 | MIT; bundled Beagle is GPL and data/resources require separate review |
| P3 | `ea7ddb8457adada99627f41c3e6e76af264ce496` | 2018-10-05 | CC0-1.0 for repository work; upstream data/tools retain their terms |
| PAnno | `0e09d28b6b819ce196677b8d832bf6779a76ae75` | 2022-12-28 | MPL-2.0 |

Commit dates describe these shallow checkouts, not a complete maintenance audit.

## Capability comparison

| Capability | PharmCAT | PyPGx 0.27.0 | P3 | PAnno 0.3.1 | Native decision |
|---|---|---|---|---|---|
| Named/star alleles | Deterministic, definition-driven; normal matcher covers a limited gene set | Small-variant definitions for a broad target set | None | Hard-coded inference for 21 genes | Build a versioned gene-agnostic schema |
| Phasing | Preserves phased/unphased and phase-set semantics; ambiguity matters | Uses existing phase or Beagle statistical phasing | Not a diplotype caller | Enumerates candidates but loses important phase semantics | Preserve GT delimiter, PS, alternatives, and uncertainty |
| Missingness | Explicit missing-position/no-call behavior and gene-specific rules | Semantic archives, imported/consolidated states | Not applicable | Unsafe absent-as-reference behavior | Never equate absent, no-call, and hom-reference without callability evidence |
| GRCh builds | Primarily GRCh38 workflow | GRCh37 and GRCh38 | Legacy feature inputs, not PGx definitions | GRCh38 only | Build is mandatory for matching and provenance |
| CNV/SV | CYP2D6 generally requires an outside call; not a general caller | Specialized models/rules for 13 genes; external bundle required | Generic cell-line CNV features, not pharmacogene structural calling | CYP2D6 handling is insufficient | Separate evidence contract; implement complex genes last |
| Phenotypes/guidance | CPIC/DPWG/FDA reporting with versioned data | Phenotype mappings for a subset; recommendation tables | Predicts in-vitro cell-line drug sensitivity | Bundled SQLite rules and annotations | Separate validated phenotype/guideline layers with source versions |
| Validation | Mature project semantics and extensive domain handling | Limited repository tests relative to scope | No meaningful patient-PGx validation suite found | Tests effectively nonfunctional | Independent authoritative fixtures required before parity claims |

## Project assessments

### PharmCAT — reference semantics, do not copy by default

The most valuable concepts are required-locus checks, strict separation of sample missingness from definition defaults, explicit ambiguity, phased versus unphased matching, gene-specific exemptions, and source/data version reporting. Its normal matcher does not solve general CNV/SV detection, and complex loci such as CYP2D6/HLA/MT-RNR1 need outside or specialized calls. Direct source/data copying would create MPL obligations and separate clinical-data provenance work, so Minki will implement the concepts independently.

### PyPGx — reference architecture and complex-locus research

PyPGx has the broadest gene target list in the reviewed projects and makes build/platform/semantic type explicit. It demonstrates that small variants, phasing, copy-number evidence, final genotype rules, and phenotype mapping must remain separate stages. The NGS path depends on a matching external `pypgx-bundle`; Beagle has its own GPL terms; trained classifier and clinical table provenance must be assessed independently. Model archives must be treated as trusted deployment assets, never accepted from untrusted uploads. Minki will not copy tables, models, JARs, or rules.

### P3 — exclude from star-allele implementation

P3 is an older Snakemake/R research workflow for predicting **in-vitro cell-line drug sensitivity**, not patient pharmacogene diplotypes. It combines expression, pathway, exome, and generic CNV features with R SuperLearner (random forest, glmnet, and mean learners). The checkout pins Python 3.4, pandas 0.16.2, Snakemake 3.4.2, and bedtools 2.24.0; some inputs require manual downloads and external datasets. The workflow saves per-compound `.RData` models but the repository does not provide a modern reproducible training manifest, model cards, held-out clinical validation, or a meaningful automated test suite. It is therefore excluded from allele/phenotype calling. Only high-level ideas—modular feature provenance and explicit training manifests—are worth retaining for any future, separately validated research prediction module.

### PAnno — design reference only

PAnno's static GRCh38 pipeline combines 21 hard-coded gene definitions, population-frequency tie-breaking, a bundled SQLite clinical database, and HTML reporting. Static review found unsafe absent-as-reference scoring and likely defects around HLA, guideline matching, ambiguity, and CLI behavior. CYP2D6 CNV/hybrid support is not robust. Its code and bundled knowledge are not a ground-truth source and will not be copied.

## Native architecture

The target layers are:

1. `modules/pgx/models.py`: immutable records for build, genotype, call state, definitions, and provenance.
2. `modules/pgx/vcf.py`: all ALT alleles, all samples, GT phasing, PS, raw FORMAT values, and explicit metadata.
3. `modules/pgx/normalization.py`: build-aware canonical representation; indels remain unmatchable until reference-backed left alignment is available.
4. `modules/pgx/definitions.py`: strict versioned schema validation; no unsupported gene is inferred.
5. Future `matcher.py`: deterministic candidate generation that retains all valid calls and evidence gaps.
6. Future `phenotype.py`: independently versioned activity-score/diplotype mappings.
7. Future `guidelines.py`: source-specific recommendations with publication/version/evidence provenance.
8. Future structural layer: depth/copy number, breakpoint, paralog/hybrid, duplication, deletion, and uncertainty evidence.

Genome-wide Ensembl VEP, UniProt, AlphaFold, and current gene-level drug evidence remain useful for genes without validated star-allele systems, but they must not be presented as genotype-based CPIC/DPWG guidance.

## Phased delivery

1. Foundation: explicit GRCh37/GRCh38/unknown, multiallelic and multisample calls, GT/PS preservation, missing/reference/no-call distinction, schema and provenance validation.
2. Initial independently verified definitions: CYP2C19, CYP2C9, SLCO1B1, TPMT.
3. Expand to DPYD, UGT1A1, NUDT15, VKORC1 and explicit contracts for HLA/MT-RNR1.
4. Add validated activity score, phenotype, and guideline layers.
5. Add complex CNV/SV genes only with suitable evidence inputs.
6. Implement CYP2D6 last, after deletion/duplication/hybrid and ambiguity models are validated.

“All genes” means all pharmacogenes with authoritative, licensed, versioned definitions—not arbitrary human genes and not invented star alleles. Every unsupported or insufficiently observed case must return an explicit no-call/unsupported result.

## Scientific and clinical limitations

Minki remains a research/education system. VCF absence does not prove a reference genotype; WES may not cover regulatory or structural loci; statistical phase is uncertain; genome builds are not interchangeable; star-allele definitions and guidelines change over time; and VEP consequence, protein impact, or gene–drug association alone does not establish treatment suitability. Results require independent validation and qualified clinical interpretation before any healthcare use.
