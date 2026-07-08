# Pharmacogenomics MCP Pipeline

A comprehensive Model Context Protocol (MCP) server for pharmacogenomics analysis. This pipeline analyzes genetic variants from VCF files and provides insights into drug interactions, protein structure impacts, and clinical significance.

## Features

- **VCF Quality Control**: Comprehensive QC metrics and variant filtering
- **Variant Annotation**: Integration with Ensembl VEP for gene and consequence annotation
- **Protein Structure Analysis**: UniProt and AlphaFold integration for 3D structure context
- **Drug Interaction Matching**: DrugBank and PharmGKB integration for pharmacogenomic insights
- **Interactive Reports**: JSON and HTML report generation with visualizations
- **MCP Server**: Full integration with Kiro development environment

## Pipeline Overview

```
VCF Input → QC Analysis → Variant Annotation (VEP) → Protein Structure (AlphaFold) → Drug Matching → Report Generation
```

### Phase 1: VCF Parsing & QC
- Parse VCF files with cyvcf2
- Quality metrics: QUAL, depth, Ts/Tv ratio, genotype statistics
- High/low confidence classification

### Phase 2: Variant Annotation
- Ensembl VEP REST API integration
- Gene symbols, consequence terms, HGVSp notation
- Canonical transcript prioritization

### Phase 3: Protein Structure Analysis
- UniProt protein information and domains
- AlphaFold structure integration with pLDDT confidence scores
- Variant position impact prediction

### Phase 4: Drug Matching
- DrugBank API integration (requires license)
- PharmGKB pharmacogenomic annotations
- Evidence-based drug ranking

### Phase 5: Report Generation
- Comprehensive HTML and JSON reports
- Interactive 3D structure visualizations
- Clinical significance assessment

## Installation

1. **Clone and setup**:
```bash
cd /home/venkatesh-g/Documents/minki
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

2. **Initialize the pipeline**:
```bash
python setup.py
```

3. **Configure API access** (edit `config.json`):
```json
{
  "drugbank": {
    "api_key": "your_drugbank_api_key_here"
  }
}
```

## Usage

### As MCP Server (with Kiro)

1. **Start the MCP server**:
```bash
python pharmacogenomics_mcp.py
```

2. **Available MCP Tools**:
- `parse_vcf`: Parse VCF file and perform QC
- `annotate_variants`: Annotate variants with VEP
- `get_protein_structure`: Analyze protein structure impact
- `match_drugs`: Find drug interactions
- `generate_report`: Create comprehensive reports
- `run_full_pipeline`: Execute complete analysis

### Standalone Usage

```python
from pharmacogenomics_mcp import PharmacogenomicsMCP

# Initialize pipeline
pipeline = PharmacogenomicsMCP()

# Run full analysis
results = await pipeline.run_full_pipeline(
    vcf_file="path/to/sample.vcf",
    sample_id="SAMPLE_001"
)

print(f"Generated reports: {results['report_result']['report_files']}")
```

### Testing

```bash
python test_pipeline.py
```

## Configuration

The `config.json` file contains all API endpoints and settings:

```json
{
  "ensembl_vep": {
    "base_url": "https://rest.ensembl.org/vep/human/region",
    "rate_limit_delay": 0.1
  },
  "uniprot": {
    "base_url": "https://rest.uniprot.org"
  },
  "alphafold": {
    "base_url": "https://alphafold.ebi.ac.uk/files"
  },
  "drugbank": {
    "api_url": "https://go.drugbank.com/api/v1",
    "api_key": "your_api_key_here"
  },
  "qc_thresholds": {
    "min_qual": 30,
    "min_depth": 10
  }
}
```

## API Limitations

### Without DrugBank License
- Limited to free DrugBank API calls (rate limited)
- Uses cached pharmacogene data for common genes (CYP2D6, CYP2C19, etc.)
- PharmGKB integration provides pharmacogenomic annotations

### With DrugBank License
- Full access to comprehensive drug-target relationships
- Batch processing capabilities
- Enhanced drug interaction analysis

## Example VCF

The pipeline includes an example VCF with known pharmacogenomic variants:
- CYP2D6 variants (codeine, tramadol metabolism)
- CYP2C19 variants (clopidogrel activation)
- CYP2C9 variants (warfarin sensitivity)
- MTHFR variants (folate metabolism)

## Output Reports

### HTML Report Sections
1. **Quality Control Summary**: Variant statistics and quality metrics
2. **Pharmacogenomic Summary**: Genes affected and drug interactions
3. **Detailed Variant Analysis**: Per-variant breakdown with structure context
4. **Clinical Recommendations**: Evidence-based guidance

### JSON Report Structure
```json
{
  "sample_information": {...},
  "quality_control": {...},
  "pharmacogenomic_summary": {...},
  "detailed_results": [
    {
      "variant_info": {...},
      "gene_annotation": {...},
      "protein_analysis": {...},
      "drug_interactions": {...},
      "clinical_significance": "high"
    }
  ]
}
```

## Known Pharmacogenes Included

- **CYP2D6**: Codeine, tramadol, metoprolol, risperidone
- **CYP2C19**: Clopidogrel, omeprazole, escitalopram  
- **CYP2C9**: Warfarin, phenytoin, celecoxib
- **TPMT**: Azathioprine, mercaptopurine
- **DPYD**: Fluorouracil, capecitabine
- **VKORC1**: Warfarin sensitivity

## Troubleshooting

### Common Issues

1. **cyvcf2 installation fails**:
   - Install system dependencies: `sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev`
   - Or use conda: `conda install -c bioconda cyvcf2`

2. **VEP API rate limits**:
   - Increase `rate_limit_delay` in config.json
   - Process variants in smaller batches

3. **Missing structure files**:
   - AlphaFold structures are downloaded automatically
   - Check internet connection for structure downloads

4. **No drug matches**:
   - Ensure DrugBank API key is configured
   - Check that gene symbols are standard HGNC symbols

## Development

### Project Structure
```
minki/
├── pharmacogenomics_mcp.py      # Main MCP server
├── config.json                  # Configuration
├── requirements.txt             # Dependencies
├── setup.py                     # Setup script
├── test_pipeline.py             # Test suite
├── modules/                     # Core modules
│   ├── vcf_parser.py           # VCF parsing and QC
│   ├── vep_annotator.py        # Ensembl VEP integration
│   ├── protein_structure.py    # UniProt/AlphaFold analysis
│   ├── drug_matcher.py         # DrugBank/PharmGKB matching
│   └── report_generator.py     # Report generation
├── examples/                    # Example files
└── results/                     # Output directory
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License and Citations

### Software License
This project is open source. See license file for details.

### Data Sources and Citations

Please cite the following when using this pipeline:

- **Ensembl VEP**: Cunningham et al. (2022) Ensembl 2022. Nucleic Acids Research.
- **AlphaFold**: Jumper et al. (2021) Highly accurate protein structure prediction with AlphaFold. Nature.
- **UniProt**: The UniProt Consortium (2021) UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Research.
- **PharmGKB**: Whirl-Carrillo et al. (2012) Pharmacogenomics knowledge for personalized medicine. Clinical Pharmacology & Therapeutics.
- **DrugBank**: Wishart et al. (2018) DrugBank 5.0: a major update to the DrugBank database for 2018. Nucleic Acids Research.

## Disclaimer

**Important**: This tool is for research purposes only and should not be used for clinical decision-making without proper validation and consultation with healthcare professionals. Pharmacogenomic predictions are based on current knowledge and may not account for all factors affecting drug response.

## Contact and Support

For questions, issues, or contributions, please open an issue in the repository or contact the development team.