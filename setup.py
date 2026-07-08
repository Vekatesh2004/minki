#!/usr/bin/env python3
"""
Setup script for Pharmacogenomics MCP Pipeline
"""

import asyncio
import sys
from pathlib import Path
from modules.drug_matcher import DrugMatcher
import json

async def setup_pipeline():
    """Setup the pharmacogenomics pipeline"""
    
    print("Setting up Pharmacogenomics MCP Pipeline...")
    
    # Load configuration
    try:
        with open("config.json", 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print("Error: config.json not found!")
        return False
        
    # Initialize drug matcher and populate known pharmacogenes
    print("Initializing drug database with known pharmacogenes...")
    drug_matcher = DrugMatcher(config)
    await drug_matcher.populate_known_pharmacogenes()
    
    # Create example directories
    print("Creating output directories...")
    Path("results").mkdir(exist_ok=True)
    Path("examples").mkdir(exist_ok=True)
    
    # Create example VCF file
    print("Creating example VCF file...")
    create_example_vcf()
    
    print("Setup completed successfully!")
    print("\nTo test the pipeline:")
    print("1. Install dependencies: pip install -r requirements.txt")
    print("2. Run: python pharmacogenomics_mcp.py")
    print("3. Test with example VCF: examples/sample_pharmacogenomics.vcf")
    
    return True

def create_example_vcf():
    """Create an example VCF file with known pharmacogenomic variants"""
    
    vcf_content = """##fileformat=VCFv4.2
##source=PharmacogenomicsExample
##reference=GRCh38
##contig=<ID=22,length=50818468>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
22	42126611	rs16947	A	G	60	PASS	DP=30	GT:DP	1/1:30
22	42127803	rs1065852	G	A	85	PASS	DP=45	GT:DP	0/1:45
22	42128936	rs3892097	C	T	92	PASS	DP=38	GT:DP	0/1:38
19	15990431	rs4149056	T	C	78	PASS	DP=42	GT:DP	1/1:42
19	15881965	rs2306283	A	G	65	PASS	DP=35	GT:DP	0/1:35
10	96522463	rs1799853	C	T	88	PASS	DP=50	GT:DP	0/1:50
10	96541616	rs1057910	A	C	95	PASS	DP=48	GT:DP	0/1:48
6	18130918	rs1801131	T	G	70	PASS	DP=40	GT:DP	0/1:40
6	18143955	rs1801133	C	T	82	PASS	DP=44	GT:DP	0/1:44
"""
    
    vcf_path = Path("examples/sample_pharmacogenomics.vcf")
    vcf_path.parent.mkdir(exist_ok=True)
    
    with open(vcf_path, 'w') as f:
        f.write(vcf_content)
        
    print(f"Created example VCF: {vcf_path}")

if __name__ == "__main__":
    asyncio.run(setup_pipeline())