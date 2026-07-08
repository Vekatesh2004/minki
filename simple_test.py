#!/usr/bin/env python3
"""
Simple test for core pipeline components (without MCP dependencies)
"""

import asyncio
import json
from pathlib import Path
import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from modules.vcf_parser import VCFParser
from modules.drug_matcher import DrugMatcher

async def test_core_components():
    """Test core components without MCP dependencies"""
    
    print("Testing Core Pharmacogenomics Pipeline Components")
    print("=" * 55)
    
    # Load configuration
    try:
        with open("config.json", 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print("Error: config.json not found!")
        return False
    
    # Test 1: VCF Parser
    print("\n1. Testing VCF Parser...")
    vcf_file = "examples/sample_pharmacogenomics.vcf"
    
    if not Path(vcf_file).exists():
        print(f"   ✗ Example VCF file not found: {vcf_file}")
        print("   Run setup.py first!")
        return False
    
    try:
        parser = VCFParser(config)
        results = await parser.parse_vcf(vcf_file, "TEST_SAMPLE")
        
        variants = results.get('variants', [])
        qc_summary = results.get('qc_summary', {})
        
        print(f"   ✓ Parsed {len(variants)} variants")
        print(f"   ✓ Pass rate: {qc_summary.get('pass_rate', 0):.1%}")
        print(f"   ✓ High confidence rate: {qc_summary.get('high_confidence_rate', 0):.1%}")
        
        if 'ts_tv_ratio' in qc_summary:
            print(f"   ✓ Ts/Tv ratio: {qc_summary['ts_tv_ratio']:.2f}")
            
        # Show first few variants
        print("   ✓ Example variants:")
        for i, variant in enumerate(variants[:3]):
            pos = f"{variant.get('chrom')}:{variant.get('pos')}"
            change = f"{variant.get('ref')}→{variant.get('alt')}"
            qual = variant.get('qual', 'N/A')
            conf = variant.get('confidence', 'unknown')
            print(f"      {i+1}. {pos} {change} (Q={qual}, {conf})")
            
    except Exception as e:
        print(f"   ✗ VCF parsing failed: {e}")
        return False
    
    # Test 2: Drug Matcher Database
    print("\n2. Testing Drug Matcher Database...")
    
    try:
        drug_matcher = DrugMatcher(config)
        
        # Test CYP2D6 drug matching
        cyp2d6_results = await drug_matcher.match_drugs("CYP2D6")
        matched_drugs = cyp2d6_results.get('matched_drugs', [])
        
        print(f"   ✓ Found {len(matched_drugs)} drugs for CYP2D6")
        
        if matched_drugs:
            print("   ✓ Example CYP2D6 drugs:")
            for i, drug in enumerate(matched_drugs[:3]):
                drug_name = drug.get('drug_name', 'Unknown')
                action = drug.get('action', 'Unknown')
                evidence = drug.get('evidence_level', 'no_evidence')
                print(f"      {i+1}. {drug_name} ({action}, {evidence})")
        
        # Test other pharmacogenes
        test_genes = ['CYP2C19', 'CYP2C9', 'TPMT']
        for gene in test_genes:
            gene_results = await drug_matcher.match_drugs(gene)
            gene_matches = len(gene_results.get('matched_drugs', []))
            print(f"   ✓ {gene}: {gene_matches} drug matches")
            
    except Exception as e:
        print(f"   ✗ Drug matching failed: {e}")
        return False
    
    # Test 3: Configuration Validation
    print("\n3. Testing Configuration...")
    
    required_sections = ['ensembl_vep', 'uniprot', 'alphafold', 'drugbank', 'qc_thresholds']
    missing_sections = []
    
    for section in required_sections:
        if section not in config:
            missing_sections.append(section)
        else:
            print(f"   ✓ {section} section present")
    
    if missing_sections:
        print(f"   ✗ Missing config sections: {', '.join(missing_sections)}")
        return False
    
    # Check API endpoints
    vep_url = config.get('ensembl_vep', {}).get('base_url', '')
    uniprot_url = config.get('uniprot', {}).get('base_url', '')
    alphafold_url = config.get('alphafold', {}).get('base_url', '')
    
    if all([vep_url, uniprot_url, alphafold_url]):
        print("   ✓ All API endpoints configured")
    else:
        print("   ! Some API endpoints missing (but not critical for testing)")
    
    # Test 4: Database Creation
    print("\n4. Testing Database Creation...")
    
    db_path = "drug_cache.db"
    if Path(db_path).exists():
        print(f"   ✓ Drug cache database exists: {db_path}")
        
        # Check database size
        db_size = Path(db_path).stat().st_size
        print(f"   ✓ Database size: {db_size} bytes")
        
        if db_size > 1000:  # Should have some data
            print("   ✓ Database appears to contain data")
        else:
            print("   ! Database is very small, may need population")
    else:
        print(f"   ✗ Drug cache database not found: {db_path}")
        return False
    
    print("\n" + "=" * 55)
    print("🎉 Core component testing completed successfully!")
    print("\nPipeline is ready for:")
    print("• VCF file parsing and quality control")
    print("• Drug-gene interaction matching")
    print("• Known pharmacogene analysis")
    print("\nNext steps:")
    print("1. Install full dependencies: pip install -r requirements.txt")
    print("2. Configure API keys in config.json")
    print("3. Test with real VCF files")
    print("4. Use as MCP server with Kiro")
    
    return True

async def show_example_usage():
    """Show example of how to use the pipeline"""
    
    print("\n" + "=" * 55)
    print("Example Usage")
    print("=" * 55)
    
    print("""
# Basic VCF analysis
from modules.vcf_parser import VCFParser
import json

# Load config
with open('config.json', 'r') as f:
    config = json.load(f)

# Parse VCF
parser = VCFParser(config)
results = await parser.parse_vcf('your_file.vcf', 'SAMPLE_ID')

print(f"Parsed {len(results['variants'])} variants")
print(f"Pass rate: {results['qc_summary']['pass_rate']:.1%}")

# Drug matching
from modules.drug_matcher import DrugMatcher

matcher = DrugMatcher(config)
drug_results = await matcher.match_drugs('CYP2D6')

for drug in drug_results['matched_drugs']:
    print(f"Drug: {drug['drug_name']} - Evidence: {drug['evidence_level']}")
""")

if __name__ == "__main__":
    try:
        success = asyncio.run(test_core_components())
        if success:
            asyncio.run(show_example_usage())
    except KeyboardInterrupt:
        print("\n\nTest interrupted by user")
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        sys.exit(1)