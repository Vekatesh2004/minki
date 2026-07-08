#!/usr/bin/env python3
"""
Pharmacogenomics MCP Server
A comprehensive pipeline for variant annotation and drug interaction analysis
"""

import asyncio
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Resource, Tool, TextContent, ImageContent
import mcp.types as types

# Import our pipeline modules
from modules.vcf_parser import VCFParser
from modules.vep_annotator import VEPAnnotator
from modules.protein_structure import ProteinStructureAnalyzer
from modules.drug_matcher import DrugMatcher
from modules.report_generator import ReportGenerator

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PharmacogenomicsMCP:
    def __init__(self, config_path: str = "config.json"):
        """Initialize the MCP server with configuration"""
        self.config = self._load_config(config_path)
        self.vcf_parser = VCFParser(self.config)
        self.vep_annotator = VEPAnnotator(self.config)
        self.structure_analyzer = ProteinStructureAnalyzer(self.config)
        self.drug_matcher = DrugMatcher(self.config)
        self.report_generator = ReportGenerator(self.config)
        
    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """Load configuration from JSON file"""
        try:
            with open(config_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            logger.error(f"Config file {config_path} not found")
            raise
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in config file: {e}")
            raise
            
    async def run_full_pipeline(self, vcf_file: str, sample_id: str, output_dir: str = "results") -> Dict[str, Any]:
        """
        Run the complete pharmacogenomics analysis pipeline
        
        Args:
            vcf_file: Path to VCF file
            sample_id: Sample identifier  
            output_dir: Output directory for reports
            
        Returns:
            Dictionary containing all analysis results and report paths
        """
        try:
            logger.info(f"Starting full pipeline analysis for {sample_id}")
            
            # Ensure output directory exists
            Path(output_dir).mkdir(exist_ok=True)
            
            # Phase 1: Parse VCF and QC
            logger.info("Phase 1: Parsing VCF and performing QC...")
            vcf_results = await self.vcf_parser.parse_vcf(vcf_file, sample_id)
            
            # Phase 2: Annotate variants with VEP
            logger.info("Phase 2: Annotating variants with VEP...")
            variants = vcf_results.get("variants", [])
            annotation_results = await self.vep_annotator.annotate_variants(variants)
            
            # Phase 3: Analyze protein structures for coding variants
            logger.info("Phase 3: Analyzing protein structures...")
            structure_results = []
            drug_results = []
            
            for variant in annotation_results:
                vep_annotation = variant.get("vep_annotation", {})
                
                # Only process coding variants
                if not vep_annotation.get("coding", False):
                    continue
                    
                gene_symbol = vep_annotation.get("gene_symbol")
                uniprot_ids = vep_annotation.get("uniprot_ids", [])
                hgvsp = vep_annotation.get("hgvsp")
                amino_acid_pos = vep_annotation.get("amino_acid_position")
                
                # Get UniProt IDs if not available from VEP
                if not uniprot_ids and gene_symbol:
                    ensembl_gene_id = vep_annotation.get("gene_id")
                    uniprot_ids = await self.structure_analyzer.map_gene_to_uniprot(gene_symbol, ensembl_gene_id)
                    
                # Analyze structure for each UniProt ID
                for uniprot_id in uniprot_ids[:1]:  # Limit to first match
                    if amino_acid_pos:
                        structure_result = await self.structure_analyzer.analyze_structure(
                            uniprot_id, amino_acid_pos, hgvsp
                        )
                        structure_results.append(structure_result)
                        
                # Phase 4: Match drugs for this gene
                if gene_symbol:
                    primary_uniprot = uniprot_ids[0] if uniprot_ids else None
                    drug_result = await self.drug_matcher.match_drugs(gene_symbol, primary_uniprot, hgvsp)
                    drug_results.append(drug_result)
                    
            # Phase 5: Generate comprehensive report
            logger.info("Phase 5: Generating comprehensive report...")
            report_result = await self.report_generator.generate_report(
                sample_id, vcf_results, annotation_results, structure_results, drug_results
            )
            
            # Compile final results
            final_results = {
                "pipeline_status": "completed",
                "sample_id": sample_id,
                "vcf_file": vcf_file,
                "output_directory": output_dir,
                "vcf_results": vcf_results,
                "annotation_results": annotation_results,
                "structure_results": structure_results,
                "drug_results": drug_results,
                "report_result": report_result,
                "summary": {
                    "total_variants": len(variants),
                    "coding_variants": len([v for v in annotation_results if v.get("vep_annotation", {}).get("coding", False)]),
                    "structure_analyses": len(structure_results),
                    "drug_matches": len([d for d in drug_results if d.get("matched_drugs")]),
                    "report_files": report_result.get("report_files", {})
                }
            }
            
            logger.info(f"Pipeline completed successfully for {sample_id}")
            return final_results
            
        except Exception as e:
            logger.error(f"Pipeline failed for {sample_id}: {str(e)}")
            return {
                "pipeline_status": "failed",
                "sample_id": sample_id,
                "vcf_file": vcf_file,
                "error": str(e),
                "summary": {}
            }

# Create the MCP server instance
app = Server("pharmacogenomics-pipeline")
pipeline = None

@app.list_tools()
async def handle_list_tools() -> List[Tool]:
    """List all available tools in the pharmacogenomics pipeline"""
    return [
        Tool(
            name="parse_vcf",
            description="Parse VCF file and perform quality control analysis",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the VCF file to analyze"
                    },
                    "sample_id": {
                        "type": "string", 
                        "description": "Optional sample identifier"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="annotate_variants",
            description="Annotate variants using Ensembl VEP for gene and consequence information",
            inputSchema={
                "type": "object",
                "properties": {
                    "variants": {
                        "type": "array",
                        "description": "Array of variant objects from parse_vcf"
                    }
                },
                "required": ["variants"]
            }
        ),
        Tool(
            name="get_protein_structure",
            description="Get protein structure information and highlight variant position",
            inputSchema={
                "type": "object",
                "properties": {
                    "uniprot_id": {
                        "type": "string",
                        "description": "UniProt accession ID"
                    },
                    "residue_pos": {
                        "type": "integer",
                        "description": "Amino acid position of the variant"
                    },
                    "amino_acid_change": {
                        "type": "string",
                        "description": "HGVSp notation of the change (e.g., p.Arg144Cys)"
                    }
                },
                "required": ["uniprot_id", "residue_pos"]
            }
        ),
        Tool(
            name="match_drugs",
            description="Find drugs that target the mutated gene/protein",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "Gene symbol (e.g., CYP2D6)"
                    },
                    "uniprot_id": {
                        "type": "string",
                        "description": "UniProt ID for additional matching"
                    },
                    "hgvsp": {
                        "type": "string",
                        "description": "HGVSp notation for variant-specific lookups"
                    }
                },
                "required": ["gene_symbol"]
            }
        ),
        Tool(
            name="generate_report",
            description="Generate comprehensive pharmacogenomics report",
            inputSchema={
                "type": "object",
                "properties": {
                    "sample_id": {
                        "type": "string",
                        "description": "Sample identifier"
                    },
                    "vcf_results": {
                        "type": "object",
                        "description": "Results from parse_vcf"
                    },
                    "annotation_results": {
                        "type": "array",
                        "description": "Results from annotate_variants"
                    },
                    "structure_results": {
                        "type": "array", 
                        "description": "Results from get_protein_structure calls"
                    },
                    "drug_results": {
                        "type": "array",
                        "description": "Results from match_drugs calls"
                    },
                    "output_format": {
                        "type": "string",
                        "enum": ["json", "html", "both"],
                        "default": "both"
                    }
                },
                "required": ["sample_id"]
            }
        ),
        Tool(
            name="run_full_pipeline",
            description="Run the complete pharmacogenomics analysis pipeline",
            inputSchema={
                "type": "object",
                "properties": {
                    "vcf_file": {
                        "type": "string",
                        "description": "Path to VCF file"
                    },
                    "sample_id": {
                        "type": "string",
                        "description": "Sample identifier"
                    },
                    "output_dir": {
                        "type": "string",
                        "description": "Output directory for reports",
                        "default": "results"
                    }
                },
                "required": ["vcf_file", "sample_id"]
            }
        )
    ]

@app.call_tool()
async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> List[types.TextContent | types.ImageContent]:
    """Handle tool execution"""
    global pipeline
    
    if pipeline is None:
        pipeline = PharmacogenomicsMCP()
    
    try:
        if name == "parse_vcf":
            result = await pipeline.vcf_parser.parse_vcf(
                arguments["file_path"],
                arguments.get("sample_id")
            )
            
        elif name == "annotate_variants":
            result = await pipeline.vep_annotator.annotate_variants(
                arguments["variants"]
            )
            
        elif name == "get_protein_structure":
            result = await pipeline.structure_analyzer.analyze_structure(
                arguments["uniprot_id"],
                arguments["residue_pos"],
                arguments.get("amino_acid_change")
            )
            
        elif name == "match_drugs":
            result = await pipeline.drug_matcher.match_drugs(
                arguments["gene_symbol"],
                arguments.get("uniprot_id"),
                arguments.get("hgvsp")
            )
            
        elif name == "generate_report":
            result = await pipeline.report_generator.generate_report(
                arguments["sample_id"],
                arguments.get("vcf_results"),
                arguments.get("annotation_results"),
                arguments.get("structure_results"),
                arguments.get("drug_results"),
                arguments.get("output_format", "both")
            )
            
        elif name == "run_full_pipeline":
            result = await pipeline.run_full_pipeline(
                arguments["vcf_file"],
                arguments["sample_id"],
                arguments.get("output_dir", "results")
            )
            
        else:
            raise ValueError(f"Unknown tool: {name}")
            
        return [types.TextContent(
            type="text",
            text=json.dumps(result, indent=2, default=str)
        )]
        
    except Exception as e:
        logger.error(f"Error executing tool {name}: {str(e)}")
        return [types.TextContent(
            type="text", 
            text=f"Error: {str(e)}"
        )]

async def main():
    """Run the MCP server"""
    async with stdio_server() as (read_stream, write_stream):
        await app.run(read_stream, write_stream, app.create_initialization_options())

if __name__ == "__main__":
    asyncio.run(main())