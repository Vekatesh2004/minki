"""
Report Generator Module - Phase 5
Handles comprehensive pharmacogenomics report generation
"""

import asyncio
import logging
from typing import Any, Dict, List, Optional
import json
from pathlib import Path
from datetime import datetime
import base64

logger = logging.getLogger(__name__)

class ReportGenerator:
    def __init__(self, config: Dict[str, Any]):
        """Initialize report generator with configuration"""
        self.config = config
        self.output_config = config.get("output", {})
        self.formats = self.output_config.get("formats", ["json", "html"])
        self.include_structures = self.output_config.get("include_structures", True)
        
    async def generate_report(self, sample_id: str, 
                             vcf_results: Optional[Dict[str, Any]] = None,
                             annotation_results: Optional[List[Dict[str, Any]]] = None,
                             structure_results: Optional[List[Dict[str, Any]]] = None,
                             drug_results: Optional[List[Dict[str, Any]]] = None,
                             output_format: str = "both") -> Dict[str, Any]:
        """
        Generate comprehensive pharmacogenomics report
        
        Args:
            sample_id: Sample identifier
            vcf_results: Results from VCF parsing
            annotation_results: Results from variant annotation
            structure_results: Results from protein structure analysis
            drug_results: Results from drug matching
            output_format: Output format ("json", "html", or "both")
            
        Returns:
            Dictionary containing report paths and summary
        """
        try:
            # Compile all data
            report_data = await self._compile_report_data(
                sample_id, vcf_results, annotation_results, 
                structure_results, drug_results
            )
            
            # Generate reports in requested formats
            report_files = {}
            
            if output_format in ["json", "both"]:
                json_path = await self._generate_json_report(report_data, sample_id)
                report_files["json"] = json_path
                
            if output_format in ["html", "both"]:
                html_path = await self._generate_html_report(report_data, sample_id)
                report_files["html"] = html_path
                
            # Generate summary
            summary = self._generate_summary(report_data)
            
            return {
                "sample_id": sample_id,
                "report_files": report_files,
                "summary": summary,
                "generated_at": datetime.now().isoformat()
            }
            
        except Exception as e:
            logger.error(f"Error generating report for {sample_id}: {str(e)}")
            return {
                "sample_id": sample_id,
                "error": str(e),
                "report_files": {},
                "summary": {},
                "generated_at": datetime.now().isoformat()
            }
            
    async def _compile_report_data(self, sample_id: str,
                                  vcf_results: Optional[Dict[str, Any]],
                                  annotation_results: Optional[List[Dict[str, Any]]],
                                  structure_results: Optional[List[Dict[str, Any]]],
                                  drug_results: Optional[List[Dict[str, Any]]]) -> Dict[str, Any]:
        """Compile all analysis results into structured report data"""
        
        # Initialize report structure
        report_data = {
            "sample_information": {
                "sample_id": sample_id,
                "analysis_date": datetime.now().isoformat(),
                "pipeline_version": "1.0.0"
            },
            "quality_control": {},
            "variants": [],
            "pharmacogenomic_summary": {
                "total_variants": 0,
                "coding_variants": 0,
                "high_confidence_variants": 0,
                "variants_with_drug_interactions": 0,
                "genes_with_variants": set(),
                "drugs_affected": set()
            },
            "detailed_results": []
        }
        
        # Add QC information
        if vcf_results:
            report_data["quality_control"] = vcf_results.get("qc_summary", {})
            report_data["variants"] = vcf_results.get("variants", [])
            report_data["pharmacogenomic_summary"]["total_variants"] = len(report_data["variants"])
            
        # Process annotated variants
        if annotation_results:
            report_data["detailed_results"] = await self._process_annotated_variants(
                annotation_results, structure_results, drug_results, report_data
            )
            
        # Convert sets to lists for JSON serialization
        report_data["pharmacogenomic_summary"]["genes_with_variants"] = list(report_data["pharmacogenomic_summary"]["genes_with_variants"])
        report_data["pharmacogenomic_summary"]["drugs_affected"] = list(report_data["pharmacogenomic_summary"]["drugs_affected"])
        
        return report_data
        
    async def _process_annotated_variants(self, annotation_results: List[Dict[str, Any]],
                                         structure_results: Optional[List[Dict[str, Any]]],
                                         drug_results: Optional[List[Dict[str, Any]]],
                                         report_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Process annotated variants and merge with structure and drug data"""
        
        detailed_results = []
        
        # Create lookup maps for structure and drug results
        structure_map = {}
        if structure_results:
            for result in structure_results:
                uniprot_id = result.get("uniprot_id")
                if uniprot_id:
                    structure_map[uniprot_id] = result
                    
        drug_map = {}
        if drug_results:
            for result in drug_results:
                gene_symbol = result.get("gene_symbol")
                if gene_symbol:
                    drug_map[gene_symbol] = result
                    
        # Process each annotated variant
        for variant in annotation_results:
            vep_annotation = variant.get("vep_annotation", {})
            
            # Skip non-coding variants for detailed analysis
            if not vep_annotation.get("coding", False):
                continue
                
            report_data["pharmacogenomic_summary"]["coding_variants"] += 1
            
            if variant.get("confidence") == "high":
                report_data["pharmacogenomic_summary"]["high_confidence_variants"] += 1
                
            gene_symbol = vep_annotation.get("gene_symbol")
            if gene_symbol:
                report_data["pharmacogenomic_summary"]["genes_with_variants"].add(gene_symbol)
                
            # Create detailed variant result
            variant_result = {
                "variant_info": {
                    "position": f"{variant.get('chrom')}:{variant.get('pos')}",
                    "reference": variant.get("ref"),
                    "alternate": variant.get("alt"),
                    "quality": variant.get("qual"),
                    "confidence": variant.get("confidence")
                },
                "gene_annotation": {
                    "gene_symbol": gene_symbol,
                    "gene_id": vep_annotation.get("gene_id"),
                    "transcript_id": vep_annotation.get("transcript_id"),
                    "consequence_terms": vep_annotation.get("consequence_terms", []),
                    "hgvsc": vep_annotation.get("hgvsc"),
                    "hgvsp": vep_annotation.get("hgvsp"),
                    "canonical": vep_annotation.get("canonical", False)
                },
                "protein_analysis": None,
                "drug_interactions": None,
                "clinical_significance": "unknown"
            }
            
            # Add structure analysis if available
            uniprot_ids = vep_annotation.get("uniprot_ids", [])
            for uniprot_id in uniprot_ids:
                if uniprot_id in structure_map:
                    structure_data = structure_map[uniprot_id]
                    variant_result["protein_analysis"] = {
                        "uniprot_id": uniprot_id,
                        "position_analysis": structure_data.get("position_analysis"),
                        "structural_impact": structure_data.get("position_analysis", {}).get("structural_impact_prediction"),
                        "confidence_score": structure_data.get("position_analysis", {}).get("plddt_score"),
                        "in_domain": structure_data.get("position_analysis", {}).get("in_domain", False),
                        "in_binding_site": structure_data.get("position_analysis", {}).get("in_binding_site", False),
                        "visualization_available": structure_data.get("visualization", {}).get("available", False)
                    }
                    break
                    
            # Add drug interaction data if available
            if gene_symbol and gene_symbol in drug_map:
                drug_data = drug_map[gene_symbol]
                matched_drugs = drug_data.get("matched_drugs", [])
                
                if matched_drugs:
                    report_data["pharmacogenomic_summary"]["variants_with_drug_interactions"] += 1
                    
                    # Extract drug names for summary
                    for drug in matched_drugs:
                        drug_name = drug.get("drug_name")
                        if drug_name:
                            report_data["pharmacogenomic_summary"]["drugs_affected"].add(drug_name)
                            
                    variant_result["drug_interactions"] = {
                        "total_matches": len(matched_drugs),
                        "high_evidence_drugs": [d for d in matched_drugs if d.get("evidence_level") in ["level_1", "level_2", "variant_specific"]],
                        "all_matches": matched_drugs,
                        "pharmgkb_annotations": drug_data.get("pharmgkb_annotations", [])
                    }
                    
                    # Determine clinical significance based on evidence
                    high_evidence_count = len(variant_result["drug_interactions"]["high_evidence_drugs"])
                    if high_evidence_count > 0:
                        variant_result["clinical_significance"] = "high"
                    elif len(matched_drugs) > 0:
                        variant_result["clinical_significance"] = "moderate"
                    else:
                        variant_result["clinical_significance"] = "low"
                        
            detailed_results.append(variant_result)
            
        return detailed_results
        
    async def _generate_json_report(self, report_data: Dict[str, Any], sample_id: str) -> str:
        """Generate JSON report"""
        
        try:
            # Ensure output directory exists
            output_dir = Path("results")
            output_dir.mkdir(exist_ok=True)
            
            # Generate filename
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            json_path = output_dir / f"{sample_id}_pharmacogenomics_{timestamp}.json"
            
            # Write JSON report
            with open(json_path, 'w') as f:
                json.dump(report_data, f, indent=2, default=str)
                
            logger.info(f"Generated JSON report: {json_path}")
            return str(json_path)
            
        except Exception as e:
            logger.error(f"Error generating JSON report: {str(e)}")
            raise
            
    async def _generate_html_report(self, report_data: Dict[str, Any], sample_id: str) -> str:
        """Generate HTML report"""
        
        try:
            # Ensure output directory exists
            output_dir = Path("results")
            output_dir.mkdir(exist_ok=True)
            
            # Generate filename
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            html_path = output_dir / f"{sample_id}_pharmacogenomics_{timestamp}.html"
            
            # Generate HTML content
            html_content = self._create_html_report(report_data)
            
            # Write HTML report
            with open(html_path, 'w') as f:
                f.write(html_content)
                
            logger.info(f"Generated HTML report: {html_path}")
            return str(html_path)
            
        except Exception as e:
            logger.error(f"Error generating HTML report: {str(e)}")
            raise
            
    def _create_html_report(self, report_data: Dict[str, Any]) -> str:
        """Create HTML report content"""
        
        sample_info = report_data.get("sample_information", {})
        qc_data = report_data.get("quality_control", {})
        summary = report_data.get("pharmacogenomic_summary", {})
        detailed_results = report_data.get("detailed_results", [])
        
        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pharmacogenomics Report - {sample_info.get('sample_id', 'Unknown')}</title>
    <style>
        body {{ 
            font-family: Arial, sans-serif; 
            line-height: 1.6; 
            margin: 0; 
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            border-bottom: 2px solid #3498db;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        .section {{
            margin-bottom: 30px;
        }}
        .section h2 {{
            color: #2c3e50;
            border-left: 4px solid #3498db;
            padding-left: 15px;
        }}
        .qc-summary, .pgx-summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .metric-card {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2980b9;
        }}
        .metric-label {{
            font-size: 0.9em;
            color: #7f8c8d;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            border: 1px solid #bdc3c7;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #34495e;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        .confidence-high {{ color: #27ae60; font-weight: bold; }}
        .confidence-moderate {{ color: #f39c12; font-weight: bold; }}
        .confidence-low {{ color: #e74c3c; font-weight: bold; }}
        .evidence-level-1, .evidence-level-2 {{ background-color: #d5f4e6; }}
        .evidence-level-3, .evidence-level-4 {{ background-color: #fff3cd; }}
        .no-evidence {{ background-color: #f8d7da; }}
        .collapsible {{
            background-color: #f1f1f1;
            color: #444;
            cursor: pointer;
            padding: 10px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 15px;
        }}
        .collapsible:hover {{
            background-color: #ddd;
        }}
        .content {{
            padding: 0 18px;
            display: none;
            overflow: hidden;
            background-color: #f9f9f9;
            margin-bottom: 10px;
        }}
    </style>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Pharmacogenomics Analysis Report</h1>
            <h2>Sample: {sample_info.get('sample_id', 'Unknown')}</h2>
            <p>Generated: {sample_info.get('analysis_date', 'Unknown')}</p>
        </div>

        <div class="section">
            <h2>Quality Control Summary</h2>
            <div class="qc-summary">
                <div class="metric-card">
                    <div class="metric-value">{qc_data.get('total_variants', 0)}</div>
                    <div class="metric-label">Total Variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{qc_data.get('pass_rate', 0):.1%}</div>
                    <div class="metric-label">Pass Rate</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{qc_data.get('high_confidence_rate', 0):.1%}</div>
                    <div class="metric-label">High Confidence</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{qc_data.get('ts_tv_ratio', 0):.2f}</div>
                    <div class="metric-label">Ts/Tv Ratio</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Pharmacogenomic Summary</h2>
            <div class="pgx-summary">
                <div class="metric-card">
                    <div class="metric-value">{summary.get('coding_variants', 0)}</div>
                    <div class="metric-label">Coding Variants</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{summary.get('variants_with_drug_interactions', 0)}</div>
                    <div class="metric-label">Variants with Drug Interactions</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(summary.get('genes_with_variants', []))}</div>
                    <div class="metric-label">Genes Affected</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(summary.get('drugs_affected', []))}</div>
                    <div class="metric-label">Drugs Potentially Affected</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Detailed Variant Analysis</h2>
            {self._generate_variant_table_html(detailed_results)}
        </div>

        <div class="section">
            <h2>Genes with Pharmacogenomic Variants</h2>
            <p><strong>Genes:</strong> {', '.join(summary.get('genes_with_variants', []))}</p>
        </div>

        <div class="section">
            <h2>Potentially Affected Drugs</h2>
            <p><strong>Drugs:</strong> {', '.join(summary.get('drugs_affected', []))}</p>
        </div>

        <div class="section">
            <h2>Methodology</h2>
            <p>This analysis used the following tools and databases:</p>
            <ul>
                <li><strong>Variant Annotation:</strong> Ensembl VEP</li>
                <li><strong>Protein Structure:</strong> UniProt and AlphaFold Database</li>
                <li><strong>Drug Interactions:</strong> DrugBank and PharmGKB</li>
                <li><strong>Quality Control:</strong> Custom metrics and filtering</li>
            </ul>
        </div>

        <div class="section">
            <h2>Disclaimer</h2>
            <p style="background-color: #fff3cd; padding: 15px; border-radius: 5px; border-left: 4px solid #ffc107;">
                <strong>Important:</strong> This report is for research purposes only and should not be used for clinical decision-making without proper validation and consultation with healthcare professionals. Pharmacogenomic predictions are based on current knowledge and may not account for all factors affecting drug response.
            </p>
        </div>
    </div>

    <script>
        // Make collapsible content work
        var coll = document.getElementsByClassName("collapsible");
        for (var i = 0; i < coll.length; i++) {{
            coll[i].addEventListener("click", function() {{
                this.classList.toggle("active");
                var content = this.nextElementSibling;
                if (content.style.display === "block") {{
                    content.style.display = "none";
                }} else {{
                    content.style.display = "block";
                }}
            }});
        }}
    </script>
</body>
</html>
        """
        
        return html
        
    def _generate_variant_table_html(self, detailed_results: List[Dict[str, Any]]) -> str:
        """Generate HTML table for detailed variant results"""
        
        if not detailed_results:
            return "<p>No coding variants with drug interactions found.</p>"
            
        table_html = """
        <table>
            <thead>
                <tr>
                    <th>Position</th>
                    <th>Gene</th>
                    <th>HGVSp</th>
                    <th>Consequence</th>
                    <th>Confidence</th>
                    <th>Drug Interactions</th>
                    <th>Clinical Significance</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for result in detailed_results:
            variant_info = result.get("variant_info", {})
            gene_annotation = result.get("gene_annotation", {})
            drug_interactions = result.get("drug_interactions", {})
            
            position = variant_info.get("position", "Unknown")
            gene_symbol = gene_annotation.get("gene_symbol", "Unknown")
            hgvsp = gene_annotation.get("hgvsp", "N/A")
            consequence = ", ".join(gene_annotation.get("consequence_terms", []))
            confidence = variant_info.get("confidence", "unknown")
            clinical_sig = result.get("clinical_significance", "unknown")
            
            # Format confidence
            confidence_class = f"confidence-{confidence}" if confidence in ["high", "moderate", "low"] else ""
            
            # Format drug interactions
            drug_info = "None"
            if drug_interactions:
                total_matches = drug_interactions.get("total_matches", 0)
                high_evidence = len(drug_interactions.get("high_evidence_drugs", []))
                if high_evidence > 0:
                    drug_info = f"{total_matches} drugs ({high_evidence} high evidence)"
                else:
                    drug_info = f"{total_matches} drugs"
                    
            # Format clinical significance
            sig_class = f"confidence-{clinical_sig}" if clinical_sig in ["high", "moderate", "low"] else ""
            
            table_html += f"""
                <tr>
                    <td>{position}</td>
                    <td>{gene_symbol}</td>
                    <td>{hgvsp}</td>
                    <td>{consequence}</td>
                    <td class="{confidence_class}">{confidence.title()}</td>
                    <td>{drug_info}</td>
                    <td class="{sig_class}">{clinical_sig.title()}</td>
                </tr>
            """
            
        table_html += """
            </tbody>
        </table>
        """
        
        return table_html
        
    def _generate_summary(self, report_data: Dict[str, Any]) -> Dict[str, Any]:
        """Generate a concise summary of the analysis"""
        
        summary = report_data.get("pharmacogenomic_summary", {})
        detailed_results = report_data.get("detailed_results", [])
        
        # Count variants by clinical significance
        sig_counts = {"high": 0, "moderate": 0, "low": 0, "unknown": 0}
        for result in detailed_results:
            sig = result.get("clinical_significance", "unknown")
            sig_counts[sig] = sig_counts.get(sig, 0) + 1
            
        # Find top genes and drugs
        gene_drug_map = {}
        for result in detailed_results:
            gene = result.get("gene_annotation", {}).get("gene_symbol")
            drug_interactions = result.get("drug_interactions", {})
            if gene and drug_interactions:
                matched_drugs = drug_interactions.get("matched_drugs", [])
                gene_drug_map[gene] = len(matched_drugs)
                
        top_genes = sorted(gene_drug_map.items(), key=lambda x: x[1], reverse=True)[:5]
        
        return {
            "total_variants_analyzed": summary.get("total_variants", 0),
            "coding_variants": summary.get("coding_variants", 0),
            "pharmacogenomically_relevant": summary.get("variants_with_drug_interactions", 0),
            "clinical_significance_breakdown": sig_counts,
            "genes_analyzed": len(summary.get("genes_with_variants", [])),
            "drugs_potentially_affected": len(summary.get("drugs_affected", [])),
            "top_pharmacogenes": dict(top_genes),
            "high_priority_variants": sig_counts["high"],
            "recommendation": self._generate_recommendation(sig_counts, summary)
        }
        
    def _generate_recommendation(self, sig_counts: Dict[str, int], summary: Dict[str, Any]) -> str:
        """Generate clinical recommendations based on findings"""
        
        high_sig = sig_counts.get("high", 0)
        moderate_sig = sig_counts.get("moderate", 0)
        
        if high_sig > 0:
            return f"High priority: {high_sig} variant(s) with strong pharmacogenomic evidence found. Consider genetic testing confirmation and potential drug therapy adjustments."
        elif moderate_sig > 0:
            return f"Moderate priority: {moderate_sig} variant(s) with potential drug interactions identified. Review with clinical pharmacist or geneticist."
        elif summary.get("variants_with_drug_interactions", 0) > 0:
            return "Low priority: Some variants with potential drug interactions found, but evidence is limited. Monitor for unusual drug responses."
        else:
            return "No significant pharmacogenomic variants identified in this analysis."