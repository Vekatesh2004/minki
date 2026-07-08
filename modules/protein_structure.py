"""
Protein Structure Analyzer Module - Phase 3
Handles UniProt mapping and AlphaFold structure analysis
"""

import asyncio
import logging
from typing import Any, Dict, List, Optional, Tuple
import aiohttp
import json
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

class ProteinStructureAnalyzer:
    def __init__(self, config: Dict[str, Any]):
        """Initialize protein structure analyzer with configuration"""
        self.config = config
        self.uniprot_config = config.get("uniprot", {})
        self.alphafold_config = config.get("alphafold", {})
        
        self.uniprot_base_url = self.uniprot_config.get("base_url", "https://rest.uniprot.org")
        self.uniprot_mapping_url = self.uniprot_config.get("id_mapping_url", "https://rest.uniprot.org/idmapping")
        self.alphafold_base_url = self.alphafold_config.get("base_url", "https://alphafold.ebi.ac.uk/files")
        
    async def analyze_structure(self, uniprot_id: str, residue_pos: int, 
                               amino_acid_change: Optional[str] = None) -> Dict[str, Any]:
        """
        Analyze protein structure for a given variant
        
        Args:
            uniprot_id: UniProt accession ID
            residue_pos: Amino acid position of the variant
            amino_acid_change: HGVSp notation of the change
            
        Returns:
            Dictionary containing structure analysis results
        """
        try:
            # Get UniProt entry information
            uniprot_data = await self._get_uniprot_entry(uniprot_id)
            
            # Get AlphaFold structure
            alphafold_data = await self._get_alphafold_structure(uniprot_id)
            
            # Analyze variant position
            position_analysis = await self._analyze_variant_position(
                uniprot_data, alphafold_data, residue_pos, amino_acid_change
            )
            
            # Generate structure visualization
            visualization = await self._generate_structure_visualization(
                uniprot_id, residue_pos, amino_acid_change
            )
            
            return {
                "uniprot_id": uniprot_id,
                "residue_position": residue_pos,
                "amino_acid_change": amino_acid_change,
                "uniprot_data": uniprot_data,
                "alphafold_data": alphafold_data,
                "position_analysis": position_analysis,
                "visualization": visualization
            }
            
        except Exception as e:
            logger.error(f"Error analyzing structure for {uniprot_id} at position {residue_pos}: {str(e)}")
            return {
                "uniprot_id": uniprot_id,
                "residue_position": residue_pos,
                "amino_acid_change": amino_acid_change,
                "error": str(e),
                "uniprot_data": None,
                "alphafold_data": None,
                "position_analysis": None,
                "visualization": None
            }
            
    async def map_gene_to_uniprot(self, gene_symbol: str, ensembl_id: Optional[str] = None) -> List[str]:
        """
        Map gene symbol or Ensembl ID to UniProt accessions
        
        Args:
            gene_symbol: Gene symbol (e.g., CYP2D6)
            ensembl_id: Optional Ensembl gene ID
            
        Returns:
            List of UniProt accession IDs
        """
        try:
            uniprot_ids = []
            
            # Try mapping by gene name first
            if gene_symbol:
                ids_from_gene = await self._map_by_gene_name(gene_symbol)
                uniprot_ids.extend(ids_from_gene)
                
            # Try mapping by Ensembl ID
            if ensembl_id:
                ids_from_ensembl = await self._map_by_ensembl_id(ensembl_id)
                uniprot_ids.extend(ids_from_ensembl)
                
            # Remove duplicates while preserving order
            return list(dict.fromkeys(uniprot_ids))
            
        except Exception as e:
            logger.error(f"Error mapping gene {gene_symbol} to UniProt: {str(e)}")
            return []
            
    async def _get_uniprot_entry(self, uniprot_id: str) -> Dict[str, Any]:
        """Get UniProt entry data"""
        
        url = f"{self.uniprot_base_url}/uniprotkb/{uniprot_id}.json"
        
        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=30)) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        return self._parse_uniprot_entry(data)
                    else:
                        logger.warning(f"UniProt entry not found: {uniprot_id}")
                        return {}
                        
        except Exception as e:
            logger.error(f"Error fetching UniProt entry {uniprot_id}: {str(e)}")
            return {}
            
    def _parse_uniprot_entry(self, uniprot_data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse UniProt entry to extract relevant information"""
        
        parsed = {
            "accession": uniprot_data.get("primaryAccession"),
            "protein_name": None,
            "gene_names": [],
            "organism": None,
            "sequence_length": 0,
            "sequence": None,
            "domains": [],
            "binding_sites": [],
            "active_sites": [],
            "other_features": []
        }
        
        # Protein names
        protein_existence = uniprot_data.get("proteinExistence", {})
        if "recommendedName" in uniprot_data.get("proteinDescription", {}):
            parsed["protein_name"] = uniprot_data["proteinDescription"]["recommendedName"].get("fullName", {}).get("value")
        
        # Gene names
        gene_names = uniprot_data.get("genes", [])
        for gene in gene_names:
            if "geneName" in gene:
                parsed["gene_names"].append(gene["geneName"]["value"])
                
        # Organism
        organism = uniprot_data.get("organism", {})
        if "scientificName" in organism:
            parsed["organism"] = organism["scientificName"]
            
        # Sequence information
        sequence = uniprot_data.get("sequence", {})
        parsed["sequence_length"] = sequence.get("length", 0)
        parsed["sequence"] = sequence.get("value", "")
        
        # Features (domains, binding sites, etc.)
        features = uniprot_data.get("features", [])
        for feature in features:
            feature_type = feature.get("type")
            description = feature.get("description", "")
            
            # Get position information
            location = feature.get("location", {})
            start_pos = None
            end_pos = None
            
            if "start" in location:
                start_pos = location["start"].get("value")
            if "end" in location:
                end_pos = location["end"].get("value")
                
            feature_info = {
                "type": feature_type,
                "description": description,
                "start": start_pos,
                "end": end_pos
            }
            
            # Categorize features
            if feature_type == "Domain":
                parsed["domains"].append(feature_info)
            elif feature_type in ["Binding site", "Metal binding"]:
                parsed["binding_sites"].append(feature_info)
            elif feature_type == "Active site":
                parsed["active_sites"].append(feature_info)
            else:
                parsed["other_features"].append(feature_info)
                
        return parsed
        
    async def _get_alphafold_structure(self, uniprot_id: str) -> Dict[str, Any]:
        """Get AlphaFold structure data"""

        alphafold_data = {
            "uniprot_id": uniprot_id,
            "pdb_available": False,
            "confidence_available": False,
            "pdb_path": None,
            "confidence_scores": {},
            "model_quality": None
        }

        # AlphaFold periodically bumps its model version (v4, v6, ...). Rather
        # than hardcode a version in the URL, query the prediction API to get
        # the current PDB URL for this accession.
        pdb_url = None
        try:
            api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=30)) as session:
                async with session.get(api_url) as api_resp:
                    if api_resp.status == 200:
                        api_data = await api_resp.json()
                        if api_data:
                            pdb_url = api_data[0].get("pdbUrl")
                            alphafold_data["model_quality"] = api_data[0].get("globalMetricValue")
        except Exception as e:
            logger.warning(f"AlphaFold API lookup failed for {uniprot_id}: {str(e)}")

        # Fallback to legacy versioned URL if the API gave us nothing
        if not pdb_url:
            pdb_url = f"{self.alphafold_base_url}/AF-{uniprot_id}-F1-model_v4.pdb"

        try:
            # Download PDB file
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=60)) as session:
                async with session.get(pdb_url) as response:
                    if response.status == 200:
                        pdb_content = await response.text()
                        
                        # Save to temporary file
                        temp_pdb = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
                        temp_pdb.write(pdb_content)
                        temp_pdb.close()
                        
                        alphafold_data["pdb_available"] = True
                        alphafold_data["pdb_path"] = temp_pdb.name
                        
                        # Parse B-factors (pLDDT scores) from PDB
                        confidence_scores = self._parse_plddt_from_pdb(pdb_content)
                        alphafold_data["confidence_scores"] = confidence_scores
                        alphafold_data["confidence_available"] = True
                        
                    else:
                        logger.warning(f"AlphaFold structure not available for {uniprot_id}")
                        
        except Exception as e:
            logger.error(f"Error downloading AlphaFold structure for {uniprot_id}: {str(e)}")
            
        return alphafold_data
        
    def _parse_plddt_from_pdb(self, pdb_content: str) -> Dict[int, float]:
        """Parse pLDDT confidence scores from AlphaFold PDB file"""
        
        confidence_scores = {}
        
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':  # CA atoms only
                try:
                    residue_num = int(line[22:26].strip())
                    plddt_score = float(line[60:66].strip())
                    confidence_scores[residue_num] = plddt_score
                except (ValueError, IndexError):
                    continue
                    
        return confidence_scores
        
    async def _analyze_variant_position(self, uniprot_data: Dict[str, Any], 
                                       alphafold_data: Dict[str, Any], 
                                       residue_pos: int, 
                                       amino_acid_change: Optional[str]) -> Dict[str, Any]:
        """Analyze the variant position in the protein structure context"""
        
        analysis = {
            "residue_position": residue_pos,
            "amino_acid_change": amino_acid_change,
            "in_domain": False,
            "domain_info": [],
            "in_binding_site": False,
            "binding_site_info": [],
            "in_active_site": False,
            "active_site_info": [],
            "plddt_score": None,
            "confidence_level": None,
            "structural_impact_prediction": None
        }
        
        # Check if position falls within annotated domains
        for domain in uniprot_data.get("domains", []):
            start = domain.get("start")
            end = domain.get("end")
            if start and end and start <= residue_pos <= end:
                analysis["in_domain"] = True
                analysis["domain_info"].append(domain)
                
        # Check binding sites
        for binding_site in uniprot_data.get("binding_sites", []):
            start = binding_site.get("start")
            end = binding_site.get("end")
            if start and end and start <= residue_pos <= end:
                analysis["in_binding_site"] = True
                analysis["binding_site_info"].append(binding_site)
                
        # Check active sites
        for active_site in uniprot_data.get("active_sites", []):
            start = active_site.get("start")
            end = active_site.get("end")
            if start and end and start <= residue_pos <= end:
                analysis["in_active_site"] = True
                analysis["active_site_info"].append(active_site)
                
        # Get AlphaFold confidence score
        confidence_scores = alphafold_data.get("confidence_scores", {})
        if residue_pos in confidence_scores:
            analysis["plddt_score"] = confidence_scores[residue_pos]
            analysis["confidence_level"] = self._classify_confidence_level(confidence_scores[residue_pos])
            
        # Predict structural impact
        analysis["structural_impact_prediction"] = self._predict_structural_impact(analysis, amino_acid_change)
        
        return analysis
        
    def _classify_confidence_level(self, plddt_score: float) -> str:
        """Classify AlphaFold confidence level based on pLDDT score"""
        
        if plddt_score >= 90:
            return "very_high"
        elif plddt_score >= 70:
            return "confident"
        elif plddt_score >= 50:
            return "low"
        else:
            return "very_low"
            
    def _predict_structural_impact(self, analysis: Dict[str, Any], amino_acid_change: Optional[str]) -> str:
        """Predict structural impact of the variant"""
        
        impact_factors = []
        
        # Domain involvement
        if analysis["in_domain"]:
            impact_factors.append("domain")
            
        # Binding site involvement
        if analysis["in_binding_site"]:
            impact_factors.append("binding_site")
            
        # Active site involvement  
        if analysis["in_active_site"]:
            impact_factors.append("active_site")
            
        # Confidence level
        confidence = analysis.get("confidence_level")
        if confidence in ["very_high", "confident"]:
            impact_factors.append("high_confidence_structure")
            
        # Amino acid change type
        if amino_acid_change:
            if "nonsense" in amino_acid_change.lower() or "*" in amino_acid_change:
                return "high_impact"
            elif "synonymous" in amino_acid_change.lower():
                return "low_impact"
                
        # Overall prediction
        if "active_site" in impact_factors:
            return "high_impact"
        elif "binding_site" in impact_factors:
            return "moderate_impact"
        elif "domain" in impact_factors:
            return "moderate_impact"
        elif confidence in ["very_low", "low"]:
            return "uncertain_impact"
        else:
            return "low_impact"
            
    async def _generate_structure_visualization(self, uniprot_id: str, residue_pos: int, 
                                              amino_acid_change: Optional[str]) -> Dict[str, Any]:
        """Generate structure visualization highlighting the variant position"""
        
        visualization = {
            "py3dmol_script": None,
            "html_viewer": None,
            "image_path": None,
            "available": False
        }
        
        try:
            # Generate py3Dmol visualization script
            py3dmol_script = self._generate_py3dmol_script(uniprot_id, residue_pos, amino_acid_change)
            visualization["py3dmol_script"] = py3dmol_script
            
            # Generate HTML viewer
            html_viewer = self._generate_html_viewer(uniprot_id, residue_pos, amino_acid_change)
            visualization["html_viewer"] = html_viewer
            
            visualization["available"] = True
            
        except Exception as e:
            logger.error(f"Error generating structure visualization: {str(e)}")
            
        return visualization
        
    def _generate_py3dmol_script(self, uniprot_id: str, residue_pos: int, 
                                amino_acid_change: Optional[str]) -> str:
        """Generate py3Dmol JavaScript for structure visualization"""
        
        script = f"""
        // Load AlphaFold structure for {uniprot_id}
        var viewer = $3Dmol.createViewer('viewer', {{width: 800, height: 600}});
        
        // Load PDB structure
        $3Dmol.download('pdb:AF-{uniprot_id}-F1-model_v4', viewer, {{}}, function() {{
            
            // Style the protein
            viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
            
            // Highlight the variant position in red
            viewer.addStyle({{resi: {residue_pos}}}, {{sphere: {{color: 'red', radius: 2.0}}}});
            
            // Add label for the variant
            viewer.addLabel('{amino_acid_change or f"Position {residue_pos}"}', 
                          {{position: {{resi: {residue_pos}}}, 
                            backgroundColor: 'black', 
                            fontColor: 'white'}});
            
            // Set camera and render
            viewer.zoomTo();
            viewer.render();
        }});
        """
        
        return script
        
    def _generate_html_viewer(self, uniprot_id: str, residue_pos: int, 
                             amino_acid_change: Optional[str]) -> str:
        """Generate HTML viewer with embedded py3Dmol"""
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Protein Structure Viewer - {uniprot_id}</title>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                #viewer {{ width: 800px; height: 600px; border: 1px solid #ccc; }}
                .info {{ margin-top: 20px; }}
            </style>
        </head>
        <body>
            <h2>Protein Structure: {uniprot_id}</h2>
            <p><strong>Variant:</strong> {amino_acid_change or f"Position {residue_pos}"}</p>
            
            <div id="viewer"></div>
            
            <div class="info">
                <h3>Information</h3>
                <p>Red sphere indicates the position of the variant.</p>
                <p>Structure from AlphaFold Database (alphafold.ebi.ac.uk)</p>
            </div>
            
            <script>
                {self._generate_py3dmol_script(uniprot_id, residue_pos, amino_acid_change)}
            </script>
        </body>
        </html>
        """
        
        return html
        
    async def _map_by_gene_name(self, gene_symbol: str) -> List[str]:
        """Map gene symbol to UniProt IDs"""
        
        # Use UniProt API to search by gene name
        search_url = f"{self.uniprot_base_url}/uniprotkb/search"
        params = {
            'query': f'gene:{gene_symbol} AND organism_id:9606',  # Human only
            'format': 'json',
            'fields': 'accession,gene_names'
        }
        
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(search_url, params=params) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get('results', [])
                        return [result['primaryAccession'] for result in results]
                        
        except Exception as e:
            logger.error(f"Error mapping gene {gene_symbol}: {str(e)}")
            
        return []
        
    async def _map_by_ensembl_id(self, ensembl_id: str) -> List[str]:
        """Map Ensembl ID to UniProt IDs using ID mapping service"""
        
        mapping_data = {
            'from': 'Ensembl',
            'to': 'UniProtKB',
            'ids': ensembl_id
        }
        
        try:
            async with aiohttp.ClientSession() as session:
                # Submit mapping job
                async with session.post(f"{self.uniprot_mapping_url}/run", data=mapping_data) as response:
                    if response.status == 200:
                        job_data = await response.json()
                        job_id = job_data['jobId']
                        
                        # Poll for results
                        for _ in range(10):  # Max 10 attempts
                            await asyncio.sleep(1)
                            async with session.get(f"{self.uniprot_mapping_url}/status/{job_id}") as status_response:
                                if status_response.status == 200:
                                    status_data = await status_response.json()
                                    if status_data['jobStatus'] == 'FINISHED':
                                        # Get results
                                        async with session.get(f"{self.uniprot_mapping_url}/results/{job_id}") as results_response:
                                            if results_response.status == 200:
                                                results_data = await results_response.json()
                                                return [result['to'] for result in results_data.get('results', [])]
                                                
        except Exception as e:
            logger.error(f"Error mapping Ensembl ID {ensembl_id}: {str(e)}")
            
        return []