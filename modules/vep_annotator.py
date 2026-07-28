"""
VEP Annotator Module - Phase 2
Handles variant annotation using Ensembl VEP REST API
"""

import asyncio
import logging
from typing import Any, Dict, List, Optional
import aiohttp
import json
import time

logger = logging.getLogger(__name__)

class VEPAnnotator:
    def __init__(self, config: Dict[str, Any]):
        """Initialize VEP annotator with configuration"""
        self.config = config
        self.vep_config = config.get("ensembl_vep", {})
        self.base_url = self.vep_config.get("base_url", "https://rest.ensembl.org/vep/human/region")
        self.rate_limit_delay = self.vep_config.get("rate_limit_delay", 0.1)
        
        # Coding consequence terms we're interested in
        self.coding_consequences = {
            "missense_variant",
            "stop_gained", 
            "stop_lost",
            "frameshift_variant",
            "inframe_insertion",
            "inframe_deletion", 
            "synonymous_variant",
            "start_lost",
            "stop_retained_variant"
        }
        
    async def annotate_variants(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Annotate variants using Ensembl VEP
        
        Args:
            variants: List of variant dictionaries from VCF parser
            
        Returns:
            List of annotated variants with gene and consequence information
        """
        try:
            # Process variants in batches to respect rate limits
            batch_size = 200  # VEP REST API batch limit
            annotated_variants = []
            
            for i in range(0, len(variants), batch_size):
                batch = variants[i:i + batch_size]
                batch_results = await self._annotate_batch(batch)
                annotated_variants.extend(batch_results)
                
                # Rate limiting between batches
                if i + batch_size < len(variants):
                    await asyncio.sleep(self.rate_limit_delay)
                    
            return annotated_variants
            
        except Exception as e:
            logger.error(f"Error in variant annotation: {str(e)}")
            raise
            
    async def _annotate_batch(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Annotate a batch of variants using VEP REST API"""
        
        # Prepare VEP input format
        vep_input = []
        for variant in variants:
            # VEP format: "chromosome:start-end:strand/alleles"
            # For SNPs: "21:26960070-26960070:1/C"
            # For indels: "22:17662833-17662834:1/-" (deletion)
            
            chrom = str(variant["chrom"]).replace("chr", "")  # Remove chr prefix if present
            pos = variant["pos"]
            # Uppercase alleles (some VCFs, e.g. VCFv3.x, use lowercase bases)
            ref = str(variant["ref"]).upper()
            alt = str(variant["alt"]).upper()
            
            # Handle different variant types
            if len(ref) == 1 and len(alt) == 1:
                # SNP
                vep_variant = f"{chrom}:{pos}-{pos}:1/{alt}"
            elif len(ref) > len(alt):
                # Deletion
                start_pos = pos + len(alt)
                end_pos = pos + len(ref) - 1
                vep_variant = f"{chrom}:{start_pos}-{end_pos}:1/-"
            elif len(alt) > len(ref):
                # Insertion
                insert_seq = alt[len(ref):]
                vep_variant = f"{chrom}:{pos}-{pos+1}:1/{insert_seq}"
            else:
                # Complex variant - use original notation
                vep_variant = f"{chrom}:{pos}-{pos+len(ref)-1}:1/{alt}"
                
            vep_input.append(vep_variant)
            
        # Make VEP API request
        vep_results = await self._query_vep_api(vep_input)
        
        # Process and merge results
        annotated_variants = []
        for i, variant in enumerate(variants):
            annotated_variant = variant.copy()
            
            # Find corresponding VEP result
            vep_result = vep_results[i] if i < len(vep_results) else {}
            
            # Extract annotation information
            annotation = self._extract_annotations(vep_result)
            annotated_variant.update(annotation)
            
            annotated_variants.append(annotated_variant)
            
        return annotated_variants
        
    async def _query_vep_api(self, vep_input: List[str]) -> List[Dict[str, Any]]:
        """Query VEP with bounded retries for transient service failures."""
        headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
        data = {
            'variants': vep_input,
            'canonical': 1,
            'hgvs': 1,
            'protein': 1,
            'xref_refseq': 1,
            'uniprot': 1,
            'domains': 1,
        }
        retryable = {429, 500, 502, 503, 504}

        for attempt in range(3):
            try:
                timeout = aiohttp.ClientTimeout(total=60)
                async with aiohttp.ClientSession(timeout=timeout) as session:
                    async with session.post(self.base_url, headers=headers, json=data) as response:
                        if response.status == 200:
                            result = await response.json()
                            if isinstance(result, list):
                                return result
                            logger.error("VEP returned an unexpected response shape")
                            return []
                        response_text = await response.text()
                        if response.status not in retryable:
                            logger.error("VEP API error %s: %s", response.status, response_text[:500])
                            return []
                        logger.warning(
                            "Temporary VEP API error %s (attempt %s/3): %s",
                            response.status, attempt + 1, response_text[:200],
                        )
            except (aiohttp.ClientError, asyncio.TimeoutError) as exc:
                logger.warning("VEP request failed (attempt %s/3): %s", attempt + 1, exc)
            except Exception as exc:
                logger.error("Unexpected VEP request error: %s", exc)
                return []

            if attempt < 2:
                await asyncio.sleep(2 ** attempt)

        logger.error("VEP remained unavailable after 3 attempts")
        return []
            
    def _extract_annotations(self, vep_result: Dict[str, Any]) -> Dict[str, Any]:
        """Extract relevant annotations from VEP result"""
        
        annotation = {
            "vep_annotation": {
                "gene_symbol": None,
                "gene_id": None,
                "transcript_id": None,
                "consequence_terms": [],
                "hgvsc": None,
                "hgvsp": None,
                "canonical": False,
                "protein_id": None,
                "uniprot_ids": [],
                "domains": [],
                "coding": False,
                "amino_acid_position": None,
                "amino_acid_change": None
            }
        }
        
        if not vep_result or 'transcript_consequences' not in vep_result:
            return annotation
            
        # Process transcript consequences
        transcript_consequences = vep_result.get('transcript_consequences', [])
        
        # Prioritize canonical transcript, then coding consequences
        best_consequence = None
        for consequence in transcript_consequences:
            
            # Check if this consequence involves coding changes
            cons_terms = consequence.get('consequence_terms', [])
            has_coding_consequence = any(term in self.coding_consequences for term in cons_terms)
            
            if not best_consequence:
                best_consequence = consequence
            elif consequence.get('canonical') and not best_consequence.get('canonical'):
                best_consequence = consequence
            elif has_coding_consequence and not any(term in self.coding_consequences 
                                                 for term in best_consequence.get('consequence_terms', [])):
                best_consequence = consequence
                
        if best_consequence:
            ann = annotation["vep_annotation"]
            
            # Basic gene information
            ann["gene_symbol"] = best_consequence.get('gene_symbol')
            ann["gene_id"] = best_consequence.get('gene_id')
            ann["transcript_id"] = best_consequence.get('transcript_id')
            ann["consequence_terms"] = best_consequence.get('consequence_terms', [])
            ann["canonical"] = best_consequence.get('canonical', False)
            
            # HGVS notation
            ann["hgvsc"] = best_consequence.get('hgvsc')
            ann["hgvsp"] = best_consequence.get('hgvsp')
            
            # Protein information
            ann["protein_id"] = best_consequence.get('protein_id')
            
            # Extract UniProt IDs if available
            if 'xref_ids' in best_consequence:
                uniprot_ids = [xref for xref in best_consequence['xref_ids'] 
                             if xref.startswith('UniProt')]
                ann["uniprot_ids"] = uniprot_ids
                
            # Protein domains
            ann["domains"] = best_consequence.get('domains', [])
            
            # Check if coding consequence
            ann["coding"] = any(term in self.coding_consequences 
                              for term in ann["consequence_terms"])
            
            # Extract amino acid position and change from HGVSp
            if ann["hgvsp"]:
                aa_info = self._parse_hgvsp(ann["hgvsp"])
                ann["amino_acid_position"] = aa_info.get("position")
                ann["amino_acid_change"] = aa_info.get("change")
                
        # Add additional VEP information
        annotation["vep_annotation"]["vep_version"] = vep_result.get('assembly_name')
        annotation["vep_annotation"]["most_severe_consequence"] = vep_result.get('most_severe_consequence')
        
        return annotation
        
    def _parse_hgvsp(self, hgvsp: str) -> Dict[str, Any]:
        """Parse HGVSp notation to extract amino acid position and change"""
        
        aa_info = {
            "position": None,
            "change": None,
            "ref_aa": None,
            "alt_aa": None
        }
        
        if not hgvsp:
            return aa_info

        # VEP returns HGVSp with a transcript prefix, e.g.
        # "ENSP00000414462.2:p.Lys120Gln". Strip everything up to and
        # including "p." so we are left with just the change string.
        if ':p.' in hgvsp:
            change_str = hgvsp.split(':p.', 1)[1]
        elif hgvsp.startswith('p.'):
            change_str = hgvsp[2:]
        else:
            return aa_info
        
        # Handle different types of changes
        if '=' in change_str:
            # Synonymous - p.Leu54=
            aa_info["change"] = "synonymous"
            # Extract position
            import re
            pos_match = re.search(r'(\d+)', change_str)
            if pos_match:
                aa_info["position"] = int(pos_match.group(1))
                
        elif '*' in change_str:
            # Stop codon - p.Arg1456*
            aa_info["change"] = "nonsense"
            import re
            pos_match = re.search(r'(\d+)', change_str)
            if pos_match:
                aa_info["position"] = int(pos_match.group(1))
                
        else:
            # Missense - p.Arg1456Cys
            import re
            # Match pattern like Arg1456Cys
            match = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3})', change_str)
            if match:
                aa_info["ref_aa"] = match.group(1)
                aa_info["position"] = int(match.group(2))
                aa_info["alt_aa"] = match.group(3)
                aa_info["change"] = f"{aa_info['ref_aa']}{aa_info['position']}{aa_info['alt_aa']}"
                
        return aa_info