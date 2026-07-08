"""
VCF Parser Module - Phase 1
Handles VCF file parsing and quality control analysis
"""

import asyncio
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import statistics

try:
    import cyvcf2
except ImportError:
    cyvcf2 = None
    logging.warning("cyvcf2 not installed, falling back to basic VCF parsing")

logger = logging.getLogger(__name__)

class VCFParser:
    def __init__(self, config: Dict[str, Any]):
        """Initialize VCF parser with configuration"""
        self.config = config
        self.qc_thresholds = config.get("qc_thresholds", {})
        self.min_qual = self.qc_thresholds.get("min_qual", 30)
        self.min_depth = self.qc_thresholds.get("min_depth", 10)
        self.max_missing_rate = self.qc_thresholds.get("max_missing_rate", 0.1)
        
    async def parse_vcf(self, file_path: str, sample_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Parse VCF file and perform quality control analysis
        
        Args:
            file_path: Path to VCF file
            sample_id: Optional sample identifier
            
        Returns:
            Dictionary containing variants and QC summary
        """
        try:
            if cyvcf2:
                return await self._parse_with_cyvcf2(file_path, sample_id)
            else:
                return await self._parse_basic_vcf(file_path, sample_id)
        except Exception as e:
            logger.error(f"Error parsing VCF {file_path}: {str(e)}")
            raise
            
    async def _parse_with_cyvcf2(self, file_path: str, sample_id: Optional[str]) -> Dict[str, Any]:
        """Parse VCF using cyvcf2 library (preferred method)"""
        variants = []
        qc_metrics = {
            "total_variants": 0,
            "pass_variants": 0,
            "high_confidence_variants": 0,
            "low_confidence_variants": 0,
            "qual_scores": [],
            "depth_values": [],
            "ts_count": 0,  # transitions
            "tv_count": 0,  # transversions
            "missing_genotypes": 0,
            "het_count": 0,
            "hom_alt_count": 0,
            "hom_ref_count": 0
        }
        
        vcf = cyvcf2.VCF(file_path)
        
        for variant in vcf:
            qc_metrics["total_variants"] += 1
            
            # Extract basic variant information
            var_info = {
                "chrom": variant.CHROM,
                "pos": variant.POS,
                "id": variant.ID if variant.ID else f"{variant.CHROM}_{variant.POS}",
                "ref": variant.REF,
                "alt": variant.ALT[0] if variant.ALT else ".",
                "qual": variant.QUAL,
                "filter": variant.FILTER,
                "info": dict(variant.INFO) if hasattr(variant, 'INFO') else {},
                "format": {}
            }
            
            # Quality metrics
            if variant.QUAL is not None:
                qc_metrics["qual_scores"].append(variant.QUAL)
                
            # FILTER status
            if variant.FILTER is None or variant.FILTER == "PASS":
                qc_metrics["pass_variants"] += 1
                
            # Depth information (if available)
            # cyvcf2's variant.format('DP') returns a 2D numpy array
            # (samples x values); flatten to a single plain int.
            dp_values = variant.format('DP') if 'DP' in variant.FORMAT else None
            if dp_values is not None and len(dp_values) > 0:
                try:
                    import numpy as _np
                    flat = _np.ravel(dp_values)
                    depth = int(flat[0]) if flat.size > 0 and flat[0] is not None else 0
                except Exception:
                    depth = 0
                var_info["depth"] = depth
                qc_metrics["depth_values"].append(depth)
                
            # Genotype information
            genotypes = variant.genotypes
            if genotypes and len(genotypes) > 0:
                gt = genotypes[0][:2]  # First sample, first two alleles
                var_info["genotype"] = gt
                
                if None in gt or -1 in gt:
                    qc_metrics["missing_genotypes"] += 1
                elif gt[0] == gt[1] == 0:
                    qc_metrics["hom_ref_count"] += 1
                elif gt[0] == gt[1] and gt[0] > 0:
                    qc_metrics["hom_alt_count"] += 1
                elif gt[0] != gt[1]:
                    qc_metrics["het_count"] += 1
                    
            # Transition/Transversion classification
            if len(variant.REF) == 1 and len(var_info["alt"]) == 1:
                ts_tv = self._classify_ts_tv(variant.REF, var_info["alt"])
                var_info["ts_tv"] = ts_tv
                if ts_tv == "transition":
                    qc_metrics["ts_count"] += 1
                elif ts_tv == "transversion":
                    qc_metrics["tv_count"] += 1
                    
            # Confidence classification
            confidence = self._classify_confidence(var_info)
            var_info["confidence"] = confidence
            if confidence == "high":
                qc_metrics["high_confidence_variants"] += 1
            else:
                qc_metrics["low_confidence_variants"] += 1
                
            variants.append(var_info)
            
        # Calculate summary statistics
        qc_summary = self._calculate_qc_summary(qc_metrics)
        
        return {
            "sample_id": sample_id or Path(file_path).stem,
            "file_path": file_path,
            "variants": variants,
            "qc_summary": qc_summary,
            "total_variants": len(variants)
        }
        
    async def _parse_basic_vcf(self, file_path: str, sample_id: Optional[str]) -> Dict[str, Any]:
        """Fallback VCF parser without cyvcf2"""
        variants = []
        qc_metrics = {
            "total_variants": 0,
            "pass_variants": 0,
            "high_confidence_variants": 0,
            "low_confidence_variants": 0,
            "qual_scores": [],
            "depth_values": [],
            "ts_count": 0,
            "tv_count": 0,
            "missing_genotypes": 0,
            "het_count": 0,
            "hom_alt_count": 0,
            "hom_ref_count": 0
        }
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                    
                qc_metrics["total_variants"] += 1
                
                var_info = {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "id": fields[2] if fields[2] != '.' else f"{fields[0]}_{fields[1]}",
                    "ref": fields[3],
                    "alt": fields[4],
                    "qual": float(fields[5]) if fields[5] != '.' else None,
                    "filter": fields[6],
                    "info": self._parse_info_field(fields[7])
                }
                
                if var_info["qual"] is not None:
                    qc_metrics["qual_scores"].append(var_info["qual"])
                    
                if fields[6] == "PASS" or fields[6] == ".":
                    qc_metrics["pass_variants"] += 1
                    
                confidence = self._classify_confidence(var_info)
                var_info["confidence"] = confidence
                if confidence == "high":
                    qc_metrics["high_confidence_variants"] += 1
                else:
                    qc_metrics["low_confidence_variants"] += 1
                    
                variants.append(var_info)
                
        qc_summary = self._calculate_qc_summary(qc_metrics)
        
        return {
            "sample_id": sample_id or Path(file_path).stem,
            "file_path": file_path,
            "variants": variants,
            "qc_summary": qc_summary,
            "total_variants": len(variants)
        }
        
    def _parse_info_field(self, info_str: str) -> Dict[str, Any]:
        """Parse VCF INFO field into dictionary"""
        info_dict = {}
        if info_str == '.' or not info_str:
            return info_dict
            
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                # Try to convert to appropriate type
                try:
                    if ',' in value:
                        info_dict[key] = value.split(',')
                    elif '.' in value:
                        info_dict[key] = float(value)
                    else:
                        info_dict[key] = int(value)
                except ValueError:
                    info_dict[key] = value
            else:
                info_dict[item] = True
                
        return info_dict
        
    def _classify_ts_tv(self, ref: str, alt: str) -> str:
        """Classify SNV as transition or transversion"""
        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
        
        if (ref, alt) in transitions:
            return "transition"
        elif len(ref) == 1 and len(alt) == 1:
            return "transversion"
        else:
            return "indel"
            
    def _classify_confidence(self, variant: Dict[str, Any]) -> str:
        """Classify variant as high or low confidence based on QC metrics"""
        qual = variant.get("qual")
        depth = variant.get("depth", 0)
        filter_status = variant.get("filter", "")
        
        # High confidence criteria
        qual_pass = qual is None or qual >= self.min_qual
        depth_pass = depth >= self.min_depth
        filter_pass = filter_status in ["PASS", ".", None]
        
        if qual_pass and depth_pass and filter_pass:
            return "high"
        else:
            return "low"
            
    def _calculate_qc_summary(self, metrics: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate summary QC statistics"""
        total = metrics["total_variants"]
        
        summary = {
            "total_variants": total,
            "pass_rate": metrics["pass_variants"] / total if total > 0 else 0,
            "high_confidence_rate": metrics["high_confidence_variants"] / total if total > 0 else 0,
        }
        
        # Quality score statistics.
        # Coerce to plain Python floats first: values may be numpy scalars
        # (from cyvcf2) which statistics.mean/median cannot handle.
        qual_scores = [float(x) for x in metrics["qual_scores"] if x is not None]
        if qual_scores:
            summary["qual_stats"] = {
                "mean": statistics.mean(qual_scores),
                "median": statistics.median(qual_scores),
                "min": min(qual_scores),
                "max": max(qual_scores)
            }

        # Depth statistics
        depth_values = [float(x) for x in metrics.get("depth_values", []) if x is not None]
        if depth_values:
            summary["depth_stats"] = {
                "mean": statistics.mean(depth_values),
                "median": statistics.median(depth_values),
                "min": min(depth_values),
                "max": max(depth_values)
            }
            
        # Ts/Tv ratio
        ts_count = int(metrics.get("ts_count", 0))
        tv_count = int(metrics.get("tv_count", 0))
        if tv_count > 0:
            summary["ts_tv_ratio"] = ts_count / tv_count
            
        # Genotype statistics
        total_genotyped = metrics.get("het_count", 0) + metrics.get("hom_alt_count", 0) + metrics.get("hom_ref_count", 0)
        if total_genotyped > 0:
            summary["genotype_stats"] = {
                "heterozygous_rate": metrics.get("het_count", 0) / total_genotyped,
                "homozygous_alt_rate": metrics.get("hom_alt_count", 0) / total_genotyped,
                "homozygous_ref_rate": metrics.get("hom_ref_count", 0) / total_genotyped,
                "missing_rate": metrics.get("missing_genotypes", 0) / (total_genotyped + metrics.get("missing_genotypes", 0))
            }
            
        return summary