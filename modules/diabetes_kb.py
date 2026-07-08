"""
Diabetes gene knowledge base.

Curated from well-established sources (OMIM, PharmGKB, DIAGRAM/GWAS Catalog,
CPIC). Each entry summarises why the gene matters in diabetes and how variants
may influence risk, progression, or treatment response.

This is a curated summary for interpretation support only - it is NOT a
diagnostic tool. Always confirm on the linked databases.
"""

from typing import Dict, Any, Optional

# gene_symbol -> knowledge record
DIABETES_GENES: Dict[str, Dict[str, Any]] = {
    "TCF7L2": {
        "category": "Type 2 diabetes (strongest common risk locus)",
        "role": "Transcription factor in the Wnt pathway; affects insulin secretion.",
        "significance": "Intronic variants (e.g. rs7903146) carry the largest common-variant effect on T2D risk, impairing beta-cell insulin secretion.",
        "treatment": "Risk-allele carriers may respond less well to sulfonylureas and relatively better to metformin.",
    },
    "PPARG": {
        "category": "Type 2 diabetes / insulin resistance",
        "role": "Nuclear receptor regulating adipocyte differentiation and insulin sensitivity.",
        "significance": "Pro12Ala (rs1801282) is protective against T2D; rare loss-of-function variants cause familial partial lipodystrophy with severe insulin resistance.",
        "treatment": "Molecular target of thiazolidinediones (pioglitazone, rosiglitazone).",
    },
    "KCNJ11": {
        "category": "Neonatal diabetes / Type 2 diabetes",
        "role": "Kir6.2 subunit of the beta-cell K-ATP channel controlling insulin release.",
        "significance": "Activating mutations cause neonatal diabetes; E23K (rs5219) is a common T2D risk variant.",
        "treatment": "K-ATP neonatal diabetes often responds to sulfonylureas, allowing a switch away from insulin (CPIC/clinical guidance).",
    },
    "ABCC8": {
        "category": "Neonatal diabetes / MODY / Type 2 diabetes",
        "role": "SUR1 subunit of the beta-cell K-ATP channel.",
        "significance": "Activating mutations cause neonatal diabetes; inactivating mutations cause hyperinsulinism. Also a T2D susceptibility gene.",
        "treatment": "Sulfonylurea responsiveness is a hallmark of K-ATP-channel neonatal diabetes.",
    },
    "HNF1A": {
        "category": "MODY3 (most common MODY)",
        "role": "Transcription factor essential for beta-cell function.",
        "significance": "Loss-of-function mutations cause MODY3, a monogenic, progressive form of diabetes.",
        "treatment": "Highly sensitive to low-dose sulfonylureas - often preferred over insulin/metformin.",
    },
    "HNF4A": {
        "category": "MODY1",
        "role": "Nuclear transcription factor regulating beta-cell development.",
        "significance": "Mutations cause MODY1 with progressive beta-cell dysfunction; may present with macrosomia/neonatal hyperinsulinism.",
        "treatment": "Also typically sulfonylurea-sensitive.",
    },
    "GCK": {
        "category": "MODY2",
        "role": "Glucokinase - the beta-cell glucose sensor.",
        "significance": "Heterozygous inactivating mutations cause MODY2: mild, stable, non-progressive fasting hyperglycaemia.",
        "treatment": "Usually requires no pharmacologic treatment; distinguishing it avoids unnecessary therapy.",
    },
    "HNF1B": {
        "category": "MODY5",
        "role": "Transcription factor; renal and pancreatic development.",
        "significance": "Mutations/deletions cause MODY5 with renal cysts and pancreatic hypoplasia.",
        "treatment": "Often insulin-requiring; sulfonylureas less effective than in HNF1A/HNF4A.",
    },
    "WFS1": {
        "category": "Wolfram syndrome / Type 2 diabetes",
        "role": "Wolframin, ER membrane protein supporting beta-cell survival.",
        "significance": "Biallelic mutations cause Wolfram syndrome (diabetes + optic atrophy); common variants modestly raise T2D risk.",
        "treatment": "Wolfram diabetes is insulin-dependent.",
    },
    "SLC30A8": {
        "category": "Type 2 diabetes",
        "role": "Zinc transporter ZnT8 in beta-cell insulin granules.",
        "significance": "rs13266634 (R325W) affects T2D risk; notably, rare loss-of-function alleles are protective.",
        "treatment": "Emerging therapeutic target; no established pharmacogenomic guidance yet.",
    },
    "CDKAL1": {
        "category": "Type 2 diabetes",
        "role": "tRNA modification enzyme affecting proinsulin-to-insulin processing.",
        "significance": "Intronic variants (rs7754840) are established T2D risk loci influencing insulin secretion.",
        "treatment": "No direct treatment implication established.",
    },
    "CDKN2A": {
        "category": "Type 2 diabetes",
        "role": "Cell-cycle regulation affecting beta-cell mass.",
        "significance": "9p21 locus variants (near CDKN2A/2B) are among the replicated T2D risk signals.",
        "treatment": "No direct treatment implication established.",
    },
    "IGF2BP2": {
        "category": "Type 2 diabetes",
        "role": "RNA-binding protein influencing beta-cell development.",
        "significance": "rs4402960 is a confirmed T2D susceptibility variant affecting insulin secretion.",
        "treatment": "No direct treatment implication established.",
    },
    "FTO": {
        "category": "Obesity-mediated Type 2 diabetes",
        "role": "Fat-mass and obesity-associated gene.",
        "significance": "rs9939609 raises BMI and, secondarily, T2D risk largely through adiposity.",
        "treatment": "Lifestyle/weight management is the primary lever.",
    },
    "INS": {
        "category": "Neonatal diabetes / MODY10",
        "role": "Encodes insulin.",
        "significance": "Mutations cause permanent neonatal diabetes and rare MODY10 via misfolded proinsulin.",
        "treatment": "Insulin-dependent.",
    },
    "GLIS3": {
        "category": "Neonatal diabetes",
        "role": "Transcription factor in beta-cell development.",
        "significance": "Mutations cause neonatal diabetes with congenital hypothyroidism.",
        "treatment": "Insulin-dependent.",
    },
    "PDX1": {
        "category": "MODY4 / pancreatic agenesis",
        "role": "Master regulator of pancreas development (a.k.a. IPF1).",
        "significance": "Mutations cause MODY4; homozygous loss causes pancreatic agenesis.",
        "treatment": "Management depends on residual beta-cell function.",
    },
    "NEUROD1": {
        "category": "MODY6",
        "role": "Beta-cell transcription factor.",
        "significance": "Mutations cause MODY6 with variable presentation.",
        "treatment": "May require insulin.",
    },
    "PTPN22": {
        "category": "Type 1 diabetes (autoimmune)",
        "role": "Regulates T-cell receptor signalling.",
        "significance": "R620W (rs2476601) is a strong non-HLA T1D autoimmune risk variant.",
        "treatment": "No direct pharmacogenomic treatment guidance; relevant to autoimmune risk.",
    },
    "INSR": {
        "category": "Severe insulin resistance",
        "role": "Insulin receptor.",
        "significance": "Mutations cause syndromes of severe insulin resistance (e.g. Donohue, Rabson-Mendenhall).",
        "treatment": "Often refractory to standard therapy; specialist management.",
    },
}


def lookup(gene_symbol: Optional[str]) -> Optional[Dict[str, Any]]:
    """Return the diabetes knowledge record for a gene, or None."""
    if not gene_symbol:
        return None
    rec = DIABETES_GENES.get(gene_symbol.upper())
    if not rec:
        return None
    out = dict(rec)
    out["gene_symbol"] = gene_symbol.upper()
    # Confirmation links for the user to validate
    out["pharmgkb_url"] = f"https://www.pharmgkb.org/search?query={gene_symbol}"
    out["omim_url"] = f"https://www.omim.org/search?search={gene_symbol}"
    out["gwas_url"] = f"https://www.ebi.ac.uk/gwas/genes/{gene_symbol}"
    return out


def is_diabetes_gene(gene_symbol: Optional[str]) -> bool:
    return bool(gene_symbol) and gene_symbol.upper() in DIABETES_GENES
