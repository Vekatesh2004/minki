# All-Gene Disease Badge Feature

## Overview
Extended the disease badge feature to show disease associations for ALL genes in the Coding Mutations table, not just diabetes genes.

## Implementation Date
October 12, 2026 (19:33)

## What Changed

### Before
- Only genes in the diabetes gene list (KCNJ11, PPARG, SLC30A8, etc.) showed disease badges
- Badge only said "Diabetes" in teal color

### After
- ALL genes with ClinVar disease associations now show disease badges
- Multiple badges per gene (up to 3 diseases shown)
- Color-coded badges by disease category:
  - **Teal/Cyan** (#0891b2): Diabetes, Metabolic disorders
  - **Red** (#dc2626): Cancer, Tumors, Carcinoma
  - **Purple** (#9333ea): Alzheimer's, Parkinson's, Neurological
  - **Orange** (#f97316): Obesity, Lipid disorders
  - **Cyan** (#06b6d4): Cystic Fibrosis, Other conditions

### Example Display
- **KCNJ11** [T2D] [Neonatal Diabetes] [Hyperinsulinism]
- **PPARG** [T2D] [Obesity] [Lipid Disorder]
- **BRCA1** [Breast Cancer] [Ovarian Cancer]
- **CFTR** [Cystic Fibrosis]
- **APOE** [Alzheimer's]

## Technical Implementation

### Data Source
- Disease information is extracted from ClinVar findings
- ClinVar `conditions` field contains disease associations for each gene
- Multiple conditions are split by semicolons/commas

### Code Location
`simple_backend.py` - Lines ~1504-1580 (renderResults function)

### Algorithm
1. Build `geneDiseasesMap` from ClinVar findings
2. For each mutation in the table:
   - Check if gene has ClinVar disease associations
   - Show up to 3 disease badges
   - Apply color coding based on disease type
   - Shorten long disease names

### Disease Name Abbreviations
- "Type 2 diabetes mellitus" → "T2D"
- "Type 1 diabetes mellitus" → "T1D"  
- "Neonatal diabetes" → "Neonatal Diabetes"
- "Breast cancer" → "Breast Cancer"
- Long names (>25 chars) → Truncated with "..."

## Benefits
1. **Comprehensive**: Shows disease info for ALL genes, not just diabetes
2. **Visual**: Color-coded badges make it easy to spot disease categories
3. **Informative**: Multiple diseases per gene provide richer context
4. **Non-intrusive**: Limited to 3 badges per gene to avoid clutter
5. **Accessible**: Works with existing ClinVar data, no new API calls needed

## Files Modified
- `simple_backend.py` - Updated mutation table rendering

## Testing
Test with:
1. **Example file** - Should show diabetes badges for KCNJ11, PPARG, SLC30A8
2. **Sample with BRCA1** - Should show cancer badges
3. **Sample with CFTR** - Should show cystic fibrosis badge
4. **Sample with APOE** - Should show Alzheimer's badge

## Future Enhancements
- Add tooltip with full disease name on hover
- Include more disease-specific colors
- Add evidence level indicators
- Link badges to disease information resources

## Notes
- Python warnings about `\s` escape sequences in JavaScript regex are harmless
- Badges only appear if gene has ClinVar disease associations
- Works for any VCF file, not just diabetes-focused samples
- Preserves all existing functionality
