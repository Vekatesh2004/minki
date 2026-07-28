#!/usr/bin/env python3
"""Generate a multi-sample VCF for the genotype-only (Hardy-Weinberg) Manhattan plot.

No phenotype file is needed for this mode. The output is a synthetic
multi-sample VCF where most variants sit in Hardy-Weinberg equilibrium and a
handful are deliberately distorted (excess/deficit of heterozygotes) so the
HWE Manhattan plot shows a few genuine peaks. Values are simulated, not real
patient data.

Usage:
    python scripts/make_hwe_example.py
    python scripts/make_hwe_example.py --samples 400 --variants 600
"""

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=300)
    parser.add_argument("--variants", type=int, default=500)
    parser.add_argument("--seed", type=int, default=24)
    parser.add_argument(
        "--out", type=Path, default=Path("examples/hwe_example_multisample.vcf"),
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    n, m = args.samples, args.variants
    samples = [f"IND{i:04d}" for i in range(n)]

    # Pick a few variants to violate HWE so the plot has visible peaks.
    distorted = set(rng.choice(m, size=max(3, m // 60), replace=False).tolist())

    header = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##source=MinkiSyntheticHWEExample",
        '##INFO=<ID=SIM,Number=1,Type=String,Description="Synthetic HWE simulation label (equilibrium|distorted)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    code = {0: "0/0", 1: "0/1", 2: "1/1"}
    bases = "ACGT"
    rows = []
    for j in range(m):
        chrom = str((j % 22) + 1)
        pos = 100000 + (j // 22) * 75000 + (j % 22) * 811
        freq = float(rng.uniform(0.1, 0.5))

        if j in distorted:
            # Force an excess of heterozygotes: most samples become 0/1.
            dosage = rng.choice([0, 1, 2], size=n, p=[0.05, 0.90, 0.05])
            label = "distorted"
        else:
            # Draw genotypes directly under HWE for this allele frequency.
            dosage = rng.binomial(2, freq, size=n)
            label = "equilibrium"

        ref = bases[rng.integers(0, 4)]
        alt = bases[(bases.index(ref) + 1 + rng.integers(0, 3)) % 4]
        gts = "\t".join(code[int(d)] for d in dosage)
        rows.append(
            f"{chrom}\t{pos}\trs{8000000 + j}\t{ref}\t{alt}\t99\tPASS\tSIM={label}\tGT\t{gts}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(header + rows) + "\n")
    print(f"Wrote {args.out}")
    print(f"  samples : {n}")
    print(f"  variants: {m} across chromosomes 1-22")
    print(f"  distorted (HWE-violating) variants: {sorted(distorted)}")
    print("Upload this file in the GWAS panel and leave the phenotype field empty.")


if __name__ == "__main__":
    main()
