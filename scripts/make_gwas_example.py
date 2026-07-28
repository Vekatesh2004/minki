#!/usr/bin/env python3
"""Generate a demonstration multi-sample VCF + phenotype table for GWAS.

This produces synthetic (not real patient) data with a small number of true
association signals, so the True Manhattan plot feature can be exercised
end-to-end. Values are simulated; they are not study results.

Usage:
    python scripts/make_gwas_example.py            # quantitative trait
    python scripts/make_gwas_example.py --binary   # case/control trait
"""

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--variants", type=int, default=120)
    parser.add_argument("--binary", action="store_true", help="Simulate a case/control trait")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", type=Path, default=Path("examples"))
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    n, m = args.samples, args.variants
    # A few causal variants scattered across chromosomes.
    causal = sorted(rng.choice(m, size=3, replace=False).tolist())
    freqs = rng.uniform(0.1, 0.5, size=m)
    geno = np.stack([rng.binomial(2, freqs[j], size=n) for j in range(m)], axis=1)

    effect = np.zeros(m)
    for index in causal:
        effect[index] = rng.uniform(0.8, 1.6) * rng.choice([-1, 1])
    linear = geno @ effect

    samples = [f"SIM{i:04d}" for i in range(n)]
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Phenotype file.
    if args.binary:
        prob = 1 / (1 + np.exp(-(linear - linear.mean())))
        trait = rng.binomial(1, prob)
        pheno_header = "sample_id,phenotype"
        pheno_rows = [f"{samples[i]},{int(trait[i])}" for i in range(n)]
    else:
        trait = linear + rng.normal(0, 1.0, size=n)
        pheno_header = "sample_id,phenotype"
        pheno_rows = [f"{samples[i]},{trait[i]:.4f}" for i in range(n)]
    pheno_path = args.outdir / "gwas_example_phenotype.csv"
    pheno_path.write_text("\n".join([pheno_header, *pheno_rows]) + "\n")

    # VCF file.
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##source=MinkiSyntheticGWASExample",
        '##MinkiNotice="Synthetic simulation for GWAS demonstration; not real study data"',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    code = {0: "0/0", 1: "0/1", 2: "1/1"}
    for j in range(m):
        chrom = str((j % 22) + 1)
        pos = 100000 + (j // 22) * 50000 + (j % 22) * 911
        gts = "\t".join(code[int(geno[i, j])] for i in range(n))
        marker = "CAUSAL" if j in causal else "null"
        lines.append(f"{chrom}\t{pos}\trs{9000000 + j}\tA\tG\t99\tPASS\tSIM={marker}\tGT\t{gts}")
    vcf_path = args.outdir / "gwas_example_genotypes.vcf"
    vcf_path.write_text("\n".join(lines) + "\n")

    print(f"Wrote {vcf_path} ({n} samples x {m} variants)")
    print(f"Wrote {pheno_path} ({'binary' if args.binary else 'quantitative'} trait)")
    print(f"Causal variant indices: {causal}")


if __name__ == "__main__":
    main()
