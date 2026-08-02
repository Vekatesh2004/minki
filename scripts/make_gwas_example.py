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
    parser.add_argument("--target-mib", type=float, help="Approximate VCF size in MiB using non-analytic INFO padding")
    parser.add_argument("--prefix", default="gwas_example", help="Output filename prefix")
    parser.add_argument("--binary", action="store_true", help="Simulate a case/control trait")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", type=Path, default=Path("examples"))
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    n, m = args.samples, args.variants
    if n < 10:
        parser.error("--samples must be at least 10 for association testing")
    if m < 3:
        parser.error("--variants must be at least 3")
    if args.target_mib is not None and args.target_mib <= 0:
        parser.error("--target-mib must be positive")

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
    pheno_path = args.outdir / f"{args.prefix}_phenotype.csv"
    pheno_path.write_text("\n".join([pheno_header, *pheno_rows]) + "\n")

    # VCF file.
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##source=MinkiSyntheticGWASExample",
        '##MinkiNotice="Synthetic simulation for GWAS demonstration; not real study data"',
        '##INFO=<ID=SIM,Number=1,Type=String,Description="Synthetic causal marker label">',
        '##INFO=<ID=PAD,Number=1,Type=String,Description="Non-analytic padding used only to create a requested test file size">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    code = {0: "0/0", 1: "0/1", 2: "1/1"}
    row_data = []
    for j in range(m):
        chrom = str((j % 22) + 1)
        pos = 100000 + (j // 22) * 50000 + (j % 22) * 911
        gts = "\t".join(code[int(geno[i, j])] for i in range(n))
        marker = "CAUSAL" if j in causal else "null"
        row_data.append((chrom, pos, j, marker, gts))

    # When a target size is requested, add standards-compliant INFO padding.
    # This keeps variant/sample counts computationally manageable while
    # producing a realistic upload-size test fixture. PAD is ignored by GWAS.
    padding_per_row = 0
    if args.target_mib:
        target_bytes = int(args.target_mib * 1024 * 1024)
        base_header = "\n".join(lines) + "\n"
        base_rows = [
            f"{chrom}\t{pos}\trs{9000000 + j}\tA\tG\t99\tPASS\tSIM={marker}\tGT\t{gts}\n"
            for chrom, pos, j, marker, gts in row_data
        ]
        base_size = len(base_header.encode("utf-8")) + sum(len(row.encode("utf-8")) for row in base_rows)
        if target_bytes > base_size:
            padding_per_row = max(0, (target_bytes - base_size) // m - len(";PAD="))

    rows = []
    for chrom, pos, j, marker, gts in row_data:
        info = f"SIM={marker}"
        if padding_per_row:
            info += ";PAD=" + ("X" * padding_per_row)
        rows.append(f"{chrom}\t{pos}\trs{9000000 + j}\tA\tG\t99\tPASS\t{info}\tGT\t{gts}")

    vcf_path = args.outdir / f"{args.prefix}_genotypes.vcf"
    vcf_path.write_text("\n".join(lines + rows) + "\n")

    size_mib = vcf_path.stat().st_size / (1024 * 1024)
    print(f"Wrote {vcf_path} ({n} samples x {m} variants; {size_mib:.2f} MiB)")
    print(f"Wrote {pheno_path} ({'binary' if args.binary else 'quantitative'} trait)")
    print(f"Causal variant indices: {causal}")


if __name__ == "__main__":
    main()
