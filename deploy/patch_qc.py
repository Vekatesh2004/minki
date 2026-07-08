#!/usr/bin/env python3
"""
Idempotent patch: make QC stats numpy-safe in modules/vcf_parser.py.
Run on the server if you can't git-pull the latest code:
    python deploy/patch_qc.py
"""
import re
from pathlib import Path

p = Path("modules/vcf_parser.py")
src = p.read_text()

old_qual = '''        # Quality score statistics
        if metrics["qual_scores"]:
            summary["qual_stats"] = {
                "mean": statistics.mean(metrics["qual_scores"]),
                "median": statistics.median(metrics["qual_scores"]),
                "min": min(metrics["qual_scores"]),
                "max": max(metrics["qual_scores"])
            }
            
        # Depth statistics
        if metrics.get("depth_values"):
            summary["depth_stats"] = {
                "mean": statistics.mean(metrics["depth_values"]),
                "median": statistics.median(metrics["depth_values"]),
                "min": min(metrics["depth_values"]),
                "max": max(metrics["depth_values"])
            }'''

new_qual = '''        # Quality score statistics (coerce numpy scalars to float)
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
            }'''

old_ts = '''        ts_count = metrics.get("ts_count", 0)
        tv_count = metrics.get("tv_count", 0)'''
new_ts = '''        ts_count = int(metrics.get("ts_count", 0))
        tv_count = int(metrics.get("tv_count", 0))'''

changed = False
if old_qual in src:
    src = src.replace(old_qual, new_qual)
    changed = True
if old_ts in src:
    src = src.replace(old_ts, new_ts)
    changed = True

if changed:
    p.write_text(src)
    print("Patched modules/vcf_parser.py")
else:
    print("Already patched (or patterns not found) - no changes made")
