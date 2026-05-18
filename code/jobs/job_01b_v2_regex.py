#!/usr/bin/env python3
"""Job 01b: v2 regex categorization. Rescue RNase P → ribozyme, phages → viral_rna."""
import re, pandas as pd, sys
from pathlib import Path

ROOT = Path("/workspace/rna3d")
DATA = ROOT / "data" / "kaggle_raw"
OUT = ROOT / ".research/30_experiments/runs/01b_description_categorization_v2"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "results").mkdir(exist_ok=True)

# Read v1 assignment as base
v1 = pd.read_csv(ROOT / ".research/30_experiments/runs/01_description_categorization/results/category_assignment.csv")

# v2 patches: rescue rows currently in 'other' (and a few in v1 'complex_with_protein' that are phages)
def patch(row):
    desc_lower = str(row["description"]).lower()
    cat = row["category"]
    # Rescue RNase P / Ribonuclease P → ribozyme
    if re.search(r"\bribonuclease\s+p\b|\brnase\s+p\b", desc_lower) and cat in ("other", "complex_with_protein"):
        return "ribozyme", "v2_rnasep"
    # Rescue picornaviral / bacteriophage / coliphage → viral_rna
    if re.search(r"\bpicornavir|bacteriophage|coliphage|capsid", desc_lower) and cat == "other":
        return "viral_rna", "v2_phage"
    return cat, row.get("matched_rule", "")

v2 = v1.copy()
patched_cat, patched_rule = zip(*v2.apply(patch, axis=1))
v2["category"] = patched_cat
v2["matched_rule"] = patched_rule
v2.to_csv(OUT / "results/category_assignment.csv", index=False)

# Counts
counts = v2.groupby(["split", "category"]).size().unstack(fill_value=0).T
counts["total"] = counts.sum(axis=1)
counts = counts.sort_values("total", ascending=False)
counts.to_csv(OUT / "results/counts_by_split.csv")

with open(OUT / "results/SUMMARY.md", "w") as f:
    f.write(f"# Run 01b — v2 regex categorization (RNase P → ribozyme, phages → viral_rna)\n\n")
    f.write(f"Diff from v1: rescued from 'other'/'complex_with_protein' bucket.\n\n")
    f.write("## Per-category × split counts\n\n```\n")
    f.write(counts.to_string())
    f.write("\n```\n\n")
    diff = (v2["category"] != v1["category"]).sum()
    f.write(f"\nRows reassigned vs v1: **{diff}**\n")

print(f"OK: 01b done. {diff if (diff := (v2['category'] != v1['category']).sum()) else 0} reassignments.")
