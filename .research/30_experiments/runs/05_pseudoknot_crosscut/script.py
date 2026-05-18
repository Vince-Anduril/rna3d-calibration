"""Job 5 — Pseudoknot crosscut by category.

Joins:
  - run 01b category_assignment.csv  (sequence -> category)
  - run 04 ground_truth.csv          (pdb chain -> has_pseudoknot)

Matching strategy: target_id in Stanford CSVs is typically `<pdb_id>_<chain>` (e.g.
`1SCL_A`). We extract pdb_id (lowercased) and chain to join with run 04. If the
join fails, we still emit per-category counts of how many have *any* pseudoknot-positive
chain by pdb_id.

Outputs:
  results/SUMMARY.md  (category x pseudoknot table)
  results/crosscut.csv (per-target join)
"""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True, parents=True)
LOG = HERE / "log.txt"

CAT_CSV = Path("/workspace/rna3d/.research/30_experiments/runs/01b_description_categorization_v2/results/category_assignment.csv")
GT_CSV = Path("/workspace/rna3d/.research/30_experiments/runs/04_ground_truth_structural_metrics/results/ground_truth.csv")


def log(m):
    line = f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {m}"
    print(line, flush=True)
    LOG.open("a").write(line + "\n")


def main():
    log("=== Job 5 pseudoknot crosscut: START ===")
    if not CAT_CSV.exists():
        log(f"[abort] missing {CAT_CSV} (run 01b must succeed first)")
        (RES / "SUMMARY.md").write_text(
            "# Job 5 — SKIPPED (Job 1 / 01b output missing)\n"
        )
        return
    if not GT_CSV.exists():
        log(f"[abort] missing {GT_CSV} (run 04 must succeed first)")
        (RES / "SUMMARY.md").write_text(
            "# Job 5 — SKIPPED (Job 4 ground-truth output missing)\n"
        )
        return

    cat = pd.read_csv(CAT_CSV)
    gt = pd.read_csv(GT_CSV)
    log(f"loaded cat={len(cat)} gt={len(gt)}")

    # Parse target_id -> pdb_id + chain
    parts = cat["target_id"].astype(str).str.split("_", n=1, expand=True)
    cat["pdb_id_lc"] = parts[0].str.lower()
    cat["chain_guess"] = parts[1] if 1 in parts.columns else ""

    gt["pdb_id_lc"] = gt["pdb_id"].astype(str).str.lower()

    # Join by pdb_id (any chain pseudoknot-positive => target has_pseudoknot)
    gt_agg = gt.groupby("pdb_id_lc").agg(
        has_pseudoknot_any=("has_pseudoknot", lambda s: bool(s.fillna(False).any())),
        max_crossed=("n_crossed_pairs", "max"),
        chains_analyzed=("chain_id", "count"),
    ).reset_index()

    merged = cat.merge(gt_agg, on="pdb_id_lc", how="left")
    merged["has_pseudoknot_any"] = merged["has_pseudoknot_any"].fillna(False)
    merged.to_csv(RES / "crosscut.csv", index=False)

    # Pivot: category x pseudoknot
    tbl = (merged.groupby(["category", "has_pseudoknot_any"])
                 .size().unstack(fill_value=0))
    if True not in tbl.columns:
        tbl[True] = 0
    if False not in tbl.columns:
        tbl[False] = 0
    tbl = tbl[[False, True]]
    tbl.columns = ["pk_absent", "pk_present"]
    tbl["total"] = tbl["pk_absent"] + tbl["pk_present"]
    tbl["pk_rate"] = (tbl["pk_present"] / tbl["total"].clip(lower=1)).round(3)
    tbl = tbl.sort_values("pk_present", ascending=False)

    lines = [
        "# Job 5 — Pseudoknot crosscut by category — SUMMARY",
        "",
        f"Joined {len(merged)} Stanford sequences to {len(gt_agg)} unique PDB ids in ground-truth.",
        f"Pseudoknot-positive targets (any chain): {int(merged['has_pseudoknot_any'].sum())} / {len(merged)}",
        "",
        "## Category x pseudoknot",
        "",
        "| category | pk_absent | pk_present | total | pk_rate |",
        "|---|---:|---:|---:|---:|",
    ]
    for cat_name, row in tbl.iterrows():
        lines.append(f"| {cat_name} | {int(row['pk_absent'])} | {int(row['pk_present'])} | {int(row['total'])} | {row['pk_rate']:.3f} |")
    lines += [
        "",
        "Pseudoknot definition: PDB ground-truth has at least one pair (i,j) crossed by another pair (k,l) with i<k<j<l, on any chain.",
        "Note: join is pdb_id-level (case-insensitive). Multi-chain targets inherit the OR across chains.",
    ]
    (RES / "SUMMARY.md").write_text("\n".join(lines))
    log("=== Job 5: DONE ===")


if __name__ == "__main__":
    main()
