"""Cross-family figure for the paper.

Bar chart: DRfold2 vs AF3 (default + best of 5) per target, with crystal
noise floor reference. Color-coded by training-cutoff status.
"""
from __future__ import annotations
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
CSV = ROOT / ".research/30_experiments/runs/multi_rna/multi_rna_results.csv"
OUT = ROOT / "paper/figures"
OUT.mkdir(parents=True, exist_ok=True)

rows = list(csv.DictReader(open(CSV)))
rows = [r for r in rows if r.get("drfold2_rmsd_vs_crystal") not in (None, "", "None")]
# parse
ids = [r["id"] for r in rows]
names = [r["name"] for r in rows]
drf = [float(r["drfold2_rmsd_vs_crystal"]) for r in rows]
def _val(r, k):
    v = r.get(k)
    if v in (None, "", "None"): return np.nan
    try: return float(v)
    except (TypeError, ValueError): return np.nan
af3_def = [_val(r, "af3_rmsd_default") for r in rows]
af3_best = [_val(r, "af3_rmsd_best") for r in rows]
def _seen(s: str) -> bool:
    s = s.lower()
    return ("likely yes" in s) or s.strip() == "yes" or s.strip().startswith("yes")
drf_seen = [_seen(r["drfold_seen"]) for r in rows]
af3_seen = [_seen(r["af3_seen"]) for r in rows]

fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(ids))
w = 0.27

bars1 = ax.bar(x - w, drf, w, color="#0ea5e9", edgecolor="#0369a1", label="DRfold2 (single-seq LM)")
bars2 = ax.bar(x,     af3_def, w, color="#fca5a5", edgecolor="#b91c1c", label="AF3 default seed")
bars3 = ax.bar(x + w, af3_best, w, color="#fde68a", edgecolor="#a16207", label="AF3 best of 5 seeds")

# crystal noise floor for 7QR3 only (we have C↔D dimer)
ax.axhline(2.60, color="#10b981", linestyle="--", linewidth=1, alpha=0.7)
ax.text(len(ids)-0.5, 2.7, "crystal noise floor (2.60 Å, 7QR3 C↔D)", color="#10b981",
        fontsize=8, ha="right", alpha=0.85)

# annotations
for i, v in enumerate(drf):
    ax.text(x[i]-w, v+0.2, f"{v:.1f}", ha="center", fontsize=8, color="#0369a1", fontweight="bold")
for i, v in enumerate(af3_def):
    if not np.isnan(v):
        ax.text(x[i], v+0.2, f"{v:.1f}", ha="center", fontsize=8, color="#b91c1c")
for i, v in enumerate(af3_best):
    if not np.isnan(v):
        ax.text(x[i]+w, v+0.2, f"{v:.1f}", ha="center", fontsize=8, color="#a16207")

# x labels with cutoff status
labels = []
for r, ds, as_ in zip(rows, drf_seen, af3_seen):
    tag = ""
    if not as_ and not ds:
        tag = "\n(both blind)"
    elif ds and not as_:
        tag = "\n(DRfold2 may have seen)"
    elif ds and as_:
        tag = "\n(both may have seen)"
    labels.append(f"{r['id']}\n{r['family']}{tag}")
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Kabsch C1' RMSD vs crystal (Å)", fontsize=11)
ax.set_title(f"Cross-family ground-truth comparison (N = {len(ids)} RNAs, single-chain AF3)",
             fontsize=11, fontweight="bold", loc="left")
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim(0, max(max(af3_def + drf + af3_best, default=0), 12) * 1.1)

# aggregate stats footer
import numpy as np
def safemean(xs): return np.nanmean(xs) if len(xs) else np.nan
footer = (f"Means: DRfold2={safemean(drf):.2f} Å, "
          f"AF3 default={safemean(af3_def):.2f} Å, AF3 best={safemean(af3_best):.2f} Å. "
          f"DRfold2 wins {sum(1 for d,a in zip(drf, af3_def) if not np.isnan(a) and d < a)}/{len(drf)} comparisons.")
fig.text(0.04, 0.005, footer, fontsize=8, color="#475569", style="italic")
plt.tight_layout(rect=(0, 0.02, 1, 1))
plt.savefig(OUT / "multi_rna_rmsd.pdf", bbox_inches="tight")
plt.savefig(OUT / "multi_rna_rmsd.png", bbox_inches="tight", dpi=200)
print(f"Written: {OUT / 'multi_rna_rmsd.pdf'}")
print(f"Written: {OUT / 'multi_rna_rmsd.png'}")
