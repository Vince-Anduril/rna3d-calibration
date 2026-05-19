"""Generate per-residue calibration figure for the paper.

Reads calibration_per_residue.csv produced by exp_f_calibration.py.
Outputs: paper/figures/per_residue_calibration.pdf  (and .png for slides).
"""
from __future__ import annotations
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
CSV = ROOT / ".research/30_experiments/runs/cpeb3_focused/calibration_per_residue.csv"
OUT = ROOT / "paper/figures"
OUT.mkdir(parents=True, exist_ok=True)

CASCADE = [9, 22, 24, 51, 60]
MUT = 30

rows = list(csv.DictReader(open(CSV)))
res = np.array([int(r["residue"]) for r in rows])
drf_err = np.array([float(r["drfold2_err_vs_crystal"]) for r in rows])
af3_err = np.array([float(r["af3_err_vs_crystal"]) for r in rows])
drf_std = np.array([float(r["drfold2_ens_std"]) for r in rows])
af3_plddt = np.array([float(r["af3_plddt"]) for r in rows])

fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True,
                         gridspec_kw={"height_ratios": [1.4, 1.0]})

# --- top: per-residue ERROR vs crystal ---
ax = axes[0]
ax.fill_between(res, 0, af3_err, color="#ef4444", alpha=0.35, label="AF3 error vs 7QR3:C")
ax.fill_between(res, 0, drf_err, color="#0ea5e9", alpha=0.55, label="DRfold2 error vs 7QR3:C")
ax.plot(res, af3_err, color="#b91c1c", linewidth=1.0)
ax.plot(res, drf_err, color="#0369a1", linewidth=1.4)
for r in CASCADE:
    ax.axvline(r, color="#94a3b8", linestyle=":", linewidth=0.8, alpha=0.7)
ax.axvline(MUT, color="#f59e0b", linestyle="-", linewidth=1.4, alpha=0.7,
           label=f"mutation (pos {MUT})")
for r, lbl in [(9, "P1 anchor"), (60, "P1 end")]:
    ax.annotate(lbl, xy=(r, max(af3_err[r-1], drf_err[r-1])), xytext=(r, 35),
                ha="center", fontsize=8, color="#475569",
                arrowprops=dict(arrowstyle="-", color="#94a3b8", lw=0.6))
ax.set_ylabel("Per-residue error vs crystal (Å)")
ax.set_xlim(1, 69)
ax.set_ylim(0, max(af3_err.max(), 35) * 1.1)
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.set_title("Per-residue calibration on R1108 (chimp CPEB3) vs 7QR3:C crystal",
             loc="left", fontsize=11, fontweight="bold")

# --- bottom: per-residue UNCERTAINTY signals ---
ax2 = axes[1]
ax2.bar(res - 0.2, drf_std, width=0.4, color="#0ea5e9", alpha=0.85,
        label="DRfold2 ensemble std (high = uncertain)")
ax2.set_ylabel("DRfold2 ensemble std (Å)", color="#0369a1")
ax2.tick_params(axis="y", labelcolor="#0369a1")
ax2.set_ylim(0, drf_std.max() * 1.15)

ax3 = ax2.twinx()
ax3.bar(res + 0.2, af3_plddt, width=0.4, color="#ef4444", alpha=0.55,
        label="AF3 pLDDT (low = uncertain)")
ax3.set_ylabel("AF3 pLDDT (0–100, low = uncertain)", color="#b91c1c")
ax3.tick_params(axis="y", labelcolor="#b91c1c")
ax3.set_ylim(0, max(50, af3_plddt.max() * 1.15))
ax3.axhline(50, color="#b91c1c", linestyle=":", linewidth=0.8, alpha=0.6)
ax3.annotate("pLDDT 50 = AF3 'low confidence' threshold", xy=(50, 51), fontsize=7,
             color="#b91c1c", alpha=0.7)
for r in CASCADE:
    ax2.axvline(r, color="#94a3b8", linestyle=":", linewidth=0.8, alpha=0.7)
ax2.axvline(MUT, color="#f59e0b", linestyle="-", linewidth=1.4, alpha=0.7)
ax2.set_xlabel("Residue position (1–69)")
# combined legend
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax3.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper right",
           frameon=False, fontsize=8)

# annotate cascade residues at top
for r in CASCADE + [MUT]:
    axes[0].annotate(str(r), xy=(r, -2.5), ha="center", va="top", fontsize=7,
                     color="#475569", xycoords=("data", "data"),
                     annotation_clip=False)

plt.tight_layout()
out_pdf = OUT / "per_residue_calibration.pdf"
out_png = OUT / "per_residue_calibration.png"
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, bbox_inches="tight", dpi=200)
print(f"Written: {out_pdf}\nWritten: {out_png}")
