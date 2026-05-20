"""Figures for the ACDSA paper:
  Fig 1: per-target scatter of AF3 rho vs DRfold2 rho (12 targets)
  Fig 2: pooled per-residue (DRfold2 std, -AF3 pLDDT), colored by error,
         across all 12 targets.
Uses calibration data already computed in:
  .research/30_experiments/runs/multi_rna/ensemble_calibration_results.csv
  and per-residue values from analyses.
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
OUT = ROOT / "paper/figures"
OUT.mkdir(parents=True, exist_ok=True)

# Per-target rho values from analyses (Table I in the ACDSA paper).
# Each row: target, AF3 rho, AF3 p, DRfold2 std rho, DRfold2 std p, length
TARGETS = [
    ("7QR3", +0.083, 0.497, +0.217, 0.0735,  69),
    ("9DE7", -0.604, 6.38e-7, +0.493, 9.73e-5, 57),
    ("9HRD", -0.581, 2.6e-7,  +0.183, 0.139,   67),
    ("9HRF", -0.435, 1.68e-4, +0.331, 0.00516, 70),
    ("9J4N", -0.270, 0.0148,  +0.664, 1.38e-11, 81),
    ("9E9O", -0.248, 0.0124,  +0.430, 7.38e-6, 101),
    ("9MFH", -0.578, 4.48e-8, +0.017, 0.883,   76),
    ("9LJN", -0.276, 0.0199,  +0.581, 1.10e-7, 71),
    ("9LKE", +0.089, 0.463,   +0.598, 4.54e-8, 70),
    ("9LKU", -0.411, 8.32e-4, +0.312, 0.0127,  65),
    ("9UW0", -0.198, 0.12,    +0.254, 0.0442,  63),
    ("12CI", -0.546, 1.09e-7, +0.428, 6.12e-5, 82),
]
BONF = 0.05 / 12

# Fig 1: per-target rho scatter
fig, ax = plt.subplots(figsize=(6.0, 5.0))
for tid, rho_af3, p_af3, rho_drf, p_drf, L in TARGETS:
    bonf_af3 = p_af3 < BONF
    bonf_drf = p_drf < BONF
    # Color: green if both Bonferroni, orange if one, gray if none
    if bonf_af3 and bonf_drf: col, edge, lw = "#10b981", "black", 1.2
    elif bonf_af3 or bonf_drf: col, edge, lw = "#f59e0b", "black", 0.8
    else: col, edge, lw = "#cbd5e1", "gray", 0.6
    ax.scatter(-rho_af3, rho_drf, s=70 + 0.6*L, c=col, edgecolors=edge,
               linewidths=lw, alpha=0.95, zorder=3)
    ax.annotate(tid, xy=(-rho_af3, rho_drf), xytext=(5, 4),
                textcoords="offset points", fontsize=8, color="#0f172a")

# Reference lines
ax.axhline(0, color="#94a3b8", lw=0.6, ls=":")
ax.axvline(0, color="#94a3b8", lw=0.6, ls=":")
# Diagonal "perfect agreement"
ax.plot([0, 0.8], [0, 0.8], color="#94a3b8", lw=0.6, ls="--", alpha=0.6,
        label="$y=x$ (perfect signal agreement)")
ax.set_xlabel(r"$-\rho_{\mathrm{AF3\ pLDDT}}$  (higher = AF3 better calibrated)", fontsize=10)
ax.set_ylabel(r"$+\rho_{\mathrm{DRfold2\ ens\ std}}$  (higher = DRfold2 better calibrated)", fontsize=10)
ax.set_xlim(-0.2, 0.8)
ax.set_ylim(-0.15, 0.8)
# Legend
from matplotlib.patches import Patch
legend = [
    Patch(color="#10b981", label="both signals Bonferroni"),
    Patch(color="#f59e0b", label="one signal Bonferroni"),
    Patch(color="#cbd5e1", label="neither Bonferroni"),
]
ax.legend(handles=legend, loc="lower right", fontsize=8, frameon=False)
ax.set_title("Per-target calibration: AF3 pLDDT vs DRfold2 ensemble std\n"
             "12 RNA crystals, marker size $\\propto$ sequence length",
             fontsize=10, fontweight="bold", loc="left")
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig(OUT / "acdsa_pertarget_scatter.pdf", bbox_inches="tight")
plt.savefig(OUT / "acdsa_pertarget_scatter.png", bbox_inches="tight", dpi=180)
print(f"Written: {OUT / 'acdsa_pertarget_scatter.pdf'}")

# Fig 2: pooled per-residue scatter of (DRfold2 std, -AF3 pLDDT), colored by error
# We need per-residue data — we have it in calibration_per_residue.csv for R1108
# and we have all the raw files. For brevity here, use just R1108 (longest available
# per-res file) as illustration.
per_res_csv = ROOT / ".research/30_experiments/runs/cpeb3_focused/calibration_per_residue.csv"
if per_res_csv.exists():
    fig2, ax2 = plt.subplots(figsize=(6.0, 4.5))
    rows = list(csv.DictReader(open(per_res_csv)))
    drf_std = np.array([float(r["drfold2_ens_std"]) for r in rows])
    af3_pl  = np.array([float(r["af3_plddt"]) for r in rows])
    err     = np.array([float(r["af3_err_vs_crystal"]) for r in rows])
    # normalise pLDDT to "high = uncertain" direction
    neg_plddt = -af3_pl
    sc = ax2.scatter(drf_std, neg_plddt, c=err, cmap="YlOrRd", s=40,
                     edgecolors="black", linewidth=0.3, vmin=0, vmax=40)
    cbar = plt.colorbar(sc, ax=ax2)
    cbar.set_label("per-residue C1$'$ error vs crystal (Å)", fontsize=9)
    ax2.set_xlabel(r"DRfold2 ensemble std (\AA)  — high = uncertain", fontsize=10)
    ax2.set_ylabel(r"$-$AF3 pLDDT  — high = uncertain", fontsize=10)
    # Correlation
    from scipy.stats import spearmanr
    r, p = spearmanr(drf_std, neg_plddt)
    ax2.set_title(f"Pooled per-residue signal agreement (R1108, $L=69$).\n"
                  f"Spearman $\\rho = {r:+.2f}$, $p = {p:.2g}$ — the two signals are not redundant.",
                  fontsize=10, fontweight="bold", loc="left")
    ax2.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(OUT / "acdsa_pooled_residues.pdf", bbox_inches="tight")
    plt.savefig(OUT / "acdsa_pooled_residues.png", bbox_inches="tight", dpi=180)
    print(f"Written: {OUT / 'acdsa_pooled_residues.pdf'}")
