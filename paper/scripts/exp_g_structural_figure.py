"""Structural figure: DRfold2 predicted C1' trace overlaid on crystal, with
per-residue error encoded as point color.

Generates a 2-panel figure: one for a target where DRfold2 wins
(low RMSD, e.g. 9LKU at 2.06 Å), and one for the 12CI failure case to
show the multi-chain artifact contrast.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
OUT = ROOT / "paper/figures"
OUT.mkdir(parents=True, exist_ok=True)


def parse_cif_c1(path: Path, chain: str):
    lines = path.read_text().splitlines()
    cols, data_start = [], None
    in_header = False
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s); in_header = True
        elif in_header and not s.startswith("_atom_site."):
            data_start = i; break
    by_res = {}
    for ln in lines[data_start:]:
        if not ln.startswith("ATOM"): continue
        parts = ln.split()
        if len(parts) < len(cols): continue
        d = dict(zip(cols, parts))
        aid = d["_atom_site.label_atom_id"].strip('"')
        if aid != "C1'": continue
        if d["_atom_site.label_asym_id"] != chain: continue
        comp = d["_atom_site.label_comp_id"]
        if comp not in {"A", "U", "G", "C"}: continue
        try:
            ri = int(d["_atom_site.label_seq_id"])
            xyz = np.array([float(d["_atom_site.Cartn_x"]),
                            float(d["_atom_site.Cartn_y"]),
                            float(d["_atom_site.Cartn_z"])], dtype=np.float32)
        except (KeyError, ValueError):
            continue
        if ri not in by_res:
            by_res[ri] = xyz
    return np.stack([by_res[k] for k in sorted(by_res.keys())], axis=0)


def pdb_c1(path: Path):
    coords = []
    for ln in path.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "C1'":
            try:
                coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError: continue
    return np.array(coords, dtype=np.float32)


def kabsch_align(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p - pc, q - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ pp.T).T + qc


def plot_panel(ax, crystal, pred_aligned, title, rmsd):
    """Plot crystal trace (gray) + prediction colored by per-residue error."""
    err = np.linalg.norm(pred_aligned - crystal, axis=-1)
    # crystal trace
    ax.plot(crystal[:, 0], crystal[:, 1], crystal[:, 2],
            color="#94a3b8", linewidth=1.5, alpha=0.8, label="crystal C1'")
    ax.scatter(crystal[:, 0], crystal[:, 1], crystal[:, 2],
               color="#94a3b8", s=18, alpha=0.6)
    # prediction trace
    ax.plot(pred_aligned[:, 0], pred_aligned[:, 1], pred_aligned[:, 2],
            color="#0ea5e9", linewidth=1.5, alpha=0.7, label="DRfold2 prediction")
    sc = ax.scatter(pred_aligned[:, 0], pred_aligned[:, 1], pred_aligned[:, 2],
                    c=err, cmap="YlOrRd", s=24, vmin=0, vmax=max(8.0, err.max()),
                    edgecolors="black", linewidth=0.3)
    ax.set_title(f"{title}  (Kabsch RMSD = {rmsd:.2f} Å)", fontsize=10, fontweight="bold")
    ax.set_xlabel("X (Å)", fontsize=8)
    ax.set_ylabel("Y (Å)", fontsize=8)
    ax.tick_params(labelsize=6)
    ax.legend(loc="upper left", fontsize=7, frameon=False)
    return sc


def render(target_id: str, crystal_path: Path, chain: str, drf_pdb: Path,
           title: str):
    cry = parse_cif_c1(crystal_path, chain)
    drf = pdb_c1(drf_pdb)
    n = min(cry.shape[0], drf.shape[0])
    cry, drf = cry[:n], drf[:n]
    aligned = kabsch_align(drf, cry)
    rmsd = float(np.sqrt(((aligned - cry) ** 2).sum(-1).mean()))
    return cry, aligned, rmsd, title


# Pick two contrasting cases
A = render("9LKU", ROOT / "data/multi_rna/9LKU.cif", "A",
           sorted((ROOT / ".research/30_experiments/runs/multi_rna/9LKU/drfold2/folds").glob("opt_0_*.pdb"))[0],
           "9LKU: 2'-dG-III riboswitch (monomer, blind, $L$=63)")
B = render("12CI", ROOT / "data/multi_rna/12CI.cif", "A",
           sorted((ROOT / ".research/30_experiments/runs/multi_rna/12CI/drfold2/folds").glob("opt_0_*.pdb"))[0],
           "12CI: dopamine aptamer (\\textit{dimer crystal}, blind, $L$=82)")

fig = plt.figure(figsize=(11, 5))
ax1 = fig.add_subplot(121, projection="3d")
ax2 = fig.add_subplot(122, projection="3d")
sc1 = plot_panel(ax1, A[0], A[1], A[3], A[2])
sc2 = plot_panel(ax2, B[0], B[1], B[3], B[2])
fig.colorbar(sc2, ax=[ax1, ax2], shrink=0.6, pad=0.02,
             label="per-residue C1' deviation (Å)")
fig.suptitle("DRfold2 predicted backbone (blue) overlaid on experimental crystal (grey),\n"
             "Kabsch-aligned, points colored by per-residue C1' deviation.",
             fontsize=11, y=0.97)
out_pdf = OUT / "structural_overlay.pdf"
out_png = OUT / "structural_overlay.png"
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, bbox_inches="tight", dpi=180)
print(f"Written: {out_pdf}")
print(f"Written: {out_png}")
