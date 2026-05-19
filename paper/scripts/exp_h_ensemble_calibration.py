"""Experiment H — cross-target DRfold2 ensemble-std calibration.

For each of the 11 multi-RNA targets, load the 80-model DRfold2 ensemble
(C4' coords from rets_dir/*.pdb), compute per-residue ensemble standard
deviation, and correlate it with per-residue distance from the experimental
crystal (using the FINAL selected DRfold2 structure to align both).

Output: ensemble_calibration_results.{csv,md}.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
MULTI = ROOT / ".research/30_experiments/runs/multi_rna"
CPEB3 = ROOT / ".research/30_experiments/runs/cpeb3_focused/R1108_chimp"

TARGETS = [
    ("7QR3", ROOT / "data/crystal/7QR3.cif", "C", CPEB3 / "drfold2"),
    ("9DE7", ROOT / "data/multi_rna/9DE7.cif", "A", MULTI / "9DE7/drfold2"),
    ("9HRD", ROOT / "data/multi_rna/9HRD.cif", "A", MULTI / "9HRD/drfold2"),
    ("9HRF", ROOT / "data/multi_rna/9HRF.cif", "A", MULTI / "9HRF/drfold2"),
    ("9J4N", ROOT / "data/multi_rna/9J4N.cif", "A", MULTI / "9J4N/drfold2"),
    ("9E9O", ROOT / "data/multi_rna/9E9O.cif", "A", MULTI / "9E9O/drfold2"),
    ("9MFH", ROOT / "data/multi_rna/9MFH.cif", "A", MULTI / "9MFH/drfold2"),
    ("9LJN", ROOT / "data/multi_rna/9LJN.cif", "A", MULTI / "9LJN/drfold2"),
    ("9LKE", ROOT / "data/multi_rna/9LKE.cif", "A", MULTI / "9LKE/drfold2"),
    ("9LKU", ROOT / "data/multi_rna/9LKU.cif", "A", MULTI / "9LKU/drfold2"),
    ("9UW0", ROOT / "data/multi_rna/9UW0.cif", "A", MULTI / "9UW0/drfold2"),
    ("12CI", ROOT / "data/multi_rna/12CI.cif", "A", MULTI / "12CI/drfold2"),
]


def parse_cif_c4(path: Path, chain: str):
    lines = path.read_text().splitlines()
    cols, data_start = [], None
    in_header = False
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s); in_header = True
        elif in_header and not s.startswith("_atom_site."):
            data_start = i; break
    by_res: dict[int, np.ndarray] = {}
    for ln in lines[data_start:]:
        if not ln.startswith("ATOM"): continue
        parts = ln.split()
        if len(parts) < len(cols): continue
        d = dict(zip(cols, parts))
        aid = d["_atom_site.label_atom_id"].strip('"')
        if aid != "C4'": continue
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
    return np.stack([by_res[k] for k in sorted(by_res.keys())], axis=0) if by_res else None


def pdb_c4(path: Path):
    coords = []
    for ln in path.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "C4'":
            try:
                coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError: continue
    return np.array(coords, dtype=np.float32) if coords else None


def kabsch_align_to(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p - pc, q - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (R @ pp.T).T + qc


def kabsch_per_res(p, q):
    if p is None or q is None or p.shape != q.shape: return None
    aligned = kabsch_align_to(p, q)
    return np.linalg.norm(aligned - q, axis=-1)


def load_ensemble(drf_dir: Path, ref_len: int):
    rets = drf_dir / "rets_dir"
    pdbs = sorted(rets.glob("*.pdb"))
    coords = []
    for p in pdbs:
        c = pdb_c4(p)
        if c is not None and c.shape[0] == ref_len:
            coords.append(c)
    return np.stack(coords, axis=0) if coords else None


def ensemble_std(ens):
    """Per-residue std after Kabsch-aligning each model to model 0."""
    ref = ens[0]
    aligned = np.stack([kabsch_align_to(m, ref) for m in ens], axis=0)
    mean_pos = aligned.mean(0)
    dev = np.linalg.norm(aligned - mean_pos[None, :, :], axis=-1)
    return dev.std(0), dev.mean(0)


def main():
    print("=" * 70)
    print("Experiment H — DRfold2 ensemble-std calibration across 11 targets")
    print("=" * 70)
    rows = []
    for tid, crystal_path, chain, drf_dir in TARGETS:
        cry = parse_cif_c4(crystal_path, chain)
        if cry is None:
            print(f"{tid}: no crystal C4'"); continue
        ens = load_ensemble(drf_dir, cry.shape[0])
        if ens is None:
            print(f"{tid}: ensemble length mismatch (crystal {cry.shape[0]})"); continue
        std, _ = ensemble_std(ens)
        # final structure (opt_0)
        final_pdbs = sorted((drf_dir / "folds").glob("opt_0_*.pdb"))
        if not final_pdbs:
            print(f"{tid}: no final fold"); continue
        final = pdb_c4(final_pdbs[0])
        if final is None or final.shape[0] != cry.shape[0]:
            print(f"{tid}: final fold length mismatch"); continue
        err = kabsch_per_res(final, cry)
        # Spearman: high ensemble std should correlate with high error
        rho, p = spearmanr(std, err)
        rows.append({
            "id": tid, "L": int(cry.shape[0]),
            "n_models": int(ens.shape[0]),
            "ens_std_mean": float(std.mean()), "ens_std_max": float(std.max()),
            "err_mean": float(err.mean()), "err_max": float(err.max()),
            "spearman_rho": float(rho), "spearman_p": float(p),
        })
        sig = "**" if p < 0.05 else ""
        bonf = "BONF" if p < 0.05 / 11 else ""
        print(f"  {tid:5} L={cry.shape[0]:3} | err mean={err.mean():.2f} | std mean={std.mean():.2f} | ρ={rho:+.3f}{sig} p={p:.3g} {bonf}")

    # Aggregate
    import datetime
    out_csv = MULTI / "ensemble_calibration_results.csv"
    out_md = MULTI / "ensemble_calibration_results.md"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(out_md, "w") as f:
        f.write("# Experiment H — DRfold2 ensemble-std calibration\n\n")
        f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
        f.write("Per-target Spearman rank correlation between DRfold2 80-model "
                "C4' ensemble standard deviation and per-residue error of the "
                "final selected structure against the crystal.\n\n")
        f.write("| PDB | L | n_models | err mean (Å) | std mean (Å) | Spearman ρ | p | Bonferroni (α=0.05/11=0.0045) |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            bf = "**yes**" if r["spearman_p"] < 0.05/11 else "no"
            sig = "**" if r["spearman_p"] < 0.05 else ""
            f.write(f"| {r['id']} | {r['L']} | {r['n_models']} | {r['err_mean']:.2f} | {r['ens_std_mean']:.2f} | "
                    f"{r['spearman_rho']:+.3f}{sig} | {r['spearman_p']:.3g} | {bf} |\n")
        n_sig = sum(1 for r in rows if r["spearman_p"] < 0.05)
        n_bonf = sum(1 for r in rows if r["spearman_p"] < 0.05/11)
        n_pos = sum(1 for r in rows if r["spearman_rho"] > 0)
        f.write(f"\nSummary: {n_pos}/{len(rows)} positive correlations (high std → high error). "
                f"{n_sig}/{len(rows)} significant at α=0.05 raw. "
                f"{n_bonf}/{len(rows)} after Bonferroni (α=0.05/{len(rows)}).\n")
    print(f"\nWritten: {out_csv}")
    print(f"Written: {out_md}")


if __name__ == "__main__":
    main()
