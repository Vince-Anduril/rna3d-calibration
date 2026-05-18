"""Analyze the saved focused-CPEB3 inference results.

Inputs (per predictor that succeeded):
    /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/<target>/<predictor>/coords.npy
    /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/<target>/<predictor>/confidence.npy   (optional)

Reference (chimp ground truth):
    /workspace/rna3d/data/kaggle_raw/PDB_RNA/7qr3.cif

Outputs:
    /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/SUMMARY.md
    /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/figures/*.png

Key analyses:
    A. Per-predictor RMSD vs chimp crystal (R1108 only; we have no R1107 crystal in CASP15 release).
    B. Confidence-vs-error scatter, per residue, both variants.
    C. Position-30 window zoom (residues 25-35) — does the predictor distinguish the variants?
    D. ECE bootstrap CI on confidence, focusing on the P1/P1.1 region (TBD residue range from secondary structure).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
PDB_REF = Path("/workspace/rna3d/data/kaggle_raw/PDB_RNA/7qr3.cif")
SUMMARY = ROOT / "SUMMARY.md"
FIGS = ROOT / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# Reuse the existing calibration helpers in this repo
sys.path.insert(0, str(Path(__file__).parent))
from calibration import (
    expected_calibration_error, brier_score, errors_to_correctness, bootstrap_ci,
    reliability_table,
)


def load_reference_coords() -> np.ndarray | None:
    """Load chimp CPEB3 (7QR3) C1' coordinates, chain C or D (RNA chains)."""
    try:
        import biotite.structure.io as bsio
        s = bsio.load_structure(str(PDB_REF))
        # Take chain C first (per RCSB summary: chains C,D are the RNA pair)
        for chain in ("C", "D"):
            mask = (s.chain_id == chain) & (s.atom_name == "C1'")
            if mask.any():
                return np.asarray(s.coord[mask], dtype=np.float32)
        return None
    except Exception as e:
        print(f"failed to load reference: {e}")
        return None


def kabsch_rmsd(pred: np.ndarray, true: np.ndarray) -> float:
    """Per-sample RMSD after Kabsch alignment. NaN if shapes mismatch."""
    if pred.shape[0] != true.shape[0] or pred.shape[0] < 3:
        return float("nan")
    pc, qc = pred.mean(0), true.mean(0)
    pp, qq = pred - pc, true - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    aligned = (R @ pp.T).T + qc
    return float(np.sqrt(((aligned - true) ** 2).sum(-1).mean()))


def per_residue_error(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    """Per-residue distance after Kabsch alignment."""
    pc, qc = pred.mean(0), true.mean(0)
    pp, qq = pred - pc, true - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    aligned = (R @ pp.T).T + qc
    return np.linalg.norm(aligned - true, axis=-1)


def main():
    ref = load_reference_coords()
    if ref is None:
        print("No reference structure available; skipping reference-based metrics.")
        ref_len = None
    else:
        ref_len = ref.shape[0]
        print(f"Loaded reference (chimp 7QR3): {ref_len} residues")

    targets = [d for d in ROOT.iterdir() if d.is_dir() and not d.name.startswith(".")]
    targets = [d for d in targets if d.name in ("R1107_human", "R1108_chimp")]

    rows = []
    for target_dir in targets:
        target_id = target_dir.name
        for pred_dir in [d for d in target_dir.iterdir() if d.is_dir()]:
            pname = pred_dir.name
            coords_path = pred_dir / "coords.npy"
            if not coords_path.exists():
                continue
            coords = np.load(coords_path)
            conf_path = pred_dir / "confidence.npy"
            conf = np.load(conf_path) if conf_path.exists() else None

            rmsd = float("nan")
            err = None
            # Reference comparison only valid for chimp R1108 (matches 7QR3)
            if ref is not None and target_id == "R1108_chimp" and coords.shape[0] == ref.shape[0]:
                rmsd = kabsch_rmsd(coords, ref)
                err = per_residue_error(coords, ref)

            row = {
                "target_id": target_id,
                "predictor": pname,
                "n_residues": int(coords.shape[0]),
                "rmsd_vs_7qr3": rmsd,
                "has_confidence": conf is not None,
            }
            rows.append(row)

            # Per-residue analysis
            if err is not None and conf is not None and len(err) == len(conf):
                # ECE if confidence in [0,1]; threshold 4 Å as "correct"
                c = np.asarray(conf, dtype=float)
                if c.min() < 0 or c.max() > 1:
                    # Normalize to [0,1] via min-max (best-effort; pLDDT-like scores are
                    # typically already 0-100; we divide by 100)
                    if c.max() > 1.5 and c.max() < 110:
                        c = c / 100.0
                    else:
                        c = (c - c.min()) / (c.max() - c.min() + 1e-9)
                correct = errors_to_correctness(err, threshold_A=4.0)
                ece, ece_lo, ece_hi = bootstrap_ci(
                    lambda cc, yy: expected_calibration_error(cc, yy, n_bins=8),
                    c, correct, n_boot=300, seed=42,
                )
                row.update({"ece_4A": ece, "ece_4A_lo": ece_lo, "ece_4A_hi": ece_hi})
                row.update({"brier_4A": brier_score(c, correct)})
                # P1/P1.1 window: residues 1-15 and 65-69 form the P1 helix
                # (5'-end pairs with 3'-end). The mispairing site involves residue 30
                # which is in the P3-P1.1 junction region.
                # We zoom on residues 25-35 (1-indexed → 24:35 in 0-indexed).
                lo, hi = 24, 35
                if hi <= len(c):
                    row.update({
                        "mean_conf_pos25_35": float(c[lo:hi].mean()),
                        "mean_err_pos25_35": float(err[lo:hi].mean()),
                    })

    # Write summary
    with open(SUMMARY, "w") as f:
        f.write("# CPEB3 focused study — analysis SUMMARY\n\n")
        f.write("## Per-predictor metrics\n\n")
        if rows:
            keys = sorted({k for r in rows for k in r.keys()})
            f.write("| " + " | ".join(keys) + " |\n")
            f.write("|" + "|".join(["---"] * len(keys)) + "|\n")
            for r in rows:
                f.write("| " + " | ".join(str(r.get(k, "")) for k in keys) + " |\n")
        else:
            f.write("_No predictor results available._\n")
        f.write("\n## Reference structure\n\n")
        f.write(f"Source: PDB 7QR3 (chimp CPEB3, 2.18 Å X-ray crystal). "
                f"Residue count loaded: {ref_len if ref_len is not None else 'N/A'}.\n")
        f.write("\n## Discrimination analysis (does any predictor distinguish R1107 from R1108?)\n\n")
        # Pair up by predictor name
        by_pred: dict[str, dict] = {}
        for r in rows:
            by_pred.setdefault(r["predictor"], {})[r["target_id"]] = r
        for p, pair in by_pred.items():
            if "R1107_human" in pair and "R1108_chimp" in pair:
                hc = pair["R1107_human"]
                cc = pair["R1108_chimp"]
                f.write(f"### {p}\n\n")
                f.write(f"- mean_conf_pos25_35 (human): {hc.get('mean_conf_pos25_35', 'N/A')}\n")
                f.write(f"- mean_conf_pos25_35 (chimp): {cc.get('mean_conf_pos25_35', 'N/A')}\n")
                if "mean_conf_pos25_35" in hc and "mean_conf_pos25_35" in cc:
                    d = cc["mean_conf_pos25_35"] - hc["mean_conf_pos25_35"]
                    f.write(f"- delta (chimp - human): **{d:+.4f}**\n")
                    f.write(f"  → if |delta| < 0.02 the predictor does NOT discriminate the single-nt difference.\n")
                f.write("\n")
    print(f"Wrote {SUMMARY}")


if __name__ == "__main__":
    main()
