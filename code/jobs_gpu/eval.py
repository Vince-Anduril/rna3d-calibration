"""Evaluation: per-category metrics + reliability diagrams + ribozyme focus deep-dive.

Loads a trained checkpoint, runs inference on val+test, stratifies metrics by
the v2 categorization, produces figures + a SUMMARY.md ready for FINDINGS.md
ingestion by rna-scientist.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).parent))
from model import RNAStructureModel, ModelConfig
from data import StanfordRNADataset, encode, PAD_ID
from loss import _kabsch_align
from calibration import (
    expected_calibration_error, maximum_calibration_error, brier_score,
    errors_to_correctness, bootstrap_ci, reliability_table,
)


@torch.no_grad()
def kabsch_rmsd(pred: torch.Tensor, true: torch.Tensor, valid: torch.Tensor) -> float:
    """Per-sample RMSD after Kabsch alignment. Returns NaN if <3 valid residues."""
    if valid.sum() < 3:
        return float("nan")
    p, q = pred[valid], true[valid]
    pc, qc = p.mean(0), q.mean(0)
    pp, qq = p - pc, q - qc
    H = pp.T @ qq
    U, S, Vt = torch.linalg.svd(H)
    d = torch.sign(torch.det(Vt.T @ U.T))
    D = torch.diag(torch.tensor([1.0, 1.0, d], device=p.device, dtype=p.dtype))
    R = Vt.T @ D @ U.T
    p_aligned = (R @ pp.T).T + qc
    rmsd = torch.sqrt(((p_aligned - q) ** 2).sum(-1).mean()).item()
    return rmsd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--data-root", default="/workspace/rna3d/data")
    ap.add_argument("--category-csv", default="/workspace/rna3d/.research/30_experiments/runs/01b_description_categorization_v2/results/category_assignment.csv")
    ap.add_argument("--out-dir", default="/workspace/rna3d/.research/30_experiments/runs/GPU_training_eval")
    ap.add_argument("--max-len", type=int, default=256)
    ap.add_argument("--cpu", action="store_true")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() and not args.cpu else "cpu")
    print(f"Device: {device}")

    cfg = ModelConfig(max_len=args.max_len)
    model = RNAStructureModel(cfg).to(device)
    ck = torch.load(args.ckpt, map_location=device, weights_only=False)
    model.load_state_dict(ck["state_dict"])
    model.eval()
    print(f"Loaded checkpoint: step={ck.get('step', '?')}, params={model.num_params/1e6:.2f}M")

    ds = StanfordRNADataset(data_root=Path(args.data_root), splits=("val", "test"),
                            max_len=args.max_len, category_csv=Path(args.category_csv))

    # Per-sample inference + RMSD + per-residue confidence
    rows = []
    per_res_records = []
    for i in range(len(ds)):
        s = ds[i]
        ids = s["token_ids"].unsqueeze(0).to(device)
        coords_true = s["coords"].to(device)
        valid = s["valid"].to(device)
        with torch.no_grad():
            out_ = model(ids)
        pred_coords = out_["coords"].squeeze(0)
        conf = out_["confidence"].squeeze(0)

        rmsd = kabsch_rmsd(pred_coords, coords_true, valid)
        rows.append({
            "target_id": s["target_id"], "split": s["split"], "category": s["category"],
            "length": s["length"], "n_valid": int(valid.sum().item()),
            "rmsd": rmsd, "mean_conf": float(conf[valid].mean().item()) if valid.any() else float("nan"),
        })
        # Per-residue: confidence vs |pred - true| after global Kabsch
        if valid.sum() >= 3:
            aligned_full = _kabsch_align(pred_coords.unsqueeze(0), coords_true.unsqueeze(0), valid.unsqueeze(0)).squeeze(0)
            err = torch.norm(aligned_full - coords_true, dim=-1)
            for j in torch.where(valid)[0].tolist():
                per_res_records.append({
                    "target_id": s["target_id"], "split": s["split"], "category": s["category"],
                    "residue_idx": j, "confidence": float(conf[j].item()), "error": float(err[j].item()),
                })

    per_seq = pd.DataFrame(rows)
    per_res = pd.DataFrame(per_res_records)
    per_seq.to_csv(out / "per_sequence.csv", index=False)
    per_res.to_csv(out / "per_residue.csv", index=False)

    # Defensive: if no residues had valid ground truth, write a clear empty SUMMARY and exit.
    if len(per_res) == 0 or "category" not in per_res.columns:
        with open(out / "SUMMARY.md", "w") as f:
            f.write("# GPU training — evaluation (EMPTY)\n\n")
            f.write(f"Checkpoint: `{args.ckpt}`\n\n")
            f.write("**No residues had valid ground-truth coordinates.** Usually means the\n")
            f.write("PDB file lookup failed (filenames/chain IDs don't match). Inspect\n")
            f.write("`per_sequence.csv` for `n_valid==0` rows and check `data.py`.\n\n")
            f.write(f"- per_sequence rows: {len(per_seq)}\n- per_residue rows: {len(per_res)}\n")
            if len(per_seq) > 0:
                f.write(f"- mean n_valid per sequence: {per_seq['n_valid'].mean():.1f}\n")
        print(f"Eval wrote empty SUMMARY (no per-residue ground truth). Check {out}/SUMMARY.md")
        return

    # Per-category aggregate
    cat_agg = per_seq.groupby("category").agg(
        n=("rmsd", "size"),
        rmsd_median=("rmsd", "median"),
        rmsd_mean=("rmsd", "mean"),
        conf_mean=("mean_conf", "mean"),
    ).round(3).sort_values("n", ascending=False)
    cat_agg.to_csv(out / "per_category.csv")

    # Calibration metrics per category — threshold "correct" at 5 Å
    ERR_THRESHOLD = 5.0
    cal_rows = []
    for cat in ["all"] + sorted(per_res["category"].dropna().unique().tolist()):
        sub = per_res if cat == "all" else per_res[per_res["category"] == cat]
        if len(sub) < 10:
            continue
        confs = sub["confidence"].values
        correct = errors_to_correctness(sub["error"].values, threshold_A=ERR_THRESHOLD)
        ece, ece_lo, ece_hi = bootstrap_ci(
            lambda c, y: expected_calibration_error(c, y, n_bins=10),
            confs, correct, n_boot=500,
        )
        mce = maximum_calibration_error(confs, correct, n_bins=10)
        brier = brier_score(confs, correct)
        cal_rows.append({
            "category": cat, "n_residues": len(sub),
            "ece": ece, "ece_lo": ece_lo, "ece_hi": ece_hi,
            "mce": mce, "brier": brier, "mean_conf": float(confs.mean()),
            "mean_err": float(sub["error"].mean()),
            "frac_correct_5A": float(correct.mean()),
        })
    cal_df = pd.DataFrame(cal_rows).round(4)
    cal_df.to_csv(out / "calibration_per_category.csv", index=False)

    # Reliability diagrams: global + ribozyme + top 2 other classes by n
    plot_cats = ["all"]
    if "ribozyme" in cal_df["category"].values:
        plot_cats.append("ribozyme")
    other = cal_df[~cal_df["category"].isin(["all", "ribozyme"])].nlargest(2, "n_residues")["category"].tolist()
    plot_cats.extend(other)

    fig, axes = plt.subplots(1, len(plot_cats), figsize=(5 * len(plot_cats), 5), squeeze=False)
    for ax, cat in zip(axes[0], plot_cats):
        sub = per_res if cat == "all" else per_res[per_res["category"] == cat]
        confs = sub["confidence"].values
        correct = errors_to_correctness(sub["error"].values, threshold_A=ERR_THRESHOLD)
        rt = reliability_table(confs, correct, n_bins=10)
        xs, ys, ns = [], [], []
        for r in rt["rows"]:
            if r["count"] > 0:
                xs.append(r["mean_confidence"])
                ys.append(r["mean_correctness"])
                ns.append(r["count"])
        ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="perfect calibration")
        if xs:
            ax.scatter(xs, ys, s=[max(20, n / 2) for n in ns], alpha=0.7, c="steelblue", edgecolors="navy")
            for x, y, n in zip(xs, ys, ns):
                ax.annotate(f"n={n}", (x, y), fontsize=7, alpha=0.7, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("predicted confidence (binned)")
        ax.set_ylabel(f"frac correct (err <= {ERR_THRESHOLD} Å)")
        ece_str = ""
        if cat in cal_df["category"].values:
            row = cal_df[cal_df["category"] == cat].iloc[0]
            ece_str = f"\nECE = {row['ece']:.3f} [{row['ece_lo']:.3f}, {row['ece_hi']:.3f}]"
        ax.set_title(f"{cat} (n={len(sub)} residues){ece_str}", fontsize=10)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
    plt.suptitle("Reliability diagrams — confidence vs. fraction of correct predictions", y=1.02)
    plt.tight_layout()
    plt.savefig(out / "reliability_diagrams.png", dpi=130, bbox_inches="tight")
    plt.close()

    # SUMMARY
    with open(out / "SUMMARY.md", "w") as f:
        f.write("# GPU training — evaluation\n\n")
        f.write(f"Checkpoint: `{args.ckpt}` (step {ck.get('step', '?')})\n\n")
        f.write(f"## Per-category accuracy (val + test)\n\n")
        f.write(cat_agg.to_markdown())
        f.write("\n\n## Per-category calibration\n\n")
        f.write(f"Threshold for 'correct': error <= {ERR_THRESHOLD} Å. Bootstrapped 95% CI on ECE (n_boot=500).\n\n")
        f.write(cal_df.to_markdown(index=False))
        f.write("\n\n## Reliability diagram\n\n![](reliability_diagrams.png)\n\n")
        f.write("## Ribozyme focus (H-001 evaluation)\n\n")
        if "ribozyme" in cal_df["category"].values:
            r = cal_df[cal_df["category"] == "ribozyme"].iloc[0]
            g = cal_df[cal_df["category"] == "all"].iloc[0]
            ratio = r["ece"] / g["ece"] if g["ece"] > 0 else float("inf")
            f.write(f"- n_residues (ribozyme): {int(r['n_residues'])}\n")
            f.write(f"- ECE (ribozyme): {r['ece']:.4f} [{r['ece_lo']:.4f}, {r['ece_hi']:.4f}]\n")
            f.write(f"- ECE (global):   {g['ece']:.4f} [{g['ece_lo']:.4f}, {g['ece_hi']:.4f}]\n")
            f.write(f"- Ratio ribozyme/global: **{ratio:.2f}x**\n")
            f.write(f"- H-001 hypothesis (ribozyme/global ≥ 2x): {'**SUPPORTED**' if ratio >= 2.0 else ('**FALSIFIED**' if ratio < 1.2 else '**INCONCLUSIVE** (band 1.2-2.0)')}\n")
        else:
            f.write("_No ribozyme sequences in val+test (cannot evaluate H-001 directly)._\n")
        f.write("\n## Per-sequence detail\n\n")
        f.write(per_seq.to_markdown(index=False))

    print(f"Eval done. Outputs in {out}")


if __name__ == "__main__":
    main()
