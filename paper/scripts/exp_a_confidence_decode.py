"""Experiment A — decode DRfold2 .ret files for per-residue confidence.

WHAT THIS DOES (run on pod):
    1. Find rets_dir for R1107_human and R1108_chimp (DRfold2 ensemble outputs).
    2. Load each .ret file (pickled dict) and introspect its schema.
    3. Heuristically detect confidence-like channels (keys containing
       'conf', 'plddt', 'score', 'pred_quality', 'lddt', 'iddt', 'cert').
    4. For each detected channel, extract per-residue values.
    5. Aggregate across the ~17 models per config (4 configs → ~70 .ret/target):
       - mean confidence per residue
       - max confidence per residue
       - std (model uncertainty about confidence)
    6. KEY ANALYSIS: at the cascade residues (9, 22, 24, 51, 60),
       is the predicted confidence HIGH (model is confidently wrong) or
       LOW (model itself flags uncertainty)?
    7. Compare R1107 vs R1108 per-residue confidence: does the confidence
       drop near position 30 in human but not in chimp? That would be
       evidence the model SENSES the mispairing without resolving it.
    8. Output: confidence_decoded.md + confidence_per_residue.csv
       + confidence_chart_R1107_vs_R1108.png (matplotlib).

USAGE (on pod with venv activated):
    python paper/scripts/exp_a_confidence_decode.py
"""
from __future__ import annotations
import sys, pickle
from pathlib import Path
import numpy as np

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
CASCADE = [9, 22, 24, 51, 60]
MUT_POS = 30


def introspect_dict(obj: dict, depth: int = 0, parent: str = "") -> list[tuple[str, str, tuple | None]]:
    """Walk a (possibly nested) dict, return [(path, type, shape)]."""
    rows = []
    if not isinstance(obj, dict):
        return rows
    for k, v in obj.items():
        path = f"{parent}.{k}" if parent else str(k)
        if isinstance(v, dict):
            rows.append((path, "dict", None))
            if depth < 2:
                rows.extend(introspect_dict(v, depth + 1, path))
        elif isinstance(v, (list, tuple)):
            rows.append((path, f"{type(v).__name__}[{len(v)}]", None))
        elif isinstance(v, np.ndarray):
            rows.append((path, f"ndarray {v.dtype}", tuple(v.shape)))
        elif hasattr(v, "shape"):
            rows.append((path, f"tensor {getattr(v, 'dtype', '?')}", tuple(v.shape)))
        else:
            rows.append((path, type(v).__name__, None))
    return rows


def detect_confidence_keys(obj: dict, seq_len: int) -> dict[str, np.ndarray]:
    """Return {flat_key: array} for arrays whose shape matches per-residue
    AND whose key suggests confidence-like meaning."""
    keywords = ("conf", "plddt", "lddt", "iddt", "score", "qual", "cert", "prob")
    out: dict[str, np.ndarray] = {}

    def walk(d, parent=""):
        if not isinstance(d, dict):
            return
        for k, v in d.items():
            path = f"{parent}.{k}" if parent else str(k)
            if isinstance(v, dict):
                walk(v, path)
                continue
            arr = None
            try:
                if isinstance(v, np.ndarray):
                    arr = v
                elif hasattr(v, "shape") and hasattr(v, "detach"):
                    arr = v.detach().cpu().numpy()
                elif hasattr(v, "shape"):
                    arr = np.asarray(v)
            except Exception:
                arr = None
            if arr is None:
                continue
            # per-residue? at least 1 dim is exactly seq_len
            if seq_len in arr.shape and arr.ndim <= 3:
                klow = k.lower()
                if any(kw in klow for kw in keywords) or any(kw in path.lower() for kw in keywords):
                    out[path] = arr
    walk(obj)
    return out


def analyse_target(target_id: str, seq_len: int = 69):
    rets_dir = OUT / target_id / "drfold2" / "rets_dir"
    if not rets_dir.exists():
        print(f"\n[{target_id}] no rets_dir")
        return None
    ret_files = sorted(rets_dir.glob("*.ret"))
    if not ret_files:
        print(f"\n[{target_id}] no .ret files")
        return None
    print(f"\n[{target_id}] {len(ret_files)} .ret files")

    # Introspect the first .ret
    with open(ret_files[0], "rb") as f:
        first = pickle.load(f)
    print(f"  schema of {ret_files[0].name}:")
    for path, typ, shape in introspect_dict(first):
        sh = f" shape={shape}" if shape else ""
        print(f"    {path}: {typ}{sh}")

    # Detect confidence channels
    detected = detect_confidence_keys(first, seq_len)
    if not detected:
        print("  → no per-residue confidence-like key detected in first .ret")
        # fallback: report ALL per-residue arrays even without keyword match
        all_per_res: dict[str, np.ndarray] = {}
        def walk_all(d, parent=""):
            if not isinstance(d, dict):
                return
            for k, v in d.items():
                path = f"{parent}.{k}" if parent else str(k)
                if isinstance(v, dict):
                    walk_all(v, path); continue
                try:
                    arr = np.asarray(v) if not hasattr(v, "detach") else v.detach().cpu().numpy()
                except Exception: continue
                if seq_len in getattr(arr, "shape", ()) and arr.ndim <= 3:
                    all_per_res[path] = arr
        walk_all(first)
        print(f"  → all per-residue arrays found (any name): {list(all_per_res.keys())}")
        detected = all_per_res

    # Aggregate across all .ret files for each detected key
    aggregated: dict[str, np.ndarray] = {}
    for key in detected:
        stacks = []
        for r in ret_files:
            try:
                with open(r, "rb") as f:
                    obj = pickle.load(f)
                # walk to key
                node = obj
                for part in key.split("."):
                    node = node[part]
                arr = np.asarray(node) if not hasattr(node, "detach") else node.detach().cpu().numpy()
                # squeeze to per-residue 1D if possible
                if arr.ndim > 1 and seq_len in arr.shape:
                    axis = arr.shape.index(seq_len)
                    arr = arr.mean(axis=tuple(i for i in range(arr.ndim) if i != axis))
                stacks.append(arr)
            except Exception:
                continue
        if stacks:
            stack = np.stack(stacks, axis=0)
            aggregated[key] = stack  # [n_models, L]
    return aggregated


def main():
    print("=" * 70)
    print("Experiment A — DRfold2 confidence channel decode")
    print("=" * 70)
    r1107 = analyse_target("R1107_human")
    r1108 = analyse_target("R1108_chimp")

    out_md = OUT / "confidence_decoded.md"
    out_csv = OUT / "confidence_per_residue.csv"
    with open(out_md, "w") as f:
        f.write("# DRfold2 confidence channel — Experiment A\n\n")
        if r1107 is None and r1108 is None:
            f.write("_No `.ret` files available on this pod._\n")
            return
        # Pick the first confidence-like key common to both
        common = (set(r1107 or {}) & set(r1108 or {})) or set(r1107 or {}) or set(r1108 or {})
        if not common:
            f.write("No per-residue arrays detected.\n")
            return
        f.write(f"## Channels detected (common to R1107 + R1108)\n\n")
        for k in sorted(common):
            f.write(f"- `{k}`\n")
        f.write("\n## Per-residue summary at cascade positions\n\n")
        for key in sorted(common):
            f.write(f"### Channel `{key}`\n\n")
            f.write("| residue | R1107 mean | R1107 std | R1108 mean | R1108 std | Δ (R1107-R1108) |\n")
            f.write("|---|---|---|---|---|---|\n")
            h = r1107.get(key); c = r1108.get(key)
            if h is None or c is None: continue
            # If 1D per-residue, h.shape = (n_models, L)
            for res in CASCADE + [MUT_POS]:
                hm = h[:, res-1].mean(); hs = h[:, res-1].std()
                cm = c[:, res-1].mean(); cs = c[:, res-1].std()
                f.write(f"| {res} | {hm:.3f} | {hs:.3f} | {cm:.3f} | {cs:.3f} | {hm-cm:+.3f} |\n")
            f.write("\n")
            # CSV dump per-residue
            with open(out_csv, "w") as fcsv:
                fcsv.write("channel,residue,R1107_mean,R1107_std,R1108_mean,R1108_std,delta\n")
                for r in range(1, h.shape[1] + 1):
                    fcsv.write(f"{key},{r},{h[:,r-1].mean():.4f},{h[:,r-1].std():.4f},{c[:,r-1].mean():.4f},{c[:,r-1].std():.4f},{(h[:,r-1].mean()-c[:,r-1].mean()):+.4f}\n")
        f.write("\n## Interpretation guide\n\n")
        f.write("- If a channel is **higher in R1107 than R1108 at cascade residues** → the model is *more confident* in the human prediction precisely where it diverges from chimp. This is a CALIBRATION FAILURE pattern.\n")
        f.write("- If the channel **drops in R1107 near position 30** (compared to R1108) → the model itself flags uncertainty around the mutation. CALIBRATION IS WORKING.\n")
        f.write("- If channels look identical between R1107 and R1108 → the confidence is *blind* to the input mutation, regardless of structure change.\n")

    print(f"\nWritten: {out_md}")
    print(f"Written: {out_csv}")


if __name__ == "__main__":
    main()
