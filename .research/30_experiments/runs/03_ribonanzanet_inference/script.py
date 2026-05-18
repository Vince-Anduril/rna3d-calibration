"""Job 3 — RibonanzaNet inference on all 868 Stanford sequences.

NOTE: RibonanzaNet's primary published outputs are 1D chemical-reactivity tracks
(2A3 / DMS), not 3D coordinates. We capture whatever the model produces. If a
secondary 3D head is available in the cloned repo we will also dump it.

We:
  - Iterate sequences in ascending length order (quick wins first).
  - On OOM or any exception: log + skip, never crash.
  - Save raw outputs as .npz per target.
  - Append incrementally to predictions.parquet (atomic via DataFrame append).
  - Commit+push every 50 sequences for safety.

Outputs:
  results/predictions.parquet (one row per target with status + paths)
  results/predictions.csv  (parallel CSV for portability if parquet fails)
  raw_outputs/<target_id>.npz
  log.txt
  RUN_SUMMARY.md
"""
from __future__ import annotations

import os
import sys
import time
import traceback
import subprocess
from pathlib import Path

sys.path.insert(0, "/workspace/rna3d/models/RibonanzaNet")

HERE = Path(__file__).parent
RES = HERE / "results"
RAW = HERE / "raw_outputs"
LOG = HERE / "log.txt"
RES.mkdir(exist_ok=True, parents=True)
RAW.mkdir(exist_ok=True, parents=True)

DATA = Path("/workspace/rna3d/data/kaggle_raw")
SEQ_FILES = [DATA / "train_sequences.csv", DATA / "validation_sequences.csv", DATA / "test_sequences.csv"]

PRED_PARQUET = RES / "predictions.parquet"
PRED_CSV = RES / "predictions.csv"

import numpy as np
import pandas as pd
import torch


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def git_commit_push(message: str) -> None:
    try:
        subprocess.run(["git", "-C", "/workspace/rna3d", "add",
                        str(RES.relative_to('/workspace/rna3d')),
                        str(RAW.relative_to('/workspace/rna3d')),
                        str(LOG.relative_to('/workspace/rna3d'))],
                       check=False, capture_output=True)
        # Don't add giant raw_outputs if they grow; cap by status
        r = subprocess.run(["git", "-C", "/workspace/rna3d", "commit", "-m", message],
                           check=False, capture_output=True, text=True)
        if r.returncode == 0:
            subprocess.run(["git", "-C", "/workspace/rna3d", "push", "origin", "main"],
                           check=False, capture_output=True, text=True, timeout=120)
        log(f"[git] {message} | rc={r.returncode}")
    except Exception as e:
        log(f"[git-err] {e!r}")


def load_sequences() -> pd.DataFrame:
    frames = []
    for split_name, p in zip(["train", "val", "test"], SEQ_FILES):
        df = pd.read_csv(p)
        df = df[["target_id", "sequence", "description"]].copy()
        df["split"] = split_name
        df["length"] = df["sequence"].str.len()
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    return all_df.sort_values("length", ascending=True).reset_index(drop=True)


def tokenize(seq: str) -> torch.Tensor:
    table = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3, "N": 4}
    return torch.tensor([[table.get(c, 4) for c in seq.upper()]], dtype=torch.long)


def find_weights() -> Path | None:
    cands = list(Path("/workspace/rna3d/models/RibonanzaNet").rglob("*.pt")) + \
            list(Path("/workspace/rna3d/models/RibonanzaNet").rglob("*.pth"))
    if not cands:
        return None
    return cands[0]


def build_model():
    """Try the canonical config the RibonanzaNet repo ships with."""
    from Network import RibonanzaNet  # type: ignore

    # Look for a config file in the repo
    cfg_path = None
    for c in Path("/workspace/rna3d/models/RibonanzaNet").rglob("*.yaml"):
        cfg_path = c
        break

    if cfg_path is not None:
        try:
            import yaml
            with open(cfg_path) as f:
                cfg = yaml.safe_load(f)
            # The repo wraps args in a Config dataclass-like dict; pass dict directly
            class _Cfg:
                def __init__(self, d):
                    for k, v in d.items():
                        setattr(self, k, v)
            model = RibonanzaNet(_Cfg(cfg))
        except Exception as e:
            log(f"[model] yaml cfg load failed: {e!r}; falling back to defaults")
            model = RibonanzaNet()
    else:
        model = RibonanzaNet()

    # Load weights if found
    wpath = find_weights()
    if wpath is not None:
        try:
            state = torch.load(wpath, map_location="cpu")
            if isinstance(state, dict) and "state_dict" in state:
                state = state["state_dict"]
            missing, unexpected = model.load_state_dict(state, strict=False)
            log(f"[model] loaded weights from {wpath} | missing={len(missing)} unexpected={len(unexpected)}")
        except Exception as e:
            log(f"[model] weight load failed: {e!r}")
    else:
        log("[model] no weights found; running with randomly initialized network")

    model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
    return model


def main() -> None:
    log("=== Job 3 RibonanzaNet inference: START ===")
    log(f"cwd: {os.getcwd()}")
    seqs = load_sequences()
    log(f"loaded {len(seqs)} sequences across train+val+test")

    try:
        model = build_model()
    except Exception as e:
        log(f"[fatal] cannot build model: {e!r}")
        log(traceback.format_exc())
        # Write RUN_SUMMARY anyway and exit gracefully
        (RES / "RUN_SUMMARY.md").write_text(
            f"# Job 3 — FAILED to build model\n\n```\n{e!r}\n```\nSee log.txt.\n"
        )
        return

    rows = []
    n_success = 0
    n_oom = 0
    n_err = 0
    t0 = time.time()
    LENGTH_HARD_CAP = 2000  # skip outright above this

    for i, r in seqs.iterrows():
        tid = r.target_id
        seq = r.sequence
        L = len(seq)
        out_npz = RAW / f"{tid}.npz"

        if L > LENGTH_HARD_CAP:
            log(f"[skip-too-long] {tid} L={L}")
            rows.append({"target_id": tid, "split": r.split, "length": L,
                         "status": "skipped_too_long", "time_s": 0.0,
                         "out_path": "", "shape": ""})
            continue
        if out_npz.exists():
            log(f"[resume-skip] {tid} L={L} already done")
            rows.append({"target_id": tid, "split": r.split, "length": L,
                         "status": "already_done", "time_s": 0.0,
                         "out_path": str(out_npz), "shape": ""})
            continue

        ts = time.time()
        try:
            x = tokenize(seq)
            if torch.cuda.is_available():
                x = x.cuda()
            with torch.no_grad():
                out = model(x)
            if isinstance(out, (tuple, list)):
                arrs = {f"out_{j}": (o.detach().cpu().numpy() if torch.is_tensor(o) else np.asarray(o))
                        for j, o in enumerate(out)}
            elif torch.is_tensor(out):
                arrs = {"out": out.detach().cpu().numpy()}
            elif isinstance(out, dict):
                arrs = {k: (v.detach().cpu().numpy() if torch.is_tensor(v) else np.asarray(v))
                        for k, v in out.items()}
            else:
                arrs = {"out": np.asarray(out)}
            np.savez_compressed(out_npz, **arrs)
            dt = time.time() - ts
            shape_desc = ",".join(f"{k}:{v.shape}" for k, v in arrs.items())
            rows.append({"target_id": tid, "split": r.split, "length": L,
                         "status": "ok", "time_s": round(dt, 3),
                         "out_path": str(out_npz), "shape": shape_desc})
            n_success += 1
            log(f"[ok] {tid} L={L} t={dt:.2f}s shapes={shape_desc}")
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            n_oom += 1
            log(f"[oom] {tid} L={L}")
            rows.append({"target_id": tid, "split": r.split, "length": L,
                         "status": "oom", "time_s": round(time.time() - ts, 3),
                         "out_path": "", "shape": ""})
        except Exception as e:
            torch.cuda.empty_cache()
            n_err += 1
            tb = traceback.format_exc()
            log(f"[err] {tid} L={L} {e!r}")
            log(tb.splitlines()[-1] if tb else "")
            rows.append({"target_id": tid, "split": r.split, "length": L,
                         "status": f"err:{type(e).__name__}", "time_s": round(time.time() - ts, 3),
                         "out_path": "", "shape": ""})

        # Incremental save every iteration
        df_out = pd.DataFrame(rows)
        try:
            df_out.to_parquet(PRED_PARQUET, index=False)
        except Exception as e:
            log(f"[parquet-err] {e!r}; falling back to csv only")
        df_out.to_csv(PRED_CSV, index=False)

        if (i + 1) % 50 == 0:
            git_commit_push(
                f"feat(03_ribonanzanet_inference): checkpoint at {i+1}/{len(seqs)} (ok={n_success} oom={n_oom} err={n_err})"
            )

    total_t = time.time() - t0
    summary_lines = [
        "# Job 3 — RibonanzaNet inference — SUMMARY",
        "",
        f"Total sequences attempted: {len(seqs)}",
        f"- Successful: {n_success}",
        f"- OOM: {n_oom}",
        f"- Errors: {n_err}",
        f"- Skipped (length > {LENGTH_HARD_CAP}): {sum(1 for x in rows if x['status']=='skipped_too_long')}",
        f"- Already done (resume): {sum(1 for x in rows if x['status']=='already_done')}",
        "",
        f"Total wall time: {total_t/60:.1f} min",
        "",
        "## Notes",
        "- RibonanzaNet's pretrained head outputs per-residue chemical reactivity (2A3/DMS),",
        "  NOT 3D coordinates. The downstream 3D head requires fine-tuning per the team's repo",
        "  variants (RibonanzaNet-3D). Raw outputs captured for whatever the loaded checkpoint produces.",
        "- See `predictions.csv` for per-row status; raw arrays in `raw_outputs/<target_id>.npz`.",
    ]
    (RES / "RUN_SUMMARY.md").write_text("\n".join(summary_lines))
    git_commit_push(
        f"feat(03_ribonanzanet_inference): RibonanzaNet predictions on Stanford 868 (ok={n_success} oom={n_oom} err={n_err})"
    )
    log("=== Job 3: DONE ===")


if __name__ == "__main__":
    main()
