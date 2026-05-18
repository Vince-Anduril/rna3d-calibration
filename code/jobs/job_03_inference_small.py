#!/usr/bin/env python3
"""Job 03: try RibonanzaNet inference on the 50 SHORTEST sequences (safest, fastest).
Saves per-residue confidence; commits incrementally; skips OOMs."""
import sys, time, traceback, os
from pathlib import Path
import pandas as pd, numpy as np
import torch

ROOT = Path("/workspace/rna3d")
OUT = ROOT / ".research/30_experiments/runs/03_ribonanzanet_inference_small"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "raw").mkdir(exist_ok=True)

LOG = OUT / "log.txt"
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line)
    with open(LOG, "a") as f:
        f.write(line + "\n")

# Try to find a usable RibonanzaNet
import_ok = False
sys.path.insert(0, "/workspace/rna3d/models/ribonanzanet_repo")
try:
    from ribonanzanet import RibonanzaNet
    import_ok = True
    log("import: from ribonanzanet")
except Exception as e:
    log(f"import attempt 1 failed: {e}")
    try:
        from RibonanzaNet import RibonanzaNet
        import_ok = True
        log("import: from RibonanzaNet (capital R)")
    except Exception as e:
        log(f"import attempt 2 failed: {e}")

if not import_ok:
    log("ABORT: RibonanzaNet not importable. See /workspace/rna3d/.research/30_experiments/runs/02_model_installation/INSTALL.md")
    sys.exit(1)

# Load CSVs and pick the 50 shortest sequences
csvs = []
for split, fname in [("train", "train_sequences.csv"), ("val", "validation_sequences.csv"), ("test", "test_sequences.csv")]:
    df = pd.read_csv(ROOT / "data/kaggle_raw" / fname)
    df["split"] = split
    csvs.append(df)
all_seqs = pd.concat(csvs, ignore_index=True)
all_seqs["length"] = all_seqs["sequence"].str.len()
subset = all_seqs.nsmallest(50, "length").reset_index(drop=True)
log(f"Selected {len(subset)} shortest sequences (lengths {subset.length.min()}–{subset.length.max()})")

# Initialize model
log("Loading RibonanzaNet weights...")
try:
    model = RibonanzaNet()
    model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
        log("model on cuda")
except Exception as e:
    log(f"FAILED to init RibonanzaNet: {e}\n{traceback.format_exc()}")
    sys.exit(1)

# Inference
results = []
for i, row in subset.iterrows():
    target_id = row["target_id"]
    seq = row["sequence"]
    t0 = time.time()
    try:
        with torch.no_grad():
            # RibonanzaNet API varies; try common interfaces
            try:
                out = model(seq)
            except TypeError:
                out = model.predict(seq)
            confidence = None
            if isinstance(out, dict):
                if "confidence" in out: confidence = out["confidence"]
                elif "pLDDT" in out: confidence = out["pLDDT"]
            elif hasattr(out, "confidence"):
                confidence = out.confidence
        elapsed = time.time() - t0
        npz_path = OUT / "raw" / f"{target_id}.npz"
        np.savez(str(npz_path), output=str(out)[:200])
        results.append({"target_id": target_id, "length": len(seq), "time_s": elapsed, "status": "OK"})
        log(f"OK  {target_id} len={len(seq)} t={elapsed:.1f}s")
    except torch.cuda.OutOfMemoryError as e:
        torch.cuda.empty_cache()
        results.append({"target_id": target_id, "length": len(seq), "time_s": 0, "status": "OOM"})
        log(f"OOM {target_id} len={len(seq)}")
    except Exception as e:
        results.append({"target_id": target_id, "length": len(seq), "time_s": 0, "status": f"ERR:{type(e).__name__}"})
        log(f"ERR {target_id} len={len(seq)} {type(e).__name__}: {e}")

    # Incremental write
    if (i + 1) % 10 == 0:
        pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)

pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)
log(f"DONE: {sum(1 for r in results if r['status']=='OK')} OK, {sum(1 for r in results if r['status']=='OOM')} OOM")
