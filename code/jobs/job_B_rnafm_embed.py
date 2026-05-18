#!/usr/bin/env python3
"""Job B: compute RNA-FM embeddings for all 868 Stanford sequences. GPU heavy.
Uses multimolecule/rnafm via HuggingFace transformers (cleanest install)."""
import sys, time, os
from pathlib import Path
import pandas as pd, numpy as np
import torch

ROOT = Path("/workspace/rna3d")
OUT = ROOT / ".research/30_experiments/runs/B_rnafm_embeddings"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "embeddings").mkdir(exist_ok=True)
LOG = OUT / "log.txt"

def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

# Try multimolecule (HF-hosted RNA-FM)
try:
    from multimolecule import RnaTokenizer, RnaFmModel
    log("multimolecule available")
except ImportError:
    log("Installing multimolecule...")
    os.system("pip install -q multimolecule 2>&1 | tail -3")
    from multimolecule import RnaTokenizer, RnaFmModel

log("Loading RNA-FM tokenizer + model...")
try:
    tokenizer = RnaTokenizer.from_pretrained("multimolecule/rnafm")
    model = RnaFmModel.from_pretrained("multimolecule/rnafm")
    model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
        log(f"Model on cuda. Mem: {torch.cuda.memory_allocated()/1e9:.1f} GB")
except Exception as e:
    log(f"FAILED to load RNA-FM: {e}")
    sys.exit(1)

# Load sequences
csvs = []
for split, fname in [("train", "train_sequences.csv"), ("val", "validation_sequences.csv"), ("test", "test_sequences.csv")]:
    df = pd.read_csv(ROOT / "data/kaggle_raw" / fname)
    df["split"] = split
    csvs.append(df)
all_seqs = pd.concat(csvs, ignore_index=True)
all_seqs["length"] = all_seqs["sequence"].str.len()
# Sort by length ascending (process small first, OOM resilience)
all_seqs = all_seqs.sort_values("length").reset_index(drop=True)
log(f"Processing {len(all_seqs)} sequences (lengths {all_seqs.length.min()}-{all_seqs.length.max()})")

results = []
MAX_LEN = 1022  # RNA-FM context limit (1024 with special tokens)
t0 = time.time()
for i, row in all_seqs.iterrows():
    target_id = row["target_id"]
    seq = row["sequence"][:MAX_LEN]  # truncate if too long
    truncated = len(row["sequence"]) > MAX_LEN
    try:
        with torch.no_grad():
            inputs = tokenizer(seq, return_tensors="pt")
            if torch.cuda.is_available():
                inputs = {k: v.cuda() for k, v in inputs.items()}
            out = model(**inputs, output_hidden_states=False)
            # Last hidden state: [1, L, D]
            last = out.last_hidden_state.squeeze(0).cpu().numpy()  # [L, D]
            # Save per-sequence pooled (mean) embedding + per-residue
            pooled = last.mean(axis=0)
            np.savez(OUT / "embeddings" / f"{target_id}.npz", per_residue=last, pooled=pooled)
        results.append({"target_id": target_id, "split": row["split"], "length": len(seq), "truncated": truncated, "embed_dim": last.shape[1], "status": "OK"})
        if (i + 1) % 25 == 0:
            log(f"OK {i+1}/{len(all_seqs)} ({target_id} len={len(seq)} mem={torch.cuda.memory_allocated()/1e9:.1f}GB) elapsed={time.time()-t0:.0f}s")
            pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)
    except torch.cuda.OutOfMemoryError:
        torch.cuda.empty_cache()
        results.append({"target_id": target_id, "split": row["split"], "length": len(seq), "truncated": truncated, "embed_dim": -1, "status": "OOM"})
        log(f"OOM {target_id} len={len(seq)}")
    except Exception as e:
        results.append({"target_id": target_id, "split": row["split"], "length": len(seq), "truncated": truncated, "embed_dim": -1, "status": f"ERR:{type(e).__name__}"})
        log(f"ERR {target_id} len={len(seq)} {e}")

pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)
ok = sum(1 for r in results if r["status"]=="OK")
log(f"DONE: {ok}/{len(all_seqs)} OK. Elapsed {time.time()-t0:.0f}s")

with open(OUT / "SUMMARY.md", "w") as f:
    f.write(f"# Job B — RNA-FM embeddings\n\n")
    f.write(f"- Model: multimolecule/rnafm (RNA-FM via HuggingFace)\n")
    f.write(f"- Total sequences: {len(all_seqs)}\n")
    f.write(f"- OK: {ok}\n- OOM: {sum(1 for r in results if r['status']=='OOM')}\n")
    f.write(f"- Truncated (>{MAX_LEN} nt): {sum(1 for r in results if r['truncated'])}\n")
    f.write(f"- Elapsed: {time.time()-t0:.0f}s\n")
    f.write(f"- Per-residue embeddings saved as .npz under embeddings/\n")
