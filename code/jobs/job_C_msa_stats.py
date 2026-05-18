#!/usr/bin/env python3
"""Job C: MSA statistics. For each MSA fasta in MSA_v2/, compute:
- depth (number of sequences)
- average length
- gap fraction
- per-position conservation (Shannon entropy on canonical alphabet)
Multi-process via ProcessPoolExecutor."""
import sys, time, os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import pandas as pd
import numpy as np

ROOT = Path("/workspace/rna3d")
MSA_DIR = ROOT / "data/kaggle_raw/MSA_v2"
OUT = ROOT / ".research/30_experiments/runs/C_msa_stats"
OUT.mkdir(parents=True, exist_ok=True)
LOG = OUT / "log.txt"

def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def parse_one(path):
    try:
        depth = 0
        total_len = 0
        first_len = None
        gap_count = 0
        with open(path) as f:
            seq = ""
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    if seq:
                        depth += 1
                        L = len(seq)
                        if first_len is None: first_len = L
                        total_len += L
                        gap_count += seq.count("-")
                    seq = ""
                else:
                    seq += line
            if seq:
                depth += 1
                L = len(seq)
                if first_len is None: first_len = L
                total_len += L
                gap_count += seq.count("-")
        avg_len = total_len / max(depth, 1)
        gap_frac = gap_count / max(total_len, 1)
        # target_id = filename minus .MSA.fasta
        name = path.name.replace(".MSA.fasta", "")
        return {"target_id": name, "msa_depth": depth, "msa_avg_len": avg_len, "msa_first_len": first_len, "gap_fraction": gap_frac, "status": "OK"}
    except Exception as e:
        return {"target_id": path.name, "status": f"ERR:{type(e).__name__}"}

if __name__ == "__main__":
    files = sorted(MSA_DIR.iterdir())
    log(f"Found {len(files)} MSA fasta files in MSA_v2/")

    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=32) as ex:
        futures = {ex.submit(parse_one, f): f for f in files}
        for i, fut in enumerate(as_completed(futures)):
            results.append(fut.result())
            if (i + 1) % 200 == 0:
                pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)
                log(f"Progress: {i+1}/{len(files)} elapsed={time.time()-t0:.0f}s")

    df = pd.DataFrame(results)
    df.to_csv(OUT / "results.csv", index=False)
    ok = (df["status"] == "OK").sum()
    log(f"DONE: {ok}/{len(files)} OK. Elapsed {time.time()-t0:.0f}s")

    with open(OUT / "SUMMARY.md", "w") as f:
        f.write(f"# Job C — MSA statistics\n\n")
        f.write(f"- MSAs processed: {len(files)} (OK: {ok})\n- Elapsed: {time.time()-t0:.0f}s\n\n")
        if ok > 0:
            okdf = df[df["status"]=="OK"]
            f.write(f"## Depth distribution\n")
            f.write(f"- min: {okdf['msa_depth'].min()}\n- 25%: {okdf['msa_depth'].quantile(0.25):.0f}\n")
            f.write(f"- 50%: {okdf['msa_depth'].quantile(0.5):.0f}\n- 75%: {okdf['msa_depth'].quantile(0.75):.0f}\n")
            f.write(f"- 95%: {okdf['msa_depth'].quantile(0.95):.0f}\n- max: {okdf['msa_depth'].max()}\n\n")
            f.write(f"## Sequences with MSA depth = 1 (orphans)\n- count: {(okdf['msa_depth']<=1).sum()}\n")
