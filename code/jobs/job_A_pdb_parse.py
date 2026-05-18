#!/usr/bin/env python3
"""Job A: parse all PDB files in data/kaggle_raw/PDB_RNA/, extract per-structure metrics.
Uses biotite (fast). Multi-process via concurrent.futures."""
import sys, time, os, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

ROOT = Path("/workspace/rna3d")
PDB_DIR = ROOT / "data/kaggle_raw/PDB_RNA"
OUT = ROOT / ".research/30_experiments/runs/A_pdb_parse"
OUT.mkdir(parents=True, exist_ok=True)
LOG = OUT / "log.txt"

def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

import biotite.structure.io as bsio
import biotite.structure as bs

def parse_one(path):
    try:
        # biotite handles .pdb and .cif transparently
        try:
            structure = bsio.load_structure(str(path))
        except Exception:
            return {"file": path.name, "status": f"LOAD_FAIL"}

        # take only RNA atoms (RA, RU, RG, RC) or any atom
        if hasattr(structure, "res_name"):
            res_names = structure.res_name
            rna_mask = np.isin(res_names, ["A", "U", "G", "C", "RA", "RU", "RG", "RC", "DA", "DT", "DG", "DC"])
            structure = structure[rna_mask]

        n_atoms = len(structure)
        if n_atoms == 0:
            return {"file": path.name, "n_atoms": 0, "status": "EMPTY"}

        # Sequence length = number of unique (chain, res_id) pairs
        coords = structure.coord
        chains = structure.chain_id
        res_ids = structure.res_id
        residue_keys = list(zip(chains.tolist(), res_ids.tolist()))
        n_residues = len(set(residue_keys))

        # Radius of gyration from C atoms
        c_mask = (structure.atom_name == "C1'") | (structure.atom_name == "P")
        c_coords = coords[c_mask]
        if len(c_coords) > 0:
            centroid = c_coords.mean(axis=0)
            rg = float(np.sqrt(((c_coords - centroid)**2).sum(axis=1).mean()))
        else:
            rg = -1

        # Chains
        n_chains = len(set(chains.tolist()))

        return {
            "file": path.name,
            "n_atoms": int(n_atoms),
            "n_residues": int(n_residues),
            "n_chains": int(n_chains),
            "radius_of_gyration": rg,
            "status": "OK",
        }
    except Exception as e:
        return {"file": path.name, "status": f"ERR:{type(e).__name__}"}

if __name__ == "__main__":
    files = sorted(PDB_DIR.iterdir())
    log(f"Found {len(files)} files in PDB_RNA/")

    results = []
    t0 = time.time()
    n_workers = min(64, len(files))  # 128 cores available, leave room
    log(f"Spawning {n_workers} workers")

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = {ex.submit(parse_one, f): f for f in files}
        for i, fut in enumerate(as_completed(futures)):
            results.append(fut.result())
            if (i + 1) % 500 == 0:
                pd.DataFrame(results).to_csv(OUT / "results.csv", index=False)
                log(f"Progress: {i+1}/{len(files)} ({(i+1)/len(files)*100:.1f}%) elapsed={time.time()-t0:.0f}s")

    df = pd.DataFrame(results)
    df.to_csv(OUT / "results.csv", index=False)
    ok = (df["status"] == "OK").sum()
    log(f"DONE: {ok}/{len(files)} OK. Elapsed {time.time()-t0:.0f}s")

    with open(OUT / "SUMMARY.md", "w") as f:
        f.write(f"# Job A — PDB ground-truth parse\n\n")
        f.write(f"- Total files: {len(files)}\n- OK: {ok}\n- Errors: {len(df)-ok}\n")
        f.write(f"- Elapsed: {time.time()-t0:.0f}s\n\n")
        if ok > 0:
            okdf = df[df["status"]=="OK"]
            f.write(f"## Stats on OK structures\n\n")
            f.write(f"- median residue count: {okdf['n_residues'].median():.0f}\n")
            f.write(f"- mean residue count: {okdf['n_residues'].mean():.0f}\n")
            f.write(f"- median R_g: {okdf['radius_of_gyration'].median():.2f}\n")
            f.write(f"- residue count percentiles: 25%={okdf['n_residues'].quantile(0.25):.0f}, ")
            f.write(f"50%={okdf['n_residues'].quantile(0.5):.0f}, 75%={okdf['n_residues'].quantile(0.75):.0f}, 95%={okdf['n_residues'].quantile(0.95):.0f}\n")
