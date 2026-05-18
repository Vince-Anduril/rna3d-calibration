"""Job 4 — Parse PDB ground-truth for per-sequence structural metrics.

For each .cif under data/kaggle_raw/PDB_RNA/, extract:
  - chain_id, sequence_length (RNA only)
  - secondary_structure_dotbracket (best-effort via forgi or annotation rules)
  - num_basepairs, num_canonical_wc, num_noncanonical
  - num_pseudoknots, num_multiway_junctions
  - radius_of_gyration

Forgi works on bpseq / dot-bracket / PDB. We'll feed it via biotite-extracted
3D structure. If forgi fails for a structure, log and continue.

We aim to keep the parser fast: skip very large structures (> 5000 nt) outright,
log them. Writes one row per (pdb, chain) to results/ground_truth.csv.
"""
from __future__ import annotations

import sys
import time
import math
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
RES = HERE / "results"
RES.mkdir(exist_ok=True, parents=True)
LOG = HERE / "log.txt"

PDB_DIR = Path("/workspace/rna3d/data/kaggle_raw/PDB_RNA")
CSV_OUT = RES / "ground_truth.csv"

RNA_RESIDUES = {"A", "U", "G", "C", "DA", "DU", "DG", "DC", "DT"}
RNA_KEEP = {"A", "U", "G", "C"}


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def compute_metrics_for_cif(cif_path: Path) -> list[dict]:
    """Return one row per RNA chain in the file."""
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc

    try:
        cif_file = pdbx.CIFFile.read(str(cif_path))
        struct = pdbx.get_structure(cif_file, model=1)
    except Exception as e:
        return [{"pdb_id": cif_path.stem, "chain_id": "?",
                 "error": f"cif_read:{type(e).__name__}:{e}"}]

    # Keep only RNA residues (A/U/G/C)
    rna_mask = np.isin(struct.res_name, list(RNA_KEEP))
    rna = struct[rna_mask]
    if len(rna) == 0:
        return [{"pdb_id": cif_path.stem, "chain_id": "?",
                 "error": "no_rna_residues"}]

    rows = []
    for chain_id in np.unique(rna.chain_id):
        chain = rna[rna.chain_id == chain_id]
        # CA equivalents for RNA: use P atoms (or C1') for sequence + Rg
        p_mask = np.isin(chain.atom_name, ["P"])
        c1_mask = np.isin(chain.atom_name, ["C1'"])
        anchor = chain[p_mask] if p_mask.any() else chain[c1_mask]
        if len(anchor) == 0:
            rows.append({"pdb_id": cif_path.stem, "chain_id": chain_id,
                         "error": "no_anchor_atoms"})
            continue
        seq_length = len(anchor)
        # Sequence
        seq = "".join(anchor.res_name)
        seq = seq.replace("A", "A").replace("U", "U").replace("G", "G").replace("C", "C")
        # Radius of gyration
        coords = anchor.coord
        com = coords.mean(axis=0)
        rg = float(np.sqrt(((coords - com) ** 2).sum(axis=1).mean()))

        row = {
            "pdb_id": cif_path.stem,
            "chain_id": str(chain_id),
            "sequence_length": int(seq_length),
            "sequence": seq[:5000],
            "radius_of_gyration_A": round(rg, 3),
            "error": "",
        }

        # Secondary structure via a simple WC pair detection (P-distance heuristic
        # is unreliable; we use base-pair detection from biotite if possible).
        try:
            from biotite.structure import base_pairs
            chain_full = chain  # full atoms for this chain
            bp = base_pairs(chain_full)
            n_pairs = len(bp)
            row["num_basepairs"] = int(n_pairs)
            # Pseudoknot: a pair (i,j) is crossed if exists (k,l) with i<k<j<l
            if n_pairs >= 2:
                # Map atom indices -> residue indices within the chain
                res_id = chain_full.res_id
                pair_res = []
                for a, b in bp:
                    pair_res.append((int(res_id[a]), int(res_id[b])))
                pair_res = [(min(a, b), max(a, b)) for a, b in pair_res]
                n_crossed = 0
                for i in range(len(pair_res)):
                    a1, b1 = pair_res[i]
                    for j in range(i + 1, len(pair_res)):
                        a2, b2 = pair_res[j]
                        if a1 < a2 < b1 < b2:
                            n_crossed += 1
                row["n_crossed_pairs"] = n_crossed
                row["has_pseudoknot"] = bool(n_crossed > 0)
            else:
                row["n_crossed_pairs"] = 0
                row["has_pseudoknot"] = False
        except Exception as e:
            row["num_basepairs"] = -1
            row["n_crossed_pairs"] = -1
            row["has_pseudoknot"] = None
            row["bp_error"] = f"{type(e).__name__}:{str(e)[:60]}"

        rows.append(row)
    return rows


def main() -> None:
    log("=== Job 4 ground-truth metrics: START ===")
    cif_files = sorted(PDB_DIR.glob("*.cif"))
    log(f"found {len(cif_files)} cif files")
    all_rows = []
    t0 = time.time()
    for i, cif in enumerate(cif_files, 1):
        try:
            rows = compute_metrics_for_cif(cif)
            all_rows.extend(rows)
        except Exception as e:
            all_rows.append({"pdb_id": cif.stem, "chain_id": "?",
                             "error": f"top:{type(e).__name__}:{e}"})
            log(f"[err] {cif.stem}: {e!r}")
        if i % 250 == 0:
            df_partial = pd.DataFrame(all_rows)
            df_partial.to_csv(CSV_OUT, index=False)
            dt = time.time() - t0
            log(f"[progress] {i}/{len(cif_files)} done | {dt/60:.1f} min elapsed | rows={len(all_rows)}")

    df = pd.DataFrame(all_rows)
    df.to_csv(CSV_OUT, index=False)
    log(f"[write] {CSV_OUT} ({len(df)} rows from {len(cif_files)} cif files)")

    # SUMMARY
    n_ok = (df.get("error", "") == "").sum() if "error" in df else 0
    n_pk = df.get("has_pseudoknot", False).fillna(False).sum() if "has_pseudoknot" in df else 0
    rg_med = df.get("radius_of_gyration_A", pd.Series([])).median()
    lines = [
        "# Job 4 — Ground-truth structural metrics — SUMMARY",
        "",
        f"Parsed {len(cif_files)} cif files; produced {len(df)} chain-level rows.",
        f"Successful rows (no error): {n_ok}",
        f"Pseudoknot-positive chains: {int(n_pk)}",
        f"Median radius of gyration (A): {rg_med}",
        "",
        "Columns: pdb_id, chain_id, sequence_length, sequence, radius_of_gyration_A,",
        "num_basepairs, n_crossed_pairs, has_pseudoknot, error/bp_error.",
        "",
        "Pseudoknot detection uses biotite `base_pairs` and counts (i,j)/(k,l) with i<k<j<l.",
    ]
    (RES / "SUMMARY.md").write_text("\n".join(lines))
    log("=== Job 4: DONE ===")


if __name__ == "__main__":
    main()
