"""Focused inference pipeline for CPEB3 ribozyme variants (R1107 human, R1108 chimp).

The plan (see docs/superpowers/specs/FOCUSED_CPEB3_PLAN.md):
    1. Run each available 3D RNA predictor on both sequences.
    2. Save per-residue (x,y,z) coordinates + any per-residue confidence score.
    3. A separate `focused_cpeb3_analyze.py` script computes RMSD vs the chimp
       crystal (PDB 7QR3), per-residue accuracy, ECE in the P1/P1.1 window.

This file just defines the sequences, the predictor wrappers, and the I/O contract.
Each predictor wrapper is best-effort — if a model fails to import, we log it and
move on. The downstream analysis works with whatever predictors succeeded.
"""
from __future__ import annotations
import os, sys, time, traceback
from pathlib import Path
import numpy as np
import torch

# ---------------------------------------------------------------------------
# Targets (exact sequences from CASP15 R1107 / R1108; verified against
# https://predictioncenter.org/casp15/target.cgi)
# ---------------------------------------------------------------------------
TARGETS = {
    "R1107_human": {
        "sequence": "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU",
        "organism": "Homo sapiens",
        "position_30_base": "A",
        "notes": "Slow-cleaving CPEB3 ribozyme. P1/P1.1 mispairing hypothesis explains lower activity.",
    },
    "R1108_chimp": {
        "sequence": "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU",
        "organism": "Pan troglodytes",
        "position_30_base": "G",
        "notes": "Fast-cleaving CPEB3 ribozyme (~4x faster than human). PDB 7QR3 is the crystal structure.",
    },
}

# Sanity check: the only difference between the two sequences is at position 30 (1-indexed).
def _assert_single_nt_difference():
    h = TARGETS["R1107_human"]["sequence"]
    c = TARGETS["R1108_chimp"]["sequence"]
    assert len(h) == len(c) == 69, f"unexpected length: {len(h)} vs {len(c)}"
    diffs = [i for i in range(len(h)) if h[i] != c[i]]
    assert diffs == [29], f"expected difference at 1-indexed position 30 (0-indexed 29), got {[i+1 for i in diffs]}"
    assert h[29] == "A" and c[29] == "G", f"expected A/G at pos 30, got {h[29]}/{c[29]}"
_assert_single_nt_difference()


OUT_ROOT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
LOG_FILE = OUT_ROOT / "run.log"

def _ensure_out():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

def log(msg: str):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    try:
        _ensure_out()
        with open(LOG_FILE, "a") as f:
            f.write(line + "\n")
    except OSError:
        # /workspace not present (Mac import-time sanity check), log to stdout only
        pass


# ---------------------------------------------------------------------------
# Predictor wrappers
#
# Each wrapper takes (target_id, sequence) and returns a dict:
#     {
#         "coords": np.ndarray of shape (L, 3),       # per-residue centroid (C1' or similar)
#         "confidence": np.ndarray of shape (L,) or None,
#         "predictor": "name",
#         "raw_path": Path to the raw predictor output (PDB/CIF),
#     }
# Or raises and the caller logs+skips.
# ---------------------------------------------------------------------------

def predict_with_drfold2(target_id: str, sequence: str) -> dict:
    """DRfold2 wrapper. Assumes DRfold2 is installed and reachable via the
    `DRfold2.predict` API or the `drfold2` CLI. See `focused_cpeb3_install.sh`."""
    out_dir = OUT_ROOT / target_id / "drfold2"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Try CLI first
    fasta_path = out_dir / "input.fasta"
    fasta_path.write_text(f">{target_id}\n{sequence}\n")
    pdb_out = out_dir / "model.pdb"
    # Placeholder — actual CLI args depend on DRfold2 install
    raise NotImplementedError("DRfold2 wrapper to be wired during install (see TODO in focused_cpeb3_install.sh)")


def predict_with_rhofold(target_id: str, sequence: str) -> dict:
    """RhoFold+ wrapper. Assumes the rhofold repo is cloned at
    /workspace/rna3d/models/rhofold and weights are downloaded."""
    out_dir = OUT_ROOT / target_id / "rhofold"
    out_dir.mkdir(parents=True, exist_ok=True)
    raise NotImplementedError("RhoFold+ wrapper to be wired during install")


def predict_with_af3_server(target_id: str, sequence: str) -> dict:
    """Manual fallback — submit to https://alphafoldserver.com/ and download the
    JSON job result. The wrapper just reads a pre-existing JSON file from disk."""
    json_path = OUT_ROOT / target_id / "af3_server.json"
    if not json_path.exists():
        raise FileNotFoundError(
            f"Manual AF3-server submission required for {target_id}. "
            f"Save the AlphaFold Server JSON output to {json_path}."
        )
    import json
    data = json.loads(json_path.read_text())
    # AF3 server format may differ — wire to actual schema once we have a sample
    raise NotImplementedError("AF3 server JSON parser — wire once we have a real sample")


PREDICTORS = {
    "drfold2": predict_with_drfold2,
    "rhofold": predict_with_rhofold,
    "af3_server": predict_with_af3_server,
}


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictors", nargs="+", default=list(PREDICTORS.keys()),
                    help="Which predictors to try.")
    args = ap.parse_args()

    log(f"=== focused_cpeb3 starting, predictors={args.predictors} ===")
    for target_id, target_data in TARGETS.items():
        seq = target_data["sequence"]
        log(f"--- target {target_id} ({len(seq)} nt, organism {target_data['organism']}) ---")
        for pname in args.predictors:
            fn = PREDICTORS[pname]
            try:
                t0 = time.time()
                result = fn(target_id, seq)
                elapsed = time.time() - t0
                # Save standardized output
                out = OUT_ROOT / target_id / pname
                out.mkdir(parents=True, exist_ok=True)
                np.save(out / "coords.npy", result["coords"])
                if result.get("confidence") is not None:
                    np.save(out / "confidence.npy", result["confidence"])
                with open(out / "meta.txt", "w") as f:
                    f.write(f"predictor={pname}\n")
                    f.write(f"target_id={target_id}\n")
                    f.write(f"sequence_length={len(seq)}\n")
                    f.write(f"elapsed_s={elapsed:.1f}\n")
                    f.write(f"raw_path={result.get('raw_path', '')}\n")
                log(f"  {pname}: OK ({elapsed:.1f}s) -> {out}")
            except NotImplementedError as e:
                log(f"  {pname}: NOT WIRED YET — {e}")
            except FileNotFoundError as e:
                log(f"  {pname}: missing input — {e}")
            except Exception as e:
                log(f"  {pname}: FAILED ({type(e).__name__}) — {e}")
                log(traceback.format_exc())
    log("=== focused_cpeb3 DONE ===")


if __name__ == "__main__":
    main()
