# CPEB3 focused study — DRfold2 RESULTS (first run, 2026-05-18)

> Recovered from prior pod stdout (the previous script overwrote this with an
> empty file due to a path bug — looking in `relax/` instead of `folds/`).
> The bug is fixed in the new `scripts/run_cpeb3_full.sh`.

## Setup

- Predictor: **DRfold2** (PLOS Biology 2025) — ensemble of 4 configs cfg_95/96/97/99
- Targets:
  - R1107 human CPEB3 (pos 30 = **A**)
  - R1108 chimpanzee CPEB3 (pos 30 = **G**)
  - Reference: PDB 7QR3 chain C (chimp X-ray)

## RMSD vs 7QR3 (chimp X-ray ground truth)

| target | best model | RMSD chain C (A) | RMSD chain D (A) |
|---|---|---|---|
| R1107 human pred | opt_0_from_cfg_97_model_16.ret.pdb | **4.204** | 2.784 |
| R1108 chimp pred | opt_0_from_cfg_96_model_3.ret.pdb  | **3.726** | 2.913 |

## DISCRIMINATION test

DRfold2 Kabsch-aligned RMSD between R1107 and R1108 predictions: **1.858 A**.
Modest but real sensitivity to the single-nt difference.

## Per-residue prediction deltas (R1107 vs R1108)

- mean = 1.643, max = 4.210, position 30 (mutation site) = **2.148 A**
- Top-5 most-divergent residues: pos **9** (4.21), **60** (4.12), **51** (3.98), **22** (3.96), **24** (3.70)
- Residues 9 and 60 form the P1 helix endpoints — the region implicated by
  Skilandat 2016 in the slow-cleaving human variant.

## What's still missing

This was Run 1. NULL CONTROL (`run_cpeb3_full.sh` Phase 2) will rerun DRfold2 with
a peripheral single-nt mutation at position 41 (A->C, J3-P4 region) and compare
the top-5 divergent residues. If overlap is large -> F-001 is a model prior; if
small -> F-001 is real.

See `.research/40_findings/FINDINGS.md` for the full paper-style writeup.
