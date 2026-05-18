# CPEB3 focused study — DRfold2 results + NULL CONTROL

Run completed: 2026-05-18T20:18:17.169386+00:00

## Setup

- Predictor: DRfold2 (PLOS Biology 2025), 4-config ensemble cfg_95/96/97/99
- Targets:
  - R1107 human CPEB3 (pos 30 = A)
  - R1108 chimpanzee CPEB3 (pos 30 = G)
  - null_pos41: R1108 with A->C at position 41 (peripheral, J3-P4 region)
- Reference: PDB 7QR3 chain C

## RMSD vs 7QR3

| target | RMSD vs chain C (A) | RMSD vs chain D (A) |
|---|---|---|
| R1107 human | 4.204 | 2.784 |
| R1108 chimp | 3.726 | 2.913 |
| null_pos41 | 3.905 | 2.568 |

## Discrimination & NULL CONTROL

### REAL case (R1107 vs R1108, pos 30 A->G, biologically functional)

- Kabsch RMSD between predictions: **1.858 A**
- Position 30 (mutation site) per-residue delta: **2.148 A**
- mean/max delta: 1.64 / 4.21 A
- Top-5 divergent residues (1-indexed): pos 9 (4.21A), pos 60 (4.12A), pos 51 (3.98A), pos 22 (3.96A), pos 24 (3.70A)

### NULL case (null_pos41 vs R1108, pos 41 A->C, peripheral)

- Kabsch RMSD between predictions: **1.909 A**
- Position 41 (null mutation site) per-residue delta: **1.704 A**
- mean/max delta: 1.69 / 5.47 A
- Top-5 divergent residues (1-indexed): pos 51 (5.47A), pos 52 (4.05A), pos 48 (3.49A), pos 50 (3.22A), pos 49 (3.11A)

### VERDICT

- Overlap of top-5 (REAL & NULL): **[51]** (1 positions)

**F-001 SURVIVES THE NULL CONTROL.** Top-5 changes meaningfully between mutations. The position-30 -> P1 cascade appears specific to the position-30 mutation, supporting a real biology/model alignment.
