# Experiment F — Calibration & contact analysis

Generated: 2026-05-19T12:14:38.819845+00:00

This analysis uses

- the **7QR3 chain C** crystal (Przytula-Mally 2022) as ground truth,
- the **full 80-model DRfold2 ensemble** (4 configs × 20 models) for each of R1107 and R1108,
- the **5 AF3 seeds** with per-residue atom pLDDT.

## (1) Does the published P1.1↔P1 mechanism hold geometrically in the crystal?

Distance C1$'$(pos 30) → C1$'$(target residue) in 7QR3:C and 7QR3:D.

| residue | crystal C (Å) | crystal D (Å) | interpretation |
|---|---|---|---|
| 9 | 8.27 | 8.28 | **direct 3D contact** (P1.1↔P1) |
| 22 | 32.13 | 32.28 | non-local |
| 24 | 28.92 | 28.36 | non-local |
| 51 | 47.57 | 44.53 | non-local |
| 60 | 17.78 | 17.36 | weak coupling |

If residues 9 and 60 are <15 Å from residue 30 in the ground-truth crystal, the parsimonious geometric explanation of our cascade observation is **confirmed by the experimental structure**.

## (2) DRfold2 ensemble uncertainty

R1108 ensemble per-residue std: mean = **3.46 Å**, max = 6.78 Å (at residue 50).

R1107 ensemble per-residue std: mean = **3.81 Å**, max = 5.46 Å (at residue 50).

## (3) Calibration — does the model know where it is uncertain?

Correlation between per-residue uncertainty and per-residue error vs the 7QR3:C crystal.

| signal | Pearson r | Spearman ρ | p (Spearman) | interpretation |
|---|---|---|---|---|
| DRfold2 ensemble std (R1108) | +0.333 | +0.217 | 0.0735 | positive = calibrated |
| AF3 pLDDT (R1108 best seed) | -0.339 | -0.268 | 0.0258 | negative = calibrated |
| AF3 pLDDT (R1108 default seed) | -0.400 | -0.452 | 9.74e-05 | negative = calibrated |
| DRfold2 std vs (-AF3 pLDDT) | -0.225 | -0.237 | 0.0503 | positive = cross-model agreement on hard residues |

## (4) Per-residue summary at cascade positions

| res | crystal d→30 (Å) | DRfold2 err vs crystal (Å) | DRfold2 ens std (Å) | AF3 err vs crystal (Å) | AF3 pLDDT |
|---|---|---|---|---|---|
| 9 | 8.27 | 2.92 | 1.79 | 7.94 | 22.2 |
| 22 | 32.13 | 8.09 | 4.50 | 26.44 | 18.7 |
| 24 | 28.92 | 2.58 | 4.16 | 17.87 | 20.2 |
| 30 | 0.00 | 3.48 | 2.05 | 16.68 | 27.2 |
| 51 | 47.57 | 5.66 | 6.48 | 14.70 | 29.4 |
| 60 | 17.78 | 2.60 | 2.74 | 4.06 | 23.8 |

## (5) Interpretation

- AF3's pLDDT shows weak/no significant correlation with per-residue crystal error (Spearman -0.27, p=0.026). At pLDDT~25 the AF3 prediction is uniformly low-confidence and the signal does not differentiate well across residues.
- DRfold2 ensemble std does not significantly correlate with per-residue crystal error (Spearman +0.22, p=0.074).

