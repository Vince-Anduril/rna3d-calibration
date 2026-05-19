# Experiment H — DRfold2 ensemble-std calibration

Generated: 2026-05-19T19:00:09.793309+00:00

Per-target Spearman rank correlation between DRfold2 80-model C4' ensemble standard deviation and per-residue error of the final selected structure against the crystal.

| PDB | L | n_models | err mean (Å) | std mean (Å) | Spearman ρ | p | Bonferroni (α=0.05/11=0.0045) |
|---|---|---|---|---|---|---|---|
| 7QR3 | 69 | 80 | 3.46 | 3.46 | +0.217 | 0.0735 | no |
| 9DE7 | 57 | 80 | 2.60 | 0.82 | +0.493** | 9.73e-05 | **yes** |
| 9HRD | 67 | 80 | 4.71 | 3.31 | +0.183 | 0.139 | no |
| 9HRF | 70 | 80 | 6.82 | 4.90 | +0.331** | 0.00516 | no |
| 9J4N | 81 | 80 | 2.89 | 0.56 | +0.664** | 1.38e-11 | **yes** |
| 9E9O | 101 | 80 | 3.99 | 3.01 | +0.430** | 7.38e-06 | **yes** |
| 9MFH | 76 | 80 | 6.38 | 4.03 | +0.017 | 0.883 | no |
| 9LJN | 71 | 80 | 3.23 | 0.67 | +0.581** | 1.1e-07 | **yes** |
| 9LKE | 70 | 80 | 3.32 | 0.58 | +0.598** | 4.54e-08 | **yes** |
| 9LKU | 63 | 80 | 1.61 | 0.86 | +0.312** | 0.0127 | no |
| 9UW0 | 63 | 80 | 1.78 | 0.86 | +0.254** | 0.0442 | no |
| 12CI | 82 | 80 | 15.45 | 4.70 | +0.428** | 6.12e-05 | **yes** |

Summary: 12/12 positive correlations (high std → high error). 9/12 significant at α=0.05 raw. 6/12 after Bonferroni (α=0.05/12).
