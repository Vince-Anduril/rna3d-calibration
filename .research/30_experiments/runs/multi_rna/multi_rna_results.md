# Cross-family analysis — Experiment G

Generated: 2026-05-19T20:04:53.938862+00:00

All AF3 runs in this table are **single-chain** (one job per RNA).

## Ground-truth RMSD per target

| ID | RNA | family | L | deposit | DRfold2 seen? | AF3 seen? | DRfold2 (Å) | AF3 default (Å) | AF3 best (Å) |
|---|---|---|---|---|---|---|---|---|---|
| 7QR3 | CPEB3 ribozyme (R1108) | ribozyme | 69 | 2021-12 | yes | borderline | 3.73 | 6.54 | 5.35 |
| 9LJN | Guanine-II riboswitch | riboswitch | 71 | 2025-01 | no | no | 4.06 | 6.55 | 6.09 |
| 9UW0 | 2'-dG-III riboswitch | riboswitch | 63 | 2025-05 | no | no | 2.14 | 10.05 | 9.33 |
| 9HRD | Class V GTP aptamer (tetramer) | aptamer | 67 | 2024-12 | no | no | 6.00 | 8.70 | 8.02 |
| 12CI | Dopamine aptamer DGR-1A (dimer) | aptamer | 82 | 2026-03 | no | no | 16.82 | 10.05 | 10.05 |
| 9LKE | Guanine-II + hypoxanthine | riboswitch | 70 | 2025-01 | no | no | 3.99 | 6.06 | 5.84 |
| 9LKU | 2'-dG-III + Guanosine | riboswitch | 65 | 2025-01 | no | no | 2.06 | 9.02 | 9.02 |
| 9HRF | Class V GTP UU variant (monomer) | aptamer | 70 | 2024-12 | no | no | 7.72 | 16.12 | 15.39 |
| 9MFH | env2 cobalamin riboswitch (apo) | riboswitch | 76 | 2024-12 | no | no | 6.76 | 7.19 | 4.82 |
| 9DE7 | HIV-1 TAR full-length | viral RNA | 57 | 2024-08 | no | no | 3.16 | 3.21 | 3.21 |
| 9J4N | E. coli Leucine tRNA | tRNA | 81 | 2024-08 | no | no | 3.58 | 4.51 | 4.20 |
| 9E9O | SARS-CoV-2 SL5 | viral RNA | 101 | 2024-11 | no | no | 4.42 | 17.37 | 16.63 |

## AF3 per-residue pLDDT calibration vs crystal

| ID | Spearman ρ | p |
|---|---|---|
| 7QR3 | +0.083 | 0.497 |
| 9LJN | -0.276 | 0.0199 |
| 9UW0 | -0.198 | 0.12 |
| 9HRD | -0.581 | 2.6e-07 |
| 12CI | -0.546 | 1.09e-07 |
| 9LKE | +0.089 | 0.463 |
| 9LKU | -0.411 | 0.000832 |
| 9HRF | -0.435 | 0.000168 |
| 9MFH | -0.578 | 4.48e-08 |
| 9DE7 | -0.321 | 0.0148 |
| 9J4N | -0.270 | 0.0148 |
| 9E9O | -0.248 | 0.0124 |

## Aggregate (N = 12 targets with full data)

- DRfold2 mean RMSD: **5.37 Å** (median 4.02)
- AF3 default-seed mean: **8.78 Å** (median 7.95)
- AF3 best-of-5 mean: **8.16 Å** (median 7.06)
- Wilcoxon signed-rank (DRfold2 vs AF3 default): W=8.0, p=0.0122 (N=12)
