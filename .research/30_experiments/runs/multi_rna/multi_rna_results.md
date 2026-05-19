# Cross-family analysis — Experiment G

Generated: 2026-05-19T13:01:39.952857+00:00

All AF3 runs in this table are **single-chain** (one job per RNA).

## Ground-truth RMSD per target

| ID | RNA | family | L | deposit | DRfold2 seen? | AF3 seen? | DRfold2 (Å) | AF3 default (Å) | AF3 best (Å) |
|---|---|---|---|---|---|---|---|---|---|
| 7QR3 | CPEB3 ribozyme (R1108) | HDV ribozyme | 69 | 2021-12 | likely yes (Lee 2025 trained post-2022) | unlikely (cutoff Sep 2021) | 3.73 | 6.54 | 5.35 |
| 9LJN | Guanine-II riboswitch | riboswitch | 71 | 2025-01 | unlikely | no | 4.06 | 6.55 | 6.09 |
| 9UW0 | 2'-dG-III riboswitch | riboswitch | 63 | 2025-05 | unlikely | no | 2.14 | 10.05 | 9.33 |
| 9HRD | Class V GTP aptamer | aptamer | 67 | 2024-12 | unlikely | no | 6.00 | 8.70 | 8.02 |
| 12CI | Dopamine aptamer DGR-1A | aptamer | 82 | 2026-03 | no | no | 16.79 | 10.05 | 10.05 |

## AF3 per-residue pLDDT calibration vs crystal

| ID | Spearman ρ | p |
|---|---|---|
| 7QR3 | +0.083 | 0.497 |
| 9LJN | -0.276 | 0.0199 |
| 9UW0 | -0.198 | 0.12 |
| 9HRD | -0.581 | 2.6e-07 |
| 12CI | -0.546 | 1.09e-07 |

## Aggregate (N = 5 targets with full data)

- DRfold2 mean RMSD: **6.55 Å** (median 4.06)
- AF3 default-seed mean: **8.38 Å** (median 8.70)
- AF3 best-of-5 mean: **7.77 Å** (median 8.02)
- Wilcoxon signed-rank (DRfold2 vs AF3 default): W=4.0, p=0.438 (N=5)
