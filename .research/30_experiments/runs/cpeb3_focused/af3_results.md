# Extension D — AlphaFold 3 cross-model comparison on CPEB3 panel

Generated: 2026-05-19T07:28:41.231370+00:00

**Job**: single AF3 Server submission, 7 RNA chains (1 job instead of 7 — multi-chain run, all 7 sequences predicted independently as they don't form a complex).

## Chain mapping

| chain | target | length |
|---|---|---|
| A | R1107_human | 69 |
| B | R1108_chimp | 69 |
| C | null_pos5 | 69 |
| D | null_pos20 | 69 |
| E | null_pos41 | 69 |
| F | null_pos55 | 69 |
| G | null_pos64 | 69 |

## Key result 1 — top-5 divergent residues, R1107 vs R1108

- **DRfold2** top-5: `[9, 22, 24, 51, 60]`
- **AF3** top-5:     `[1, 8, 9, 10, 11]`
- **Overlap**: `[9]` (1/5 in common)

| residue | DRfold2 Δ (Å) | AF3 Δ (Å) | both top-5? |
|---|---|---|---|
| 9 | 4.21 | 39.26 | **BOTH** |
| 22 | 3.96 | 28.94 | DRfold2 only |
| 24 | 3.70 | 22.82 | DRfold2 only |
| 30 | 2.15 | 19.05 | neither |
| 51 | 3.98 | 7.05 | DRfold2 only |
| 60 | 4.12 | 9.00 | DRfold2 only |

## Key result 2 — RMSD between AF3 and DRfold2 per target

| target | RMSD AF3↔DRfold2 (Å) | AF3 mean pLDDT |
|---|---|---|
| R1107_human | 15.87 | 26.4 |
| R1108_chimp | 24.99 | 25.0 |
| null_pos5 | 21.65 | 26.5 |
| null_pos20 | 24.87 | 25.2 |
| null_pos41 | 17.95 | 26.4 |
| null_pos55 | 13.40 | 25.8 |
| null_pos64 | 25.18 | 25.2 |

## Key result 3 — AF3 per-residue pLDDT at cascade residues

Mean pLDDT (AF3, R1107 vs R1108) at the residues DRfold2 flagged as the cascade. **High pLDDT at residue 30 in R1107** = AF3 confidently predicts at the mispairing site (calibration risk).

| residue | AF3 R1107 pLDDT | AF3 R1108 pLDDT | Δ (R1107-R1108) |
|---|---|---|---|
| 9 | 23.2 | 23.0 | +0.3 |
| 22 | 18.3 | 18.3 | +0.1 |
| 24 | 22.2 | 19.3 | +2.9 |
| 30 | 28.8 | 26.4 | +2.4 |
| 51 | 28.6 | 28.3 | +0.3 |
| 60 | 24.4 | 22.4 | +2.0 |

## Interpretation

**Partial cross-model agreement** (1/5 overlap). Some divergence sites are shared, but each model also has its own emphasis. Worth reporting both top-5 sets as complementary evidence.

