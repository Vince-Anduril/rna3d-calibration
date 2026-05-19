# Experiment #2 — n=20 random NULL controls

Generated: 2026-05-19T07:03:09.200916+00:00

## Headline

- **REAL case** (pos 30, biological): both anchors (9, 60) in top-5 = `True`
- **NULL controls** (n=20, random single-nt at uniformly sampled positions):
  - with pos 9 in top-5:  **2/20**
  - with pos 60 in top-5: **3/20**
  - with **BOTH** 9 AND 60 in top-5: **0/20**
  - with EITHER 9 OR 60 in top-5:  **5/20**

## Statistical interpretation

- Add-1-smoothed empirical p-value for ``both anchors in top-5'' under H0 (random prior): **0.0476**
- Theoretical baseline (uniform random top-5 of 5 out of 69): 0.0043

**Outcome:** No random NULL mutation reproduces the ``both anchors'' pattern. The empirical p-value upper bound is 0.048 (1/(n+1) with add-1 smoothing). The biological mutation is statistically distinguishable from the random-mutation null at this n.

## Detailed table

| target | pos | kind | mean Δ | max Δ | Δ@9 | Δ@60 | top-5 | both anchors? |
|---|---|---|---|---|---|---|---|---|
| R1107_human | 30 | REAL | 1.64 | 4.21 | 4.21 | 4.12 | {9,22,24,51,60} | **YES** |
| null_rand_02 | 2 | NULL | 0.97 | 4.18 | 1.13 | 3.23 | {22,23,24,25,60} | no |
| null_rand_03 | 3 | NULL | 1.12 | 3.98 | 0.83 | 2.25 | {1,2,24,51,52} | no |
| null_rand_04 | 4 | NULL | 1.40 | 4.26 | 0.83 | 3.20 | {1,2,3,4,60} | no |
| null_rand_06 | 6 | NULL | 2.46 | 7.40 | 2.18 | 1.22 | {24,25,30,50,51} | no |
| null_rand_07 | 7 | NULL | 2.24 | 6.65 | 4.07 | 3.84 | {9,22,24,30,50} | no |
| null_rand_09 | 9 | NULL | 0.86 | 4.35 | 3.86 | 1.31 | {9,22,24,49,50} | no |
| null_rand_14 | 14 | NULL | 2.30 | 4.87 | 1.62 | 2.12 | {25,49,50,51,61} | no |
| null_rand_15 | 15 | NULL | 1.57 | 5.52 | 1.09 | 1.85 | {24,51,61,62,63} | no |
| null_rand_28 | 28 | NULL | 2.91 | 11.41 | 1.36 | 3.54 | {24,25,26,27,28} | no |
| null_rand_29 | 29 | NULL | 3.57 | 12.29 | 2.82 | 0.82 | {24,25,26,28,29} | no |
| null_rand_33 | 33 | NULL | 1.28 | 3.16 | 0.67 | 2.76 | {22,23,24,47,60} | no |
| null_rand_36 | 36 | NULL | 0.97 | 3.60 | 0.78 | 0.52 | {21,22,25,26,51} | no |
| null_rand_37 | 37 | NULL | 1.24 | 3.70 | 0.68 | 1.65 | {25,38,39,51,59} | no |
| null_rand_39 | 39 | NULL | 2.18 | 9.09 | 1.54 | 2.27 | {21,22,23,24,51} | no |
| null_rand_45 | 45 | NULL | 0.81 | 2.72 | 0.68 | 0.73 | {22,25,26,50,51} | no |
| null_rand_49 | 49 | NULL | 1.14 | 4.60 | 1.14 | 0.88 | {22,24,48,49,51} | no |
| null_rand_57 | 57 | NULL | 9.66 | 28.45 | 2.04 | 7.74 | {47,48,49,50,51} | no |
| null_rand_59 | 59 | NULL | 1.56 | 6.49 | 0.84 | 2.24 | {22,23,24,25,51} | no |
| null_rand_63 | 63 | NULL | 1.32 | 7.60 | 0.82 | 1.87 | {22,50,61,62,63} | no |
| null_rand_69 | 69 | NULL | 0.94 | 3.08 | 0.65 | 0.53 | {22,23,24,48,51} | no |
