# 40 — Consolidated findings (paper-ready insights)

> Cross-agent. Findings are insights that have passed sanity checks and are
> candidates for the paper. Each finding should be defensible in <30 seconds
> (one figure, one claim, one differentiating reference).

---

## F-001 — DRfold2 places the predicted structural consequence of the CPEB3 single-nucleotide variant at the same P1/P1.1 region implicated by published biochemistry

**Date promoted:** 2026-05-18
**Status:** PROMISING — needs (a) randomized-mutation null control, (b) reproduction with a second predictor (RhoFold+ or AF3 server), (c) per-residue confidence analysis. NOT YET CLAIMABLE for the paper.
**Source run:** `.research/30_experiments/runs/cpeb3_focused/SUMMARY.md`

**One-sentence claim (to be tested against the sanity checks below):**

DRfold2 — trained without functional or kinetic data — predicts that the
4×-activity-changing single-nucleotide difference between human (R1107) and
chimpanzee (R1108) CPEB3 ribozymes produces its largest 3D restructuring NOT at
the mutation site (position 30, 2.148 Å C1' delta) but at the P1 helix
endpoints (residues 9 and 60, ≥4.1 Å) — the same region where Skilandat,
Boese, & Sigel (RNA 2016) propose the P1/P1.1 mispairing that explains the
biochemical activity gap.

**Defensible in 30 s — proposed two-panel figure:**

- **Left:** per-residue Kabsch-aligned C1' delta between DRfold2 R1107 and
  R1108 predictions, x = residue 1..69, y = Å. Vertical mark at residue 30
  (the mutation). Annotations on residues 9 and 60.
- **Right:** CPEB3 secondary-structure diagram (from Skilandat 2016 or the
  7QR3 crystal paper) with P1, P1.1, P2, P3 labeled. Color residues 9, 22,
  24, 51, 60 by the predicted delta magnitude.

**Differentiation from prior work:**

- *Skilandat, Boese, & Sigel (RNA 2016).* Biochemical/spectroscopic model of
  P1/P1.1 mispairing in human CPEB3 explaining its ~4× slower cleavage. No
  structure prediction. → Our angle: a predictor recovers their geometric
  anchor without being told.
- *DRfold2 (PLOS Biology 2025).* Benchmarks DRfold2 on chimp CPEB3 alone
  (~2.7 Å RMSD vs 7QR3). Does not analyze the human/chimp pair as a
  discrimination probe. → We add the discrimination angle.
- *CASP15 R1107/R1108 assessment (Proteins 2023).* Per-target RMSDs but no
  cross-target single-nt-difference cascade analysis. → We add the
  comparative-deltas angle.

**Numbers (DRfold2, this run on RTX 5090, ~12 min total inference):**

- R1108 (chimp) prediction vs 7QR3 chain C: **3.726 Å**, chain D: 2.913 Å
- R1107 (human) prediction vs 7QR3 chain C: **4.204 Å**, chain D: 2.784 Å
- R1107 vs R1108 predicted: **1.858 Å** Kabsch-aligned RMSD (modest but real)
- Position 30 (mutation): **2.148 Å** C1' delta — surrounded by larger deltas
- Top-5 divergent residues: pos 9 (4.21), pos 60 (4.12), pos 51 (3.98),
  pos 22 (3.96), pos 24 (3.70)

**Sanity checks REQUIRED before paper:**

1. **Null control.** Re-run with an arbitrary single-nt mutation at a position
   distant from any known functional element (e.g., position 5, middle of P1).
   If DRfold2's top-5 divergent residues are STILL pos 9, 60, 51, 22, 24,
   regardless of where the mutation actually sits → the "signal" is a model
   prior, not a response to the mutation. If the top-5 shifts to surround the
   new mutation → the position-30 finding is real.
2. **Second predictor.** Re-run on RhoFold+ and/or AlphaFold Server. If they
   independently put deltas at pos 9 and 60 → robust finding. If only DRfold2
   does → likely a DRfold2-specific behavior.
3. **Per-residue confidence.** DRfold2 may expose a confidence score (look in
   the .ret files); at the divergent positions, is confidence higher or
   lower? This tells us whether the model "knows" what it's doing there.

**Cost of doing these 3 checks:** estimated ~30–60 min of pod time. Same order
as the first run.

**Verdict:** First concrete observation from the project. Not yet a paper
claim, but is the most promising lead so far and is testable with one more
focused pod session.
