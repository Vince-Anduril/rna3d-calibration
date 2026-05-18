# 40 — Consolidated findings (paper-ready insights)

> Cross-agent. Findings are insights that have passed sanity checks and are
> candidates for the paper. Each finding should be defensible in <30 seconds
> (one figure, one claim, one differentiating reference).

---

## F-001 — DRfold2 places the predicted structural consequence of the CPEB3 single-nucleotide variant at the same P1/P1.1 region implicated by published biochemistry

**Date promoted:** 2026-05-18
**Last validated:** 2026-05-19 (Experiment #1 geometric distance check)
**Status:** **DOWNGRADED to a clean model-behavior observation, not a biological discovery.** After running Experiment #1, the F-001 claim reduces cleanly to: *DRfold2 correctly places the published P1.1↔P1 3D contact in its predicted structure and propagates a single-nt perturbation across that contact in a mutation-specific way.* Of the original top-5 divergent residues {9, 22, 24, 51, 60}:
- Residues **9 and 60** are simultaneously (a) within 15 Å of position 30 in the predicted structure (direct geometric coupling, consistent with the published P1.1↔P1 contact) AND (b) specific to the real mutation (0/5 NULL controls have them in their top-5). These are the cleanly defensible cascade residues.
- Residues **22, 24, 51** are at long distance from position 30 (>20 Å) AND each appears in 2/5 NULL control top-5 lists. These are DRfold2 *attractor positions* (residues the model tends to vary regardless of input), not mutation-specific responses.

So the refined finding is parsimoniously explained by "the model correctly learned the documented 3D contact and propagates locally through it" — a sanity check on model fidelity, not evidence of any deeper mechanistic understanding.
**Source runs:** initial inference `.research/30_experiments/runs/cpeb3_focused/SUMMARY.md`; NULL CONTROL was the same SUMMARY file after `scripts/run_cpeb3_full.sh` ran with the pos-41 sequence.

**One-sentence claim:**

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

**Sanity checks:**

1. **Null control. — PASSED 2026-05-18.** Re-ran DRfold2 on R1108 with a
   peripheral mutation A→C at position 41 (J3-P4 region). The top-5 divergent
   residues moved to {51, 52, 48, 50, 49} — all clustered AROUND the new
   mutation site, with only ONE residue (pos 51) overlapping the REAL case.
   The position-30 → P1 cascade is therefore NOT a model prior; it is a
   specific response to the position-30 mutation that happens to put the
   structural consequence at the P1 helix endpoints.

2. **Second predictor (RhoFold+, single-sequence mode). — CONDUCTED 2026-05-18,
   produced a DIFFERENTIAL finding rather than a simple reproduction.**
   RhoFold+ run without MSA (single-seq mode) gave:
   - Much worse RMSD vs ground truth on chimp 7QR3: **7.699 Å (RhoFold+)** vs
     3.726 Å (DRfold2). On the human variant: 13.356 Å vs 4.204 Å.
   - REAL top-5 divergent residues: {48, 49, 50, 21, 51}
   - NULL top-5 divergent residues: {48, 49, 50, 21, 51}
   - Overlap = **5/5** (the same residues dominate the variance regardless of
     which position was mutated).
   RhoFold+ in single-seq mode therefore has a STRONG structural prior that
   masks any mutation-specific response. This is the OPPOSITE pattern from
   DRfold2 (1/5 overlap). The two predictors disagree on the structural
   response to mutation. This is a richer finding than a clean reproduction:
   it shows that whether a predictor reproduces the Skilandat 2016 P1/P1.1
   anchor depends on the predictor itself (and likely on whether MSA is
   provided). DRfold2's RCLM (language-model-based, no MSA needed by design)
   does; RhoFold+ without MSA does not. Re-running RhoFold+ with a proper MSA
   would clarify whether the prior collapses or not — that is the natural
   next experiment.

3. **Per-residue confidence.** Still pending. .ret files on the pod likely
   contain per-residue scores; need to inspect their format.

**Strengthened claim, post-NULL CONTROL:**

DRfold2 predicts a **specifically allosteric-like response** to the human/chimp
CPEB3 mutation — a single nucleotide at position 30 (the P1.1 region) causes
the largest predicted structural change ~30 residues away, at the P1 helix
endpoints (residues 9 and 60). When the same predictor receives a peripheral
mutation (position 41 in J3-P4), the structural response is local and
clustered around the mutation site (residues 48-52). The long-range pattern
seen for the biologically functional mutation is therefore a *response to that
mutation*, not a fixed structural prior of the model. This pattern coincides
with the biochemical mechanism proposed by Skilandat 2016 (RNA) — the
P1.1-mediated mispairing in human CPEB3 that explains its 4× slower cleavage.

**Verdict:** **PAPER-GRADE finding on the null-control axis.** Can be
written into a methods + results section as the headline observation, with
the remaining two sanity checks listed as future work (or done in the same
study if time allows).

**Numbers (DRfold2 NULL CONTROL, 2026-05-18, no extra pod GPU time —
inference of null_pos41 was already complete on the persistent /workspace
from the previous chained run):**

- NULL Kabsch RMSD (null_pos41 vs R1108 predictions): **1.909 Å** (vs 1.858 Å
  for REAL — similar magnitude)
- NULL top-5 divergent residues: **51, 52, 48, 50, 49** (local cluster around
  the mutation site at pos 41)
- Overlap with REAL top-5 {9, 60, 51, 22, 24}: **{51}** — exactly one residue

**Numbers (RhoFold+ single-seq, 2026-05-18, ~3 sec per inference):**

- RhoFold+ R1108 vs 7QR3 chain C: 7.699 Å (DRfold2 was 3.726 Å)
- RhoFold+ R1107 vs 7QR3 chain C: 13.356 Å (DRfold2 was 4.204 Å)
- REAL top-5: {48, 49, 50, 21, 51}
- NULL top-5: {48, 49, 50, 21, 51}
- Overlap REAL/NULL = **5/5** — RhoFold+ in single-seq mode shows a strong
  fixed structural prior; no mutation-specific response detectable.

**Cross-predictor synthesis (paper-grade observation):**

The two predictors give qualitatively different answers to the same question.
DRfold2 (RCLM, no MSA needed) reproduces the Skilandat 2016 biochemical
anchor. RhoFold+ (Evoformer-style, designed for MSA input) without MSA does
not — its predicted variance is dominated by a small set of residues that do
not respond to the input mutation. This is itself a calibration finding:
single-sequence RNA prediction confidence is not uniformly meaningful across
SOTA predictors, and the *kind* of response a predictor produces to a small
input perturbation may be more informative than its raw RMSD ranking.
