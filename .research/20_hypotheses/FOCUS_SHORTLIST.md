# Sprint 0 Week 1 — Focus-candidate shortlist (ranked)

> Owner: rna-scientist.
> Created: 2026-05-18.
> Inputs synthesized: run `01_description_categorization` data signals +
> `10_literature/focus_candidates/*.md` literature briefs.
> Output: user-facing ranked list. The user picks ONE focus family at end of
> Week 1; the picked family fixes H-001's variable F.

## TL;DR

Three families pass the publishability filter at MEDIUM/STRONG strength on the
current Stanford categorization: **ribozyme** (#1, STRONG, val+test=4),
**riboswitch** (#2, STRONG, val+test=4), and **viral_rna** (#3, MEDIUM,
val+test=2). All three combine biological depth, a sharp falsifiable angle,
and 3-month feasibility; ribozyme edges ahead on active-site geometry as a
load-bearing structural feature. Regardless of pick, pseudoknot re-detection
across the chosen family (from PDB ground-truth secondary structure) is
queued as a Week-2 cross-cutting lens — pseudoknot itself is too data-thin
to stand alone but cuts cleanly across all three top candidates.

## Ranking criteria applied

The four criteria from `00_thesis.md` are: (1) biological meaning, (2) under-
study for 3D prediction, (3) 3-month feasibility without wet-lab, (4) val/test
sufficiency. Tie-breakers, in order: (a) sharpness of the falsifiable claim at
small n (favors families with a load-bearing per-residue substructure); (b)
val+test count (a hard floor — zero kills the novel-prediction story on the
Stanford benchmark); (c) sub-classification overhead (families requiring less
extra annotation work rank higher in a 3-month budget).

## Ranked candidates

### 1. ribozyme — STRONG

- Empirical: train=54, val=2, test=2, total=58 (val+test=4).
- Why #1: catalysis is the defining function and active-site geometry is the
  load-bearing structural feature — a sub-RMSD restricted to a few catalytic
  residues per class is the sharpest small-n metric in the shortlist. The
  family also overlaps naturally with the pseudoknot lens (HDV, twister) and
  with MSA-depth stratification across classical small ribozymes.
- Sharp falsifiable angle: SOTA predictors are systematically over-confident
  at the catalytic core — predicted pLDDT-equivalent is high while local
  active-site RMSD is high — on at least one of the four val+test ribozymes.
- Risks / caveats: requires per-class active-site residue annotation
  (~one extra week of literature work); val+test=4 supports case-study
  framing, not class-wide statistical claims.

### 2. riboswitch — STRONG

- Empirical: train=17, val=2, test=2, total=21 (val+test=4).
- Why #2: conformational duality (aptamer ON/OFF, ligand-bound vs apo) is a
  genuine epistemological mismatch with single-structure supervised training,
  and small Rfam families for orphan classes give MSA-driven methods limited
  signal. Biological centrality (antibiotic targets, gene regulation) gives
  the paper a clear translational hook.
- Sharp falsifiable angle: predictor confidence is mis-calibrated specifically
  at the ligand-binding pocket — high pLDDT-equivalent, high local RMSD — on
  the 4 val+test riboswitches.
- Risks / caveats: total n=21 is the smallest viable bucket; per-class
  diversity within those 21 rows is unverified (might be dominated by one
  class, which would weaken the cross-class story).

### 3. viral_rna — MEDIUM

- Empirical: train=75, val=1, test=1, total=77 (val+test=2).
- Why #3: largest train-set volume of any candidate; biomedically loaded
  (frameshifting, IRES, packaging); high overlap with the pseudoknot blind
  spot (frameshift elements, tRNA-mimic 3' UTRs). After the v2 regex rescue
  (picornaviral / bacteriophage / coliphage → viral_rna; queued as Q-007),
  the bucket will grow further from the current `other` pile.
- Sharp falsifiable angle: a within-class split — small bacteriophage stems
  (MS2, boxB) easy, complex IRES / frameshift elements hard — yields a
  clean small-vs-complex contrast.
- Risks / caveats: heterogeneous by construction; val+test=2 almost certainly
  spans two different sub-families, which forecloses single-sub-class
  statistics; the story requires sub-categorization work before any sharp
  claim, eating into the 3-month budget.

## Not on the shortlist (reasoned exclusions)

- `ribosome_subunit` (n=199 train, val+test=0): over-studied + zero val/test
  rows kill the novel-prediction angle on the Stanford benchmark.
- `complex_with_protein` (n=43, val+test=4): out-of-scope per `00_thesis.md`
  (pure-RNA only).
- `tRNA` (n=66, val+test=2): canonical tRNA is over-studied at the global
  fold; the non-canonical sub-population (mitochondrial, suppressor, mimic)
  is interesting but requires sub-classification and is conditional on those
  sub-classes actually being present in the 2 val+test rows.
- `pseudoknot` (n=11, val+test=0): strong biology but zero val/test rows
  forecloses the standalone categorical story. Repurposed as the cross-
  cutting secondary lens (see next section).
- `nmr_solution_motif` (n=131, val+test=0): methodology bucket with zero
  val/test rows; best deployed as a methodological appendix (ensemble-aware
  evaluation) on whichever biological family is chosen.
- `loop_motif` (n=52, val+test=2): heterogeneous motif bucket; weak as a
  standalone focus, best repurposed as a per-motif accuracy lens layered on
  the chosen family.
- `synthetic_designed` (n=4): fails the >=10 viability floor.

## Secondary stratification (regardless of pick)

Per rna-research's brief on `pseudoknot.md`, pseudoknot re-detection on PDB
ground-truth secondary structure should be run across **all** Stanford rows
(not just the n=11 categorical bucket). This produces a `has_pseudoknot`
boolean column that cuts cleanly across whichever family the user picks,
enabling a nested-vs-crossed base-pair recall metric and a check on whether
the pseudoknot blind spot inside the chosen family explains residual error.
Queued as a Week-2 deliverable that runs in parallel with focus-specific work
(included in Q-007 to rna-code).

## User decision needed

Pick ONE focus family from the top 3. Suggested decision frame:

- **If you care about the sharpest, most falsifiable claim**, pick
  **ribozyme**: the active-site sub-RMSD is the cleanest small-n metric and
  the failure mode (high-confidence wrong-geometry at catalytic residues) is
  biologically concrete and visually compelling in a figure.
- **If you care about biological depth and conformational reasoning**, pick
  **riboswitch**: the ON/OFF duality + ligand-pocket calibration story is
  conceptually richer and connects naturally to RNA gene-regulation biology.
- **If you want safety margin on data volume** (train-set robustness, room
  to subset and reframe), pick **viral_rna**: largest n by far, the v2 regex
  refinements grow it further, and the frameshift / IRES sub-classes overlap
  with the pseudoknot lens for free.

Please commit your pick by appending a line to this file
(`User pick: <family> — <date>`) or by replying in the orchestrator session;
rna-scientist will then fork H-001 into H-001a (focus-specific) and queue
rna-code's first real inference run for the chosen family.
