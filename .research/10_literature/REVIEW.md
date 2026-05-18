# 10 — Literature Review

> Owner: rna-research. Synthesized view of the literature relevant to the focused
> finding paper. Each entry should cite a paper (full title + venue + year + link)
> and summarize its relevance in 2–4 sentences.

## Sprint 0 Week 1 focus-candidate briefs (companion files)

Companion one-page briefs in `focus_candidates/`, produced 2026-05-18 by
cross-referencing the 11-category scheme (Q-2026-05-18-002 ANSWERED) with the
empirical category counts (run `01_description_categorization`) and the Week-1
reading plan below. Verdicts apply the publishability filter from `00_thesis.md`:

- [`riboswitch.md`](focus_candidates/riboswitch.md) — **STRONG** — biologically central, conformational duality angle, val+test=4.
- [`ribozyme.md`](focus_candidates/ribozyme.md) — **STRONG** — active-site geometry is a sharp falsifiable lens, val+test=4.
- [`pseudoknot.md`](focus_candidates/pseudoknot.md) — **MEDIUM** — strong biology but zero val/test rows; best deployed as a secondary lens on riboswitch/ribozyme.
- [`viral_rna.md`](focus_candidates/viral_rna.md) — **MEDIUM** — heterogeneous super-class; needs sub-categorization to make a sharp claim.
- [`tRNA.md`](focus_candidates/tRNA.md) — **MEDIUM** — canonical tRNA over-studied; story conditional on non-canonical sub-population being present in val+test.
- [`nmr_solution_motif.md`](focus_candidates/nmr_solution_motif.md) — **WEAK** standalone — methodology bucket, zero val/test; best as a methodological appendix (ensemble-aware evaluation).
- [`loop_motif.md`](focus_candidates/loop_motif.md) — **WEAK** standalone — heterogeneous motif bucket; best as a cross-cutting motif-accuracy lens on a biological focus family.

Skipped: `ribosome_subunit` (over-studied + zero val/test) and `complex_with_protein` (out-of-scope per `00_thesis.md`).

Handoff: rna-scientist will use these briefs plus run-01 counts to draft `20_hypotheses/FOCUS_SHORTLIST.md` (see Q-2026-05-18-006).

## Week 1 reading plan (2026-05-18)

This is a **reading plan, not a synthesis**. Deep summaries are produced
file-by-file in `papers/` after the scientist confirms priorities
(see Q-2026-05-18-001 in `questions_for_scientist.md`). Papers ordered by
reading priority for Sprint 0 focus-candidate selection.

### Bucket A — RNA 3D structure prediction (SOTA, what we are looking *with*)

1. **"RibonanzaNet: deep learning of RNA secondary structure and 3D from
   chemical mapping"** — Das lab et al., 2024 (bioRxiv / Nature Methods accepted).
   <https://www.biorxiv.org/content/10.1101/2024.02.24.581671>
   *Rationale:* Directly trained on Ribonanza/Stanford data; the most relevant
   model for the Kaggle benchmark. Must understand its inductive biases before
   discussing where it struggles.

2. **"Accurate prediction of RNA 3D structure with AlphaFold3"** — Abramson
   et al., *Nature*, 2024. <https://www.nature.com/articles/s41586-024-07487-w>
   *Rationale:* The RNA-specific sections of the AF3 paper define the current
   benchmark ceiling for non-coding RNA structure prediction; we must cite this
   carefully and humbly in framing.

3. **"RhoFold+: Accurate RNA 3D structure prediction using a language-model-based
   deep learning approach"** — Shen et al., *Nature Methods*, 2024.
   <https://www.nature.com/articles/s41592-024-02487-0>
   *Rationale:* Second key SOTA for pure-RNA prediction; expected to be one of
   the three models used for inter-model error-variance signals on Sprint 0.

4. **"Boltz-1: Democratizing Biomolecular Interaction Modeling"** — Wohlwend
   et al., MIT, 2024 (bioRxiv). <https://www.biorxiv.org/content/10.1101/2024.11.19.624167>
   *Rationale:* Open-weights AF3-class model that includes RNA. Likely the most
   reproducible "AF3-like" comparator we can actually run on the pod.

5. **"Vfold-Pipeline: a web server for RNA 3D structure prediction from
   sequences"** — Li & Chen, *Bioinformatics*, 2022 (+ subsequent Vfold3D /
   IsRNA work). <https://academic.oup.com/bioinformatics/article/38/16/4042/6633920>
   *Rationale:* Physics-informed comparator; complements the deep-learning
   models and is known to behave differently on pseudoknots and tertiary
   contacts.

### Bucket B — Calibration & uncertainty (the *tool* layer)

6. **"On Calibration of Modern Neural Networks"** — Guo, Pleiss, Sun,
   Weinberger, *ICML*, 2017. <https://arxiv.org/abs/1706.04599>
   *Rationale:* Canonical reference for ECE / reliability diagrams /
   temperature scaling. Foundation for the calibration lens we will apply to
   per-residue confidence on the focus family.

7. **"A Tutorial on Conformal Prediction"** — Angelopoulos & Bates, 2021
   (extended tutorial / *Foundations & Trends in ML*, 2023).
   <https://arxiv.org/abs/2107.07511>
   *Rationale:* Distribution-free uncertainty quantification — natural fit for
   small test sets (62 sequences) where parametric calibration is fragile.

8. **"Confidence-Aware AlphaFold2: pLDDT calibration and reliability of
   structural predictions"** — representative follow-up to AF2 pLDDT
   (e.g., Roney & Ovchinnikov 2022, or Akdel et al. 2022 structurome survey).
   <https://www.nature.com/articles/s41594-022-00849-w> (Akdel et al.)
   *Rationale:* Closest prior art for "confidence-score reliability on a
   structure predictor"; defines methodological vocabulary we should reuse
   rather than reinvent.

### Bucket C — Under-studied RNA classes / hard motifs (the *focus* layer)

9. **"RNA pseudoknots: prediction and modeling"** — Reuter & Mathews-style
   review, e.g. Legendre et al., *Bioinformatics*, 2018 ("DotKnot") or the
   more recent Sloma & Mathews review.
   <https://academic.oup.com/bioinformatics/article/34/16/2790/4953367>
   *Rationale:* Pseudoknots are the most widely-acknowledged blind spot of
   MSA-based predictors and are a strong focus-candidate; need to know what
   has *already* been said before claiming novelty.

10. **"Riboswitches: structure, function, and prediction challenges"** —
    Breaker lab review, e.g. McCown et al., *RNA*, 2017 ("Riboswitch diversity
    and distribution").
    <https://rnajournal.cshlp.org/content/23/7/995>
    *Rationale:* Riboswitches have multiple alternative folds (aptamer ON/OFF)
    — a single ground-truth structure under-represents their biology, which is
    exactly the kind of "structural challenge" angle the paper aims for.
    Strong focus-candidate.

### Stretch (read only if Week 1 closes early)

- Multi-way junctions: Bindewald et al., *NAR*, 2008 ("RNAJunction database").
- Ribozyme catalytic-core geometry: Lilley, *Phil. Trans. R. Soc. B*, 2011.

## How this plan will be executed

- Each paper, once green-lit by the scientist, gets its own file in
  `papers/` named `<lastname>_<year>_<shortslug>.md` with: full citation,
  4-bullet TL;DR, methods, datasets used, what they claim about RNA classes
  relevant to us, gaps / open questions.
- After the first 3 deep summaries land, this REVIEW.md is updated with a
  synthesis section and `gaps.md` / `novel_angles.md` are seeded.
