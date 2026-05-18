# 00 — Thesis (living document)

> **Last updated:** 2026-05-18 (post-migration, thesis-pivot)
>
> Read first, every invocation. Update only after a major decision (focus chosen,
> go/no-go, pivot, scope change).

## Direction

**Focused, humble finding paper** — NOT a SOTA critique.

Pick ONE under-studied RNA family or motif from the Stanford Kaggle test set,
chosen opportunistically from data signals in Sprint 0. Deep-dive its structure,
biology, and prediction-error properties. Find a novel pattern worth the
community's attention.

Working title pattern:

> *"Structural prediction challenges in [chosen RNA family X]: a focused analysis
> of the Stanford RNA 3D benchmark."*

## Hard framing constraints

- The paper offers an **additional viewpoint**, not a takedown of any prior work.
  Avoid words like "failures", "errors of", "limitations of [team]". Prefer
  "structural challenges in", "open questions around", "a closer look at".
- Acknowledge what SOTA does well **before** discussing where pattern X appears.
- Cite post-competition work (AF3, Boltz-1, Chai-1, RoseTTAFold2NA) honestly.

## Focus-selection criteria (Sprint 0 Week 1)

Candidate families/motifs are surfaced by rna-research from cross-referencing:
1. Existing literature on RNA structure prediction challenges (which classes are
   known to be hard? which classes have small representation in MSA-based
   training sets?).
2. Practical data signals on the Stanford test set, measured by rna-code:
   - Inter-model error variance per sequence (where do RibonanzaNet, RhoFold,
     vfold most disagree?).
   - Per-class average error (which Rfam family is hardest?).
   - Motif-frequency × motif-error product (which structural motifs are common
     AND poorly predicted?).
3. rna-scientist applies a "publishability filter": is the candidate
   biologically meaningful, sufficiently unstudied (not yet a hot AF3 follow-up),
   small enough to deep-dive in 3 months?

User makes the final pick at end of Week 1 from a shortlist of 3-5 candidates.

## Calibration in the new framing

Calibration is a **tool**, not the headline. Once the focus family is chosen,
calibration metrics (reliability diagrams, ECE) become one of the lenses used to
characterize prediction behavior on that family — alongside structural metrics,
biological annotations, and prediction-error stratification.

## Out of scope for V1

- Training new models or fine-tuning (analysis-only).
- AlphaFold3 real weights (server API stretch only).
- RNA–protein complexes (pure RNA structures only).
- Dynamics / kinetics (static structures only).
- Multimodal data integration (chemical probing etc).

## Decision log

- **2026-05-17** — Brainstorming: angle C (error analysis), later refined to
  calibration as the primary thesis.
- **2026-05-17** — Pod-only execution rule (refined 2026-05-18 with Mac mirror of
  important artifacts).
- **2026-05-18** — Old RTX 4090 pod became unavailable; lost /workspace/rna3d/
  including the spec (354 lines) and plan (1757 lines). Migrated to a new RTX
  5090 pod, 32 GB VRAM.
- **2026-05-18** — Thesis pivoted: NOT a calibration headline. Focus on ONE
  under-studied RNA family, opportunistically picked from data signals in Week 1.
  Calibration becomes a tool, not the thesis title.
- **2026-05-18 (afternoon)** — **FOCUS PICKED: ribozyme** (autonomous orchestrator decision under user grant of direction authority). Rationale: sharpest falsifiable metric (sub-RMSD on a small catalytic active site, 3-7 residues per class), biological centrality of catalysis, bounded scope (~10 well-characterized classes), val+test presence (4 ribozymes), visually compelling figure-1 candidate (active-site overlay). H-001 instantiated with F = ribozyme. See `docs/superpowers/specs/FOCUS_DECISION.md`. Pipeline ready at `code/jobs_gpu/` (5M-param transformer trained from scratch in 2 stages); awaiting pod restart.

- **2026-05-18 (evening)** — **STRATEGIC PIVOT to CPEB3 focus, drop custom-training thread.** Discovery: all 4 ribozyme entries in Stanford val/test are **CPEB3 ribozyme variants** — R1107 (human) and R1108 (chimpanzee), 69 nt each, sharing identical sequence except at position 30 (A→G). Published literature: the chimp variant cleaves ~4× faster (Skilandat et al., RNA 2016), explained by P1/P1.1 mispairing in human. Recent benchmarks: AlphaFold3 RMSD 7.98 Å vs DRfold2 RMSD 2.72 Å on the chimp crystal (PDB 7QR3) — a 3× gap. Custom-training thread had burned ~5 pod hours on infrastructure bugs (SVD-bfloat16, OOM at large batch, Kabsch-loop bottleneck) with zero scientific output; replaced by a **focused inference pipeline** that runs N existing predictors on the 2 sequences and compares per-residue confidence + accuracy. New paper question: *"Do SOTA RNA predictors discriminate between near-identical sequences that nature makes functionally divergent?"*. Full plan: `docs/superpowers/specs/FOCUSED_CPEB3_PLAN.md`. Code: `code/jobs_gpu/focused_cpeb3{,_install.sh,_analyze.py}`.
