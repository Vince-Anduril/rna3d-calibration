# Focused Analysis Plan — CPEB3 Human vs Chimpanzee

**Date:** 2026-05-18
**Status:** Strategic pivot — replaces the "train a custom 5M-param predictor" angle.
**Reason:** Custom training spent 5+ pod hours on infrastructure bugs with zero scientific output. CPEB3 gives us a concrete, narrow, high-novelty target where ~30 min of GPU time can produce paper-grade results.

---

## The thesis (one sentence)

*"Do state-of-the-art RNA structure predictors discriminate between functionally divergent single-nucleotide variants? A focused study on the CPEB3 human/chimpanzee ribozyme pair from the Stanford 3D Folding (CASP15) benchmark."*

## Why CPEB3, why now

1. **The Stanford benchmark val/test set has ONLY ribozyme entries from CPEB3** — R1107 (human) and R1108 (chimpanzee), 69 nt each. Anyone studying ribozyme prediction on Stanford is implicitly studying CPEB3.
2. **The two sequences differ at exactly ONE position (30: A→G)** — yet the chimpanzee variant cleaves **~4× faster** (Skilandat et al., RNA 2016).
3. The activity gap is explained by a **P1/P1.1 mispairing** in the human variant — a subtle tertiary distinction that a calibrated model should "see".
4. **Recent benchmarks show large divergence on this exact target**: DRfold2 ≈ 2.72 Å, DeepFoldRNA ≈ 5.68 Å, AlphaFold3 ≈ 7.98 Å (RMSD on the chimp PDB 7QR3). Models disagree by a factor of 3.
5. CPEB3 has **functional importance** (mammalian object-location memory, PMC 10274809, Feb 2024) — biology readers care.
6. CASP15 looked at CPEB3 but **did not analyze the human/chimp pair as a confidence-discrimination probe** — our angle.

## What we DO NOT claim

- We do not claim to beat any SOTA model on accuracy.
- We do not claim to fix the CPEB3 prediction problem.
- We do not claim novelty of the underlying CPEB3 biology (all published).

## What we DO claim

- A focused calibration analysis on a pair of biologically interesting variants that have not been jointly studied as a model-discrimination probe.
- Per-residue confidence and predicted-structure comparison between two near-identical sequences, with attention to the P1/P1.1 region (around position 30) where the activity difference resides.
- Plain-English implications for what "confidence" means in current RNA predictors.

## Concrete artifacts the paper needs

1. **Figure 1** — sequence alignment of R1107/R1108 highlighting position 30, plus a cartoon of the secondary structure with P1/P1.1 marked.
2. **Figure 2** — per-residue confidence (from each predictor) for both variants, overlaid. Highlight what happens around position 30.
3. **Figure 3** — predicted structures (cartoon, superposed via Kabsch) for both variants from one or two predictors. Show whether the predictor places nucleotide 30 in the same local geometry for both variants.
4. **Table 1** — per-predictor: RMSD vs 7QR3 (chimp), per-residue mean error, ECE, mean confidence in the P1/P1.1 window.
5. **Optional Table 2** — if we can resurrect a CASP15 leaderboard for R1107/R1108, contextualize our numbers.

## Models to try (in order of priority)

| Model | Why | Install path | Risk |
|---|---|---|---|
| **DRfold2** | Best CPEB3 perf (RMSD 2.72 Å, PLOS Biology 2025) | PyPI or GitHub | Medium — newish |
| **RhoFold+** | Open, well-known, Nature Methods 2024 | github.com/RFOLD/RhoFold or biomap-research | Medium |
| **trRosettaRNA** | Solid baseline | github.com/zhanggroup/trRosettaRNA | Medium |
| **AlphaFold Server** | Free tier (~20 predictions/day), gives a SOTA reference | Browser/API | Low (works) but slow |
| **vfold** | We have the notebook | local notebook | Already done |

**Plan:** install whichever 2 of the top 3 install cleanly within 30 min on the pod; fall back to AF3-Server + vfold if both fail.

## Compute budget

| Step | Where | Time |
|---|---|---|
| Install 1-2 models | pod | 15–30 min |
| Inference on 2 sequences × N models | pod, GPU | < 5 min per model |
| Per-residue analysis + figures | pod, CPU | < 5 min |
| **Total** | | **30-60 min of pod time** |

This is **6× cheaper** than the failed 3h training plan and produces actually-useful artifacts.

## Files in this pivot (in `code/jobs_gpu/`)

- `focused_cpeb3.py` — main pipeline: define sequences, run each available predictor, save per-residue confidence + xyz, compare.
- `focused_cpeb3_install.sh` — installation steps for DRfold2 / RhoFold+ on the pod.
- `focused_cpeb3_analyze.py` — analyze the saved predictions: per-residue accuracy vs 7QR3, ECE, P1/P1.1 region zoom.

## Out of scope

- Training any model from scratch (already failed at ~5 pod hours; not the right tool here).
- The 5M-param custom transformer in `code/jobs_gpu/model.py` — kept in the repo for reproducibility/audit, not in the publication.
- The Kabsch SVD bug saga — superseded by the focused-inference approach which does not need Kabsch (we have one ground-truth crystal and we can align via PyMOL or biotite).

## Decision log entry to add to `.research/00_thesis.md`

> **2026-05-18 (evening)** — Strategic pivot: drop the custom-predictor training thread, replace with a focused inference comparison on the **CPEB3 ribozyme** human/chimpanzee variants (R1107, R1108) from the Stanford val/test set. Driven by the discovery that **all 4 ribozyme entries in val+test are CPEB3 variants**, and that the **single-nucleotide difference at position 30 produces ~4× activity divergence** in vitro. The paper question becomes: "do SOTA predictors discriminate these variants?". See `docs/superpowers/specs/FOCUSED_CPEB3_PLAN.md`.
