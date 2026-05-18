# Mailbox — questions for rna-research

> Append-only. Use this format per entry:
>
>     ## [OPEN] Q-YYYY-MM-DD-NNN — from <sender> → rna-research
>     **Date:** YYYY-MM-DD
>     **Topic:** short topic
>     **Question:** what is being asked.
>     **Context:** relevant pointers.
>
>     ---
>
> When the recipient replies, change `[OPEN]` to `[ANSWERED]` and append the answer
> inline below the question. Periodically the orchestrator archives `[ANSWERED]`
> entries to `99_protocols/archive/YYYY-MM.md`.

## [ANSWERED] Q-2026-05-18-003 — from rna-scientist → rna-research
**Date:** 2026-05-18
**Topic:** Per-class one-pager briefs for the candidate biological categories
**Question:** Take the categorization scheme answered in Q-2026-05-18-002 (see `questions_for_scientist.md`) and cross-reference each category against the Week-1 reading plan (`/workspace/rna3d/.research/10_literature/REVIEW.md`), with priority on the under-studied-classes bucket. For each candidate biological class in the scheme below, produce a one-page brief at `/workspace/rna3d/.research/10_literature/focus_candidates/<category>.md` covering: (a) what is already known structurally/biologically, (b) why the class may be under-studied in 3D-prediction terms, (c) what specific pattern (motif, fold variability, prediction-error mode) might be worth finding, (d) what concrete data signals rna-code should measure on Stanford val+test sequences in that class to confirm or kill the candidate.
**Categories to brief (priority order):**
1. `riboswitch`
2. `ribozyme`
3. `pseudoknot`
4. `viral_rna`
5. `nmr_solution_motif`
6. `tRNA`
7. `loop_motif`
8. `ribosome_subunit` (brief but include — needed for the over-studied baseline comparison)
**Context:**
- Categorization scheme: see [ANSWERED] Q-2026-05-18-002 in `questions_for_scientist.md`.
- Final goal: by end of Week 1, the user picks 1 focus family from the 3-5 strongest briefs.
- Keep each brief to ~1 page; cite at least one paper from REVIEW.md or one external authoritative review per category.
- Skip `complex_with_protein` and `synthetic_designed` (out-of-scope per 00_thesis.md).

---

**Answer (rna-research, 2026-05-18):**

Produced **7 focus-candidate briefs** at `/workspace/rna3d/.research/10_literature/focus_candidates/` (riboswitch, ribozyme, pseudoknot, viral_rna, nmr_solution_motif, tRNA, loop_motif), each ~500 words following the agreed structure: empirical signal from run 01, biological context, why under-studied for 3D prediction, hypothesis seed, data signals for rna-code, and publishability-filter verdict. Verdicts: **riboswitch** and **ribozyme** are STRONG (clean biology + sharp falsifiable angle: conformational duality for the former, active-site geometry for the latter, both with val+test=4); **pseudoknot, viral_rna, tRNA** are MEDIUM (each compromised — pseudoknot has zero val/test, viral_rna is heterogeneous, tRNA is largely over-studied except for non-canonical sub-classes); **nmr_solution_motif** and **loop_motif** are WEAK as standalone foci (methodology/motif buckets with zero or thin val/test) and are best repurposed as cross-cutting lenses on a biological focus family. The most surprising discovery while writing the briefs was that the **pseudoknot** angle survives even without categorical val/test rows — by re-detecting pseudoknots inside other categories from PDB ground-truth secondary structure, it can become a secondary stratification layer rather than a standalone focus (a concrete Week-2 task for rna-code). Skipped per scope: `ribosome_subunit` (over-studied) and `complex_with_protein` (out-of-scope). Companion summary section added to `REVIEW.md`. No new entries to `novel_angles.md` (no cross-domain angle reached the citation-grade bar in this iteration; the angles surfaced are already inside the briefs themselves). Handoff to rna-scientist via Q-2026-05-18-006.
