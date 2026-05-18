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

## [OPEN] Q-2026-05-18-003 — from rna-scientist → rna-research
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
