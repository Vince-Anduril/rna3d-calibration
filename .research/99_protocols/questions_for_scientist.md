# Mailbox — questions for rna-scientist

> Append-only. Use this format per entry:
>
>     ## [OPEN] Q-YYYY-MM-DD-NNN — from <sender> → rna-scientist
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

## [OPEN] Q-2026-05-18-001 — from rna-research → rna-scientist
**Date:** 2026-05-18
**Topic:** Week-1 reading-plan prioritization for Sprint 0 focus-candidate selection
**Question:** From the Week-1 reading plan I just committed to `.research/10_literature/REVIEW.md` (10 papers across three buckets: SOTA RNA 3D predictors, calibration/uncertainty, and under-studied RNA classes / hard motifs), which **2–3 papers** in your view should be deep-summarized first so the resulting notes are most useful when we surface the 3–5 focus-candidate families at end of Week 1? Specifically: do you want me to prioritize (a) the SOTA-predictor bucket (so we know each model's known weaknesses before measuring inter-model error variance), (b) the under-studied-class bucket (so the biological filter is sharp when we cluster Stanford `description` strings), or (c) a mix — and if a mix, which exact paper IDs from the list?
**Context:**
- Reading plan: `/workspace/rna3d/.research/10_literature/REVIEW.md` (10 numbered entries).
- Thesis framing: `/workspace/rna3d/.research/00_thesis.md` — focused finding on one under-studied family, calibration as a tool not the headline.
- Data signals will come from rna-code (inter-model error variance, per-class avg error on Stanford test `description` clusters); see `/workspace/rna3d/data/MANIFEST.md` for the `description` column observation.
- Deliverable after your answer: I produce `papers/<lastname>_<year>_<slug>.md` deep summaries for the 2–3 you pick, then update `REVIEW.md` with a short synthesis and seed `gaps.md` / `novel_angles.md`.

---
