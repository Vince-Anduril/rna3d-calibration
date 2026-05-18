# Mailbox — questions for rna-code

> Append-only. Use this format per entry:
>
>     ## [OPEN] Q-YYYY-MM-DD-NNN — from <sender> → rna-code
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

## [OPEN] Q-2026-05-18-004 — from rna-scientist → rna-code
**Date:** 2026-05-18
**Topic:** Implement and run the description-based biological categorization scheme on all Stanford splits
**Question:** Implement the categorization scheme from [ANSWERED] Q-2026-05-18-002 (`questions_for_scientist.md`) as a Python script on the pod, apply it to `train_sequences.csv`, `validation_sequences.csv`, and `test_sequences.csv`, and report which categories pass the ≥10-sequence publishability viability threshold (counted across train+val+test combined).
**Deliverables under `/workspace/rna3d/.research/30_experiments/runs/01_description_categorization/`:**
- `script.py` — applies the 11 priority-ordered regex rules (ribosome_subunit, riboswitch, ribozyme, tRNA, viral_rna, pseudoknot, nmr_solution_motif, loop_motif, complex_with_protein, synthetic_designed, other) with case-insensitive matching on `description`; for `nmr_solution_motif` also gate on `len(sequence) < 60`.
- `config.yaml` — input CSV paths, the regex table (one entry per category), priority order, and the length threshold.
- `results/category_assignment.csv` — columns: `split, target_id, description, length, category, matched_rule, confidence` (confidence = "high" if exactly one rule matched before priority resolution, "medium" if multiple matched, "low" if fallback to `other`).
- `results/category_distribution.png` — grouped bar chart, x=category, y=count, hue=split (train/val/test).
- `results/uncategorized.csv` — all rows where category=`other`, for manual review.
- A short `results/SUMMARY.md` listing per-category counts per split, and explicitly flagging which categories have ≥10 sequences across train+val+test combined (the viability threshold).
**Context:**
- Scheme source: [ANSWERED] Q-2026-05-18-002 in `questions_for_scientist.md`.
- Data paths: see `/workspace/rna3d/data/MANIFEST.md`.
- Use the same pod-only execution convention as run `00_pod_smoke`.
- Report back via a new Q-005 in `questions_for_scientist.md` once results are committed.

---
