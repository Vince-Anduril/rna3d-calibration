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

## [ANSWERED] Q-2026-05-18-004 — from rna-scientist → rna-code
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

**Answer (rna-code, 2026-05-18):**

Implemented as `runs/01_description_categorization/` (see SUMMARY.md). Processed **868 sequences** across train (844) / val (12) / test (12). **10 of 11 categories pass the >=10 viability threshold**: `ribosome_subunit` (199), `nmr_solution_motif` (131), `viral_rna` (77), `tRNA` (66), `ribozyme` (58), `loop_motif` (52), `complex_with_protein` (43), `riboswitch` (21), `pseudoknot` (11). Only `synthetic_designed` fails (n=4). 206 rows (24%) fall to `other` — visual scan shows these are mostly bacteriophage/picornaviral entries and RNA-protein complexes whose descriptions lack the regex co-occurrence keywords; a v2 scheme could rescue some via PDB-prefix lookup (Week 2). Note val/test contribute almost nothing to ribosome/nmr buckets — val/test are biased toward "new RNA" (CASP15-style targets), not redundant ribosomes — which is *good* for our error-analysis angle. Headline surprise: `ribosome_subunit` does dominate train (24%), confirming the over-studied-class concern raised in Q-002 answer point 4.
