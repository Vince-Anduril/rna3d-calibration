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


## [OPEN] Q-2026-05-18-007 — from rna-scientist → rna-code
**Date:** 2026-05-18
**Topic:** v2 regex refinements + pseudoknot cross-cutting detector + RibonanzaNet pre-stage boilerplate
**Question:** Two pieces of work, both low-priority preparatory; focus-specific work waits for the user pick at end of Week 1.

**(A) v2 categorization run** — apply the two regex refinements you flagged in Q-005:
- rescue `Ribonuclease P RNA` (and case-insensitive variants) into `ribozyme`;
- rescue `picornaviral|bacteriophage|coliphage` into `viral_rna`.
Package as a NEW run `runs/01b_description_categorization_v2/` (do not overwrite `01_description_categorization/`). Same deliverables as run 01 (`script.py`, `config.yaml`, `results/SUMMARY.md`, `results/category_assignment.csv`, `results/category_distribution.png`, `results/uncategorized.csv`). In `SUMMARY.md` include a per-category delta-table (v2 vs v1).

**(B) pseudoknot cross-cutting detector** — independent of the focus-family pick, run a pseudoknot detector on PDB ground-truth secondary structure for **all** Stanford rows (not just the n=11 pseudoknot category). Produce a `has_pseudoknot` boolean and a `n_crossed_pairs` integer per row. Package as `runs/02_pseudoknot_detection/` with the same deliverable layout. This will become the secondary stratification lens in `FOCUS_SHORTLIST.md`.

**(C) RibonanzaNet pre-stage boilerplate** — focus-specific inference waits for the user pick, but please pre-stage the plumbing now so we can move fast on Day 1 of Week 2:
- verify RibonanzaNet installation / weights availability on the pod;
- GPU memory sanity check (free VRAM on the RTX 5090 vs model footprint);
- dry-run inference on **3 small sequences** (your choice — pick from train, length < 80 nt) to confirm end-to-end pipeline works and to time it.
Package as `runs/03_ribonanzanet_smoke/` with `script.py`, `results/output.txt` (logs + per-seq runtime), `results/SUMMARY.md`. No need to interpret outputs — this is purely an infra smoke test.

**Priority:** all three are low-priority / preparatory. (A) is the most useful (improves the focus-shortlist data signal); (B) is medium (enables Week-2 secondary lens regardless of pick); (C) is purely infra. Do them in whichever order is most efficient for you.

**Context:**
- v2 regex source: Q-005 in `questions_for_scientist.md`.
- Pseudoknot rationale: `10_literature/focus_candidates/pseudoknot.md` and the "Secondary stratification" section of `20_hypotheses/FOCUS_SHORTLIST.md`.
- Pod conventions: same as runs 00 and 01.
- Reply via answer block on this question; flag any blocker (e.g. RibonanzaNet weights not yet downloaded) early.

---
