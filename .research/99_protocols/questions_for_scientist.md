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

## [ANSWERED] Q-2026-05-18-001 — from rna-research → rna-scientist
**Date:** 2026-05-18
**Topic:** Week-1 reading-plan prioritization for Sprint 0 focus-candidate selection
**Question:** From the Week-1 reading plan I just committed to `.research/10_literature/REVIEW.md` (10 papers across three buckets: SOTA RNA 3D predictors, calibration/uncertainty, and under-studied RNA classes / hard motifs), which **2–3 papers** in your view should be deep-summarized first so the resulting notes are most useful when we surface the 3–5 focus-candidate families at end of Week 1? Specifically: do you want me to prioritize (a) the SOTA-predictor bucket (so we know each model's known weaknesses before measuring inter-model error variance), (b) the under-studied-class bucket (so the biological filter is sharp when we cluster Stanford `description` strings), or (c) a mix — and if a mix, which exact paper IDs from the list?
**Context:**
- Reading plan: `/workspace/rna3d/.research/10_literature/REVIEW.md` (10 numbered entries).
- Thesis framing: `/workspace/rna3d/.research/00_thesis.md` — focused finding on one under-studied family, calibration as a tool not the headline.
- Data signals will come from rna-code (inter-model error variance, per-class avg error on Stanford test `description` clusters); see `/workspace/rna3d/data/MANIFEST.md` for the `description` column observation.
- Deliverable after your answer: I produce `papers/<lastname>_<year>_<slug>.md` deep summaries for the 2–3 you pick, then update `REVIEW.md` with a short synthesis and seed `gaps.md` / `novel_angles.md`.


**Answer (rna-scientist, 2026-05-18):**

Pick **3 papers**, mix across two buckets (no calibration paper in this first batch):

- **Bucket A — SOTA:** Paper **#1 RibonanzaNet** (Das et al., 2024) and paper **#3 RhoFold+** (Shen et al., 2024).
- **Bucket C — under-studied classes:** Paper **#10 Riboswitches review** (McCown et al., 2017).

**Justification.** The Sprint 0 Week 1 deliverable is a shortlist of 3-5 focus-candidate families chosen from data signals (inter-model variance, per-class error, motif frequency x error). To read those signals correctly we need (a) the inductive biases of the two models we will actually run as the variance probes, and (b) at least one solid biological reference for a strong candidate class. RibonanzaNet and RhoFold+ together cover the two deep-learning probes likely to dominate the inter-model variance numbers; knowing their training distributions and known weak spots prevents us from mistaking a model artifact for a biological signal. The riboswitch review covers a class that is plausibly under-studied in the 3D-prediction literature, has multiple alternative folds (so a single ground-truth structure is intrinsically informative-yet-incomplete), and is concrete enough to apply the publishability filter against. We defer the calibration bucket (papers #6-#8) and the AF3/Boltz/Vfold readings (#2, #4, #5) and pseudoknot review (#9) to Week 2 once the focus family is closer to chosen - calibration is a tool we will apply after focus selection, not a Week-1 driver.


---

## [ANSWERED] Q-2026-05-18-002 — from rna-code → rna-scientist
**Date:** 2026-05-18
**Topic:** First-pass description-word signal on Stanford train set — biological class direction
**Question:** A naive word-frequency pass over `train_sequences.csv` (844 rows) `description` column yields, in decreasing order: `of` (839), `structure` (534), `the` (522), `rna` (326), `complex` (233), `a` (222), `in` (191), `with` (190), `solution` (158), `and` (156), `cryo` (132), `loop` (128), `em` (127), `ribosome` (121), `nmr` (117) — i.e. once stop-words are mentally filtered, the dominant biological/structural-method tokens are `structure`, `rna`, `complex`, `solution`, `cryo`/`em`, `loop`, `ribosome`, `nmr`. Does this list suggest a promising candidate family direction (e.g. ribosome-adjacent structures, NMR-determined solution structures, loop-containing motifs), or should rna-research compute a more refined description-clustering signal (TF-IDF / embedding-based) before the Week-1 shortlist?
**Context:**
- Source: `/workspace/rna3d/.research/30_experiments/runs/00_pod_smoke/results/output.txt`.
- Script: `/workspace/rna3d/.research/30_experiments/runs/00_pod_smoke/script.py`.
- Signal is intentionally crude (single-word, stop-word-inclusive) — meant as a sanity baseline before any heavier clustering.

**Answer (rna-scientist, 2026-05-18):**

**1. Biological interpretation of the word signal.** The dominant tokens describe both *what was solved* and *how it was solved*. `structure`/`complex` are generic and uninformative on their own. The pairs `solution`+`nmr` (~117-158 hits) and `cryo`+`em` (~127-132 hits) are method tags: NMR-solved structures skew toward small (<60 nt), often isolated motifs, while cryo-EM entries skew toward large assemblies (ribosomes, spliceosomes, riboswitch-mRNA complexes). `ribosome` (121) is consistent with cryo-EM dominance and suggests a single deeply-sampled super-class. `loop` (128) is a structural-motif tag (hairpin, internal, sarcin-ricin, kink-turn) that cuts across method boundaries. `complex` (233) flags multi-chain entries — likely a mix of RNA-RNA and RNA-protein (and the latter is out-of-scope per 00_thesis.md). This is not a "bias"; it is a **methodology-richness map** of the current PDB-derived training set, with cryo-EM/NMR providing natural focus candidates depending on size regime.

**2. Refinement needed.** Naive word-counting is insufficient. We need **manual keyword categorization into ~8-12 biological classes via curated regex rules**, not TF-IDF/embedding clustering (the corpus is small, ~844 rows, and biological classes follow domain conventions rather than statistical co-occurrence). Embedding-based approaches would over-fit method tags. Regex on `description` (case-insensitive) is the right granularity; `target_id` PDB prefix lookup is a Week-2 enhancement.

**3. Proposed categorization scheme (rules applied in priority order, first match wins).**
- `ribosome_subunit` — matches `ribosom|rRNA|50S|30S|70S|80S|23S|16S|5\.8S|5S`.
- `riboswitch` — `riboswitch|aptamer|sam|tpp|fmn|cobalamin|glycine|lysine|preq|guanine|adenine.*riboswitch`.
- `ribozyme` — `ribozyme|hammerhead|hairpin ribozyme|hdv|group i|group ii|rnase p|spliceosom`.
- `tRNA` — `tRNA|transfer rna`.
- `viral_rna` — `hiv|hcv|sars|ires|frameshift|pseudoknot.*viral|tar`.
- `pseudoknot` — `pseudoknot` (not already viral).
- `nmr_solution_motif` — `solution` AND `nmr` AND length < 60 (small isolated motif).
- `loop_motif` — `loop|hairpin|kink.turn|sarcin.ricin|tetraloop` (not already classified).
- `complex_with_protein` — `complex.*(protein|with)` (likely out-of-scope, flag for exclusion).
- `synthetic_designed` — `designed|engineered|synthetic|chimer`.
- `other` — fallback, sent to `uncategorized.csv` for manual review.

**4. Publishability filter heuristics.** Keep classes that are (a) ≥10 sequences across train+val+test (statistical floor), (b) NOT dominated by a single super-famous structure (ribosome may fail this — it's already heavily mined), (c) biologically coherent (riboswitches, ribozymes, pseudoknots are strong here), (d) doable in 3 months without wet-lab data. Most promising a priori: **riboswitches, ribozymes, pseudoknots, viral_rna, nmr_solution_motif**. Less promising: ribosome (over-studied), complex_with_protein (out-of-scope). Final pick deferred to end of Week 1 once rna-code returns counts and rna-research returns one-pagers.

---

## [OPEN] Q-2026-05-18-005 — from rna-code → rna-scientist
**Date:** 2026-05-18
**Topic:** Categorization run complete — next step on focus shortlist
**Question:** The 11-category scheme is now applied across all three Stanford splits (run `01_description_categorization`, see `runs/01_description_categorization/results/SUMMARY.md`). 868 rows processed; **10/11 categories pass the >=10 viability floor** (`ribosome_subunit` 199, `nmr_solution_motif` 131, `viral_rna` 77, `tRNA` 66, `ribozyme` 58, `loop_motif` 52, `complex_with_protein` 43, `riboswitch` 21, `pseudoknot` 11); only `synthetic_designed` fails (n=4); 206 rows fell to `other`. Do you want me to (a) immediately apply the publishability filter (drop over-studied / out-of-scope, weight by val+test presence, drop classes whose val+test contribution is 0) and produce a ranked Week-1 shortlist of 3-5 focus candidates now, or (b) wait for rna-research to deliver the focus_candidate briefs (paper #1/#3/#10 deep summaries) so the shortlist is informed by both data and literature signals? Two minor v2-regex refinements I would queue either way: rescue `Ribonuclease P RNA` -> `ribozyme` and `picornaviral|bacteriophage|coliphage` -> `viral_rna`.
**Context:**
- Results: `/workspace/rna3d/.research/30_experiments/runs/01_description_categorization/results/SUMMARY.md`
- Per-sequence: `results/category_assignment.csv` (868 rows).
- Uncategorized for review: `results/uncategorized.csv` (206 rows).
- Distribution chart: `results/category_distribution.png`.
- Scheme spec: [ANSWERED] Q-2026-05-18-002 in this file.

---

## [OPEN] Q-2026-05-18-006 — from rna-research → rna-scientist
**Date:** 2026-05-18
**Topic:** Focus-candidate briefs ready — handoff for Week-1 ranked shortlist
**Question:** The 7 focus-candidate briefs (per Q-2026-05-18-003) are now committed under `/workspace/rna3d/.research/10_literature/focus_candidates/` — one per priority category: `riboswitch.md`, `ribozyme.md`, `pseudoknot.md`, `viral_rna.md`, `nmr_solution_motif.md`, `tRNA.md`, `loop_motif.md`. Each brief contains: empirical counts (from run 01), biological context with at least one Week-1 reading-plan citation where appropriate, an under-study argument for 3D prediction, a concrete hypothesis seed, follow-up data signals for rna-code, and a publishability-filter verdict (STRONG/MEDIUM/WEAK). Two verdicts came out STRONG (`riboswitch`, `ribozyme`), three MEDIUM (`pseudoknot`, `viral_rna`, `tRNA`), and two WEAK as standalone foci but useful as cross-cutting lenses (`nmr_solution_motif`, `loop_motif`); a short companion section listing all verdicts has been added near the top of `REVIEW.md`. You now have both data (run 01 counts) and literature evidence (these briefs + the Q-002 categorization rationale) to produce the **Week-1 ranked shortlist of 3-5 focus candidates**. Suggested deliverable path: `/workspace/rna3d/.research/20_hypotheses/FOCUS_SHORTLIST.md` (your owned subtree, not mine to write). Once your shortlist is in, the user makes the final pick at the end of Week 1, and rna-code's queued v2 regex refinements (`Ribonuclease P RNA -> ribozyme`, `picornaviral|bacteriophage|coliphage -> viral_rna`, plus pseudoknot re-detection across all categories) can be prioritized accordingly. Anything you want clarified or re-scoped in any brief before you draft the shortlist?
**Context:**
- Briefs: `/workspace/rna3d/.research/10_literature/focus_candidates/*.md`
- Summary section: top of `/workspace/rna3d/.research/10_literature/REVIEW.md`
- Run-01 counts: `/workspace/rna3d/.research/30_experiments/runs/01_description_categorization/results/SUMMARY.md`
- Categorization scheme rationale: [ANSWERED] Q-2026-05-18-002 above.
- rna-code's queued v2 refinements: [OPEN] Q-2026-05-18-005 above.

---
