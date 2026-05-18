# Focus candidate: pseudoknot

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 11, val: 0, test: 0, total: 11
- val+test presence: **0** (critical limitation for a novel-prediction angle on the Stanford benchmark)
- viable per Week-1 threshold? yes (>=10 total) — but only marginally and only on train

## Biological context (what is already known)

- A pseudoknot is a secondary-structure motif where bases of a hairpin loop pair with sequence outside the hairpin, producing crossed (non-nested) base pairs.
- Functional roles include programmed -1 ribosomal frameshifting (HIV, SARS-CoV-2), telomerase template anchoring, viral IRES elements, and core architecture of many ribozymes (HDV, twister) and riboswitches (preQ1, SAM-II).
- Hallmark structural features: coaxial helical stacking, kissing-loop tertiary contacts, often single-Mg2+-dependent stability.
- The non-nested base-pair topology breaks the standard dynamic-programming assumption of nested folds — most classical SS predictors (Mfold, RNAfold default) cannot predict pseudoknots without special algorithms (DotKnot, IPknot).
- n=11 train is small but reflects an **under-studied class** rather than absence of biology — pseudoknots are everywhere in viral and regulatory RNAs but PDB-deposited pure-pseudoknot structures are scarce.

## Why it could be under-studied for 3D prediction

- MSA-based predictors inherit secondary-structure biases from co-evolution signals, which are often weaker across non-nested pairings.
- The pseudoknot literature has been actively pointing out this blind spot for years (Reuter & Mathews, Sloma & Mathews) — yet post-2024 DL predictors have not (to our knowledge) been systematically benchmarked specifically on pseudoknot accuracy as a stratified metric — flagged as uncertainty.
- Pseudoknots are typically embedded inside larger RNAs in real biology (frameshift elements, ribozyme cores) rather than appearing as standalone deposits — categorization on `description` alone undercounts them.

## What pattern might be worth finding

- Even **without** Stanford val/test rows, the n=11 train rows could be used to study what predictors learn versus miss at the pseudoknot crossing — e.g. recover pseudoknot base pairs from predicted 3D coordinates and measure base-pair recall against ground truth.
- Cross-class signal: identify pseudoknot-containing rows inside other categories (riboswitch, ribozyme, viral_rna) via secondary-structure detection on PDB ground-truth — this is the realistic path forward given the zero val/test presence.

## Data signals rna-code should measure next

- Run a pseudoknot detector (e.g. extract base pairs from PDB ground truth, classify as nested vs crossed) across **all** Stanford rows, not just the n=11 — produces a `has_pseudoknot` boolean column that subsumes the regex category.
- Per-row count of crossed base pairs (intensity of pseudoknot character).
- Once Week-2 inference is up, base-pair recall stratified by nested vs crossed pairs per row.

## Publishability filter assessment

(a) Biologically meaningful — **yes**: pseudoknots are functional units in viruses, telomerase, and ribozymes. (b) Sufficiently under-studied for 3D prediction — **yes** at the analysis-stratification level, though the pseudoknot blind spot is itself a well-known meta-claim — the novelty would have to come from a sharper per-residue or per-motif analysis, not from re-asserting the blind spot. (c) Doable in 3 months — **yes**, but requires a pseudoknot-detection pipeline. (d) Val+test presence sufficient — **NO at the categorical level**: zero val/test rows kill the direct novel-prediction story. The cross-class redefinition (pseudoknots-inside-other-categories) is the only viable angle. **Verdict: MEDIUM candidate** — strong biology, but the categorical signal is too weak; best deployed as a **secondary lens** on a richer focus family (riboswitch or ribozyme), not as the primary focus.
