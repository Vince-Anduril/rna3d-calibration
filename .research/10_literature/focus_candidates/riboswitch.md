# Focus candidate: riboswitch

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 17, val: 2, test: 2, total: 21
- val+test presence: 4 (key for novel-prediction angle)
- viable per Week-1 threshold? yes (>=10 total)

## Biological context (what is already known)

- Riboswitches are cis-acting non-coding RNA elements in 5' UTRs of mRNAs (predominantly bacterial) that bind small-molecule metabolites and gate downstream expression — e.g. TPP, FMN, SAM, cobalamin, purine, glycine, lysine, c-di-GMP classes [McCown 2017].
- Structural hallmark: a conserved **aptamer domain** (the ligand-binding pocket, typically a 3-way or 4-way junction) coupled to a **variable expression platform** (terminator, anti-terminator, sequestrator) whose fold flips upon ligand binding.
- Many families show **two mutually exclusive ground-state folds** (ON / OFF), meaning a single PDB structure represents only one functional conformer — the biology is intrinsically multi-state.
- Biological importance is large: riboswitches are antibiotic targets, gene-regulation models, and biosensor chassis. The Breaker lab catalogued 40+ validated classes with hundreds of candidate orphan classes still under bioinformatic investigation.
- Aptamer tertiary contacts (kissing loops, pseudoknot stems, A-minor motifs) are highly idiosyncratic per class — generic MSA features may not transfer well across riboswitch families.

## Why it could be under-studied for 3D prediction

- The "one sequence, one structure" supervised setup mis-represents the conformational duality of riboswitches — predictors are trained against whichever conformer happens to be in the PDB.
- Many riboswitch classes have small Rfam alignments (orphan or recently discovered), limiting MSA-driven methods (RhoFold+, AF3 RNA path).
- Riboswitches are under-represented in CASP-RNA targets (CASP15 had only a handful of pure-RNA targets, and to our knowledge riboswitches were a small subset); we have not verified post-CASP15 community benchmarks specifically for riboswitches — flagged as uncertainty.
- Pseudoknot-containing aptamer cores (e.g. preQ1, SAM-II) compound the difficulty: pseudoknots are a known blind spot of MSA-based predictors.

## What pattern might be worth finding

- Quantify whether predictor confidence (pLDDT-equivalent) is **mis-calibrated specifically at the ligand-binding pocket** of riboswitch aptamers — if confidence is high but RMSD is high in those exact residues, that is a publishable, specific finding.
- Test whether models systematically predict only one conformer regardless of context (e.g. always the ligand-bound form when training PDB skews that way) — measurable as per-residue error concentrated in the "switching helix" region.
- Stratify error by riboswitch class to see whether one specific class (e.g. SAM, TPP, cobalamin) dominates the residual error budget — under-studied classes within the under-studied family.

## Data signals rna-code should measure next

- Per-class sequence-length distribution and pairwise sequence identity within the 21 riboswitch rows (sanity check: are these mostly the same class or diverse?).
- Cross-reference `target_id` against Rfam family IDs to label each row by riboswitch class (RF00050 FMN, RF00059 TPP, etc.) — Week-2 deliverable.
- Per-sequence inter-model error variance (RibonanzaNet vs RhoFold+) on the 4 val+test riboswitch rows once Week-2 inference is online.
- MSA depth per row (from any available alignment artifact in the Stanford bundle, or recomputed with rMSA / Infernal on Rfam seed).

## Publishability filter assessment

(a) Biologically meaningful — **yes**: riboswitches are central to bacterial gene regulation and a long-standing antibiotic-target program. (b) Sufficiently under-studied for 3D prediction — **likely yes**, with caveat: scientist literature (McCown 2017) confirms classical structural diversity; we have not verified post-2024 benchmarks specifically targeting riboswitches as a class. (c) Doable in 3 months without wet-lab — **yes**: pure inference + analysis of existing PDB ground truth. (d) Val+test presence sufficient — **borderline**: 4 rows is small but non-zero, enough for case-study analysis but not for statistical claims. **Verdict: STRONG candidate** (high biological coherence and conformational-duality angle compensates for small val+test n).
