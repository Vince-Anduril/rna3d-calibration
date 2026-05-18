# Focus candidate: nmr_solution_motif

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 131, val: 0, test: 0, total: 131
- val+test presence: **0** (critical limitation)
- viable per Week-1 threshold? yes (>=10 total) — train-only

## Biological context (what is already known)

- This category is a **methodology bucket**, not a biological family: small (< 60 nt by the scheme rule) RNA fragments solved by solution NMR.
- Typical contents: tetraloops, kissing-loop fragments, kink-turns, sarcin-ricin loops, isolated stem-loops from larger RNAs studied as independent motifs.
- NMR-derived structures are often deposited as **ensembles** (10-20 conformers) rather than single coordinates — a fundamentally different ground truth than crystal or cryo-EM single models.
- Many of these motifs are recurrent across larger RNAs (the tetraloop family in particular is one of the most-studied 4-nt motifs in all of structural RNA biology).

## Why it could be under-studied for 3D prediction

- DL predictors are typically trained against single coordinate sets — an NMR ensemble averaged or arbitrarily picked may introduce label noise that the model can never resolve.
- Conformational heterogeneity is the **point** of many NMR studies (flexible loops, dynamic stems), yet predictors output one structure with a confidence score — the mismatch in epistemological framing is itself the story.
- Small motif structure may suffer from lack of MSA context (these are fragments deliberately taken out of their larger biological context).
- We have not seen a study quantifying DL-predictor accuracy specifically against NMR-ensemble spread (vs single-conformer ground truth) — flagged as uncertainty; this could be a contribution.

## What pattern might be worth finding

- For NMR-ensemble rows in train, compute predictor RMSD against **each** ensemble member and ask: does the predictor consistently land near one cluster, or does its confidence inversely correlate with ensemble spread (well-calibrated uncertainty)? A clean negative result here is publishable.
- Identify whether predicted "low confidence" residues actually correspond to NMR-flagged flexible regions (R1, R2, hetNOE-supported) — direct calibration check against an orthogonal experimental signal.

## Data signals rna-code should measure next

- For each row in this category, count the number of models in the deposited PDB / mmCIF file (1 vs >1 indicates single conformer vs ensemble).
- Compute pairwise intra-ensemble RMSD (a measure of conformational spread per row).
- Cross-reference `target_id` against PDB experimental method tag (`SOLUTION NMR` vs `X-RAY DIFFRACTION`) — sanity check on the regex categorization.
- Per-residue B-factor / ensemble RMSF distribution per row.

## Publishability filter assessment

(a) Biologically meaningful — **partial**: the category is a methodology bucket, not a biological class; biological meaning has to be re-extracted post-hoc (e.g. focus on tetraloops specifically). (b) Sufficiently under-studied for 3D prediction — **yes** at the ensemble-vs-single-structure framing; ensemble-aware evaluation of DL RNA predictors is genuinely under-explored. (c) Doable in 3 months — **yes**, requires PDB ensemble parsing infrastructure. (d) Val+test presence sufficient — **NO**: zero val/test rows kills the standard novel-prediction story on the Stanford benchmark; would require an out-of-benchmark validation set (other NMR-ensemble PDB entries). **Verdict: WEAK candidate as a standalone focus** — the methodological angle is interesting but the val/test gap forces an external dataset, which doubles project scope. Best repurposed as a **methodological appendix** (ensemble-aware evaluation) layered on whichever biological family is chosen.
