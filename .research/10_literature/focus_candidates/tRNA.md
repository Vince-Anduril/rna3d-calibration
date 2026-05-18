# Focus candidate: tRNA

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 64, val: 1, test: 1, total: 66
- val+test presence: 2
- viable per Week-1 threshold? yes (>=10 total)

## Biological context (what is already known)

- tRNA is the canonical L-shaped RNA: 4-way junction folding the cloverleaf secondary structure into the acceptor-arm / D-arm / T-arm / anticodon-arm tertiary architecture.
- Hallmark tertiary contacts: D-loop / T-loop interaction, base triple at position 15-48, post-transcriptional modifications (pseudouridine, dihydrouridine, m1A, m2G, ...) at conserved positions critical for fold stability and translation.
- Length is highly stereotyped (~76 nt, plus variable-loop class differences). Massive prior knowledge: thousands of tRNA-related PDB entries, deep Rfam alignments.
- Modified bases are typically **not present** in deposited sequences or in predictor inputs — a known limitation of the prediction setup.
- Biology is central to translation; tRNA mimicry by viral 3' UTRs, tmRNA, and selenocysteine tRNA add functional variants.

## Why it could be under-studied for 3D prediction

- Honestly, tRNA is **arguably over-studied** at the global-fold level — most predictors recover the L-shape comfortably for canonical cytosolic tRNAs.
- Where the field is thinner: prediction of **modified-base influence** on fold (typically ignored), prediction of non-canonical tRNAs (mitochondrial tRNAs with truncated/missing arms, selenocysteine tRNA, tmRNA), and prediction of tRNA-mimicking elements (viral 3' UTRs).
- Mitochondrial tRNAs in particular are structurally degenerate (some lack the D- or T-arm entirely) and we have not verified how DL predictors handle them — flagged as uncertainty.

## What pattern might be worth finding

- Stratify the 66 rows into canonical / mitochondrial / suppressor / tRNA-mimic sub-classes and ask whether prediction error concentrates in the non-canonical sub-classes — would be a sharp claim if true.
- Whether predictors trained on bulk PDB tRNA data over-impose the canonical L-shape on non-canonical inputs (a high-confidence wrong-fold failure mode).
- Quantify whether the absence of modified-base information correlates with systematic local-error patterns at canonical modification sites.

## Data signals rna-code should measure next

- Sub-classify the 66 rows by `description` into canonical/mitochondrial/suppressor/mimic (regex follow-up).
- Length distribution (canonical ~76 nt vs mitochondrial often shorter/variable).
- Per-row check on whether modified bases are reported in the PDB ground truth.
- Per-row inter-model error variance once Week-2 inference lands.

## Publishability filter assessment

(a) Biologically meaningful — **yes**: tRNA is universally essential, and the non-canonical sub-classes are biologically rich. (b) Sufficiently under-studied for 3D prediction — **partial**: the canonical class is well-studied; the **non-canonical sub-classes** are the genuinely under-studied slice, but isolating them requires sub-categorization. (c) Doable in 3 months — **yes**. (d) Val+test presence sufficient — **borderline**: 2 rows, and whether either is a non-canonical sub-class is unknown without re-inspection. **Verdict: MEDIUM candidate** — strong only conditional on a non-canonical sub-population being present in val+test; if all val+test rows are canonical, the story collapses to "models do well on tRNA" which is not novel.
