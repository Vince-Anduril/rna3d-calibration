# Focus candidate: loop_motif

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 50, val: 1, test: 1, total: 52
- val+test presence: 2
- viable per Week-1 threshold? yes (>=10 total)

## Biological context (what is already known)

- This category is a **motif bucket**: rows whose description prominently mentions loop, hairpin, kink-turn, sarcin-ricin, tetraloop, etc., not already absorbed by a higher-priority biological class.
- Tetraloops (GNRA, UNCG, CUUG) are the most-studied small RNA motifs, with strong sequence-structure rules and recurrent appearance across rRNA, ribozymes, and riboswitches.
- Kink-turns introduce a ~120 degree bend, common in box C/D snoRNAs, riboswitches (SAM, lysine), and ribosome.
- Sarcin-ricin loop in 23S rRNA is the target of ribotoxins and is highly structurally conserved.
- These motifs are typically substructures within larger RNAs, so this category captures isolated/fragment depositions of well-characterized motifs.

## Why it could be under-studied for 3D prediction

- As **isolated fragments**, these rows lose the structural context that normally constrains them in their host RNA — predictors may struggle differently on fragments vs in-context motifs.
- Tetraloop accuracy at the per-residue level has rarely been reported as a stratified DL-predictor benchmark to our knowledge — flagged as uncertainty.
- The category is heterogeneous (multiple motif types collapsed), reducing the sharpness of any class-wide claim.

## What pattern might be worth finding

- Sub-classify rows by motif type (tetraloop / kink-turn / sarcin-ricin / other) and measure per-motif accuracy; would expect tetraloops to be best-predicted and any sharp negative finding (e.g. kink-turn angle systematically mispredicted) would be publishable.
- Test fragment-vs-in-context: when the same motif appears both as a deposited fragment and as a substructure inside a larger ribosome/ribozyme deposit, do predictors agree?

## Data signals rna-code should measure next

- Sub-classify the 52 rows into motif types via secondary regex pass on `description`.
- Length distribution (tetraloops are tiny ~4-12 nt, sarcin-ricin ~30 nt).
- For each motif type, count whether a matched in-context occurrence exists in other categories.
- Per-row inter-model error variance once Week-2 inference lands.

## Publishability filter assessment

(a) Biologically meaningful — **partial**: motifs are biologically important but the category as currently defined is heterogeneous; meaning would need to come from sub-classifying to one specific motif type. (b) Sufficiently under-studied for 3D prediction — **partial**: tetraloops in particular are well-studied at the sequence-structure level; the gap is at the per-residue DL-predictor stratification, which is real but narrow. (c) Doable in 3 months — **yes**, but requires sub-classification and a clean motif-detection pipeline. (d) Val+test presence sufficient — **borderline**: 2 rows, almost certainly different motif types, limiting per-motif statistics. **Verdict: WEAK candidate as a standalone focus** — too heterogeneous, val/test too thin per sub-type. Best repurposed as a **cross-cutting motif-accuracy lens** layered onto a biological focus family (e.g. measure kink-turn accuracy within riboswitches).
