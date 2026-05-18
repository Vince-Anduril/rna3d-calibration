# Focus candidate: viral_rna

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 75, val: 1, test: 1, total: 77
- val+test presence: 2 (small but non-zero)
- viable per Week-1 threshold? yes (>=10 total)

## Biological context (what is already known)

- Heterogeneous super-class: HIV TAR/RRE, HCV IRES, SARS frameshift element, picornaviral IRES, bacteriophage operator stems (MS2, lambda boxB), influenza panhandle, viral 3' UTRs.
- Structural hallmarks vary by virus family: tRNA-mimicking 3' UTRs (turnip yellow mosaic), tetraloop-bound capsid recognition (MS2), kissing-loop dimerization (HIV DIS), highly structured IRES regions (HCV, picornavirus).
- Functional pressure on structure is intense: viruses use structure for ribosome recruitment, packaging, frameshifting, and replication priming. Conservation of structure can exceed conservation of sequence.
- High biomedical relevance (antiviral target programs, vaccine RNA design), so structural results have a clear translational hook.

## Why it could be under-studied for 3D prediction

- Viral RNAs often have unique or sparsely-populated Rfam families (one virus family per RNA element), so MSA-based methods may have limited co-evolution signal.
- Many viral structural elements are pseudoknotted (frameshift elements, 3' UTR tRNA-like structures) or contain non-canonical pairs, hitting two known prediction blind spots simultaneously.
- The class is **heterogeneous by construction** — a per-class statement risks being a per-virus statement; this is a methodological caution, not a deal-breaker.
- We have not located a recent comprehensive benchmark of DL predictors on viral RNA elements as a stratified class — flagged as uncertainty.

## What pattern might be worth finding

- A within-class stratification: do predictors do well on small bacteriophage stem-loops (MS2, boxB) but fail on functionally complex IRES / frameshift elements? A clean small-vs-complex contrast would be a clear story.
- Specific motif accuracy: kissing-loop dimers (HIV DIS), tRNA-mimic 3' UTRs, frameshift pseudoknots — measure accuracy at the motif level, not at the global RMSD level.
- Whether predictors trained on Ribonanza chemical-mapping data (RibonanzaNet) generalize better to viral elements than MSA-only methods (RhoFold+).

## Data signals rna-code should measure next

- Split the 77 rows by **virus family** via `description` regex (HIV, HCV, SARS, MS2, lambda, etc.) — produces a more honest sub-categorization before any class-wide claim.
- Length distribution within viral_rna — likely bimodal (small operator stems vs full IRES elements).
- Number of unique virus families represented (concentration vs diversity check).
- Per-sub-family inter-model error variance once Week-2 inference lands.

## Publishability filter assessment

(a) Biologically meaningful — **yes**: viral RNA elements drive frameshifting, IRES translation, and packaging — direct biomedical relevance. (b) Sufficiently under-studied for 3D prediction — **yes for sub-classes**, less so for the broad super-class. (c) Doable in 3 months — **yes**, but requires sub-categorization work to avoid making a vague heterogeneous claim. (d) Val+test presence sufficient — **borderline**: 2 rows is very small and almost certainly belong to different sub-families, limiting any single-sub-class statistical claim. **Verdict: MEDIUM candidate** — the heterogeneity is the main risk; would become STRONG only if a specific sub-class (e.g. frameshift pseudoknots) shows strong val+test presence after re-classification, in which case it merges naturally with the pseudoknot lens.
