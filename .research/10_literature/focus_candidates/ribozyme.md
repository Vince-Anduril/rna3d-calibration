# Focus candidate: ribozyme

> Owner: rna-research. Maintained as new evidence comes in.

## Empirical signal (from run 01_description_categorization)

- train: 54, val: 2, test: 2, total: 58
- val+test presence: 4 (key for novel-prediction angle)
- viable per Week-1 threshold? yes (>=10 total)

## Biological context (what is already known)

- Ribozymes are catalytic RNAs. Classical small ribozymes: hammerhead, hairpin, HDV, VS, twister, pistol, hatchet. Large ribozymes: Group I, Group II self-splicing introns, RNase P, and the ribosomal peptidyl transferase center.
- Structural hallmark: a precisely positioned **catalytic core** stabilized by tertiary contacts (e.g. tertiary stem-loop kissing in hammerhead type-II; ribose-zipper in Group I P4-P6).
- Catalytic activity depends on Angstrom-scale geometry around the scissile phosphate — local error matters more than for non-catalytic RNAs (a 2 Å displacement in the active site is functionally significant).
- Often involve divalent-cation coordination (Mg2+ or Mn2+) — most predictors ignore ions entirely, leaving an explicit modelling gap.
- HDV and twister classes contain pseudoknots at the active site, compounding prediction difficulty.

## Why it could be under-studied for 3D prediction

- Active-site accuracy is rarely separated from global-fold accuracy in standard RMSD/TM metrics — a model can score well overall while missing the catalytic geometry entirely.
- Group I/II introns are long (>200 nt) and conformationally complex; few benchmarks report per-domain or per-active-site accuracy.
- Small ribozymes (hammerhead, hairpin) have crystal structures but represent a small slice of the Rfam ribozyme universe — MSA depth varies dramatically.
- We have not located a recent (post-2024) systematic benchmark of SOTA predictors on ribozyme catalytic-core geometry specifically — flagged as uncertainty; Lilley 2011 is the natural biological reference for catalytic-core geometry.

## What pattern might be worth finding

- Define an "active-site RMSD" (sub-RMSD restricted to a small set of catalytically essential residues per ribozyme class) and show that SOTA predictors systematically over-confidence this region — predicted pLDDT high, local RMSD high.
- Show that pseudoknot-containing ribozymes (HDV, twister) suffer disproportionately versus pseudoknot-free (hairpin) — links a structural-motif gap to a functional-class gap.
- Test whether MSA depth correlates with active-site accuracy: do small Rfam families lose more than global metrics suggest?

## Data signals rna-code should measure next

- Per-row Rfam-class label (RNase P RF00010, hammerhead RF00163, HDV RF00094, Group I RF00028, etc.) — Week-2 cross-reference.
- Length distribution (small ribozymes <100 nt vs Group I/II ~300+ nt) — splits the class into two regimes for analysis.
- Active-site residue list per Rfam class (curated from PDB literature) — needed before any "active-site RMSD" metric.
- Inter-model error variance on the 4 val+test rows once inference is online.

## Publishability filter assessment

(a) Biologically meaningful — **yes**: catalysis is the defining function of these RNAs and active-site geometry is the load-bearing structural feature. (b) Sufficiently under-studied for 3D prediction — **likely yes**: active-site-stratified accuracy is rarely reported; we have not confirmed a comprehensive post-2024 benchmark, flagged as uncertainty. (c) Doable in 3 months — **yes**, but requires per-class active-site annotation work (extra week of literature mining vs riboswitches). (d) Val+test presence sufficient — **borderline**: 4 rows; useful as case studies, not for class-wide statistics. **Verdict: STRONG candidate** — the active-site-accuracy angle is a sharp, falsifiable, biologically motivated story even at small n.
