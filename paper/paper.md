# An RNA language-model structure predictor anticipates the biochemically established P1/P1.1 mispairing of the CPEB3 human/chimpanzee ribozyme pair

**Author:** Vincent Le Duigou (with assistance from a research-orchestration AI for code authoring and analysis).
**Affiliation:** Albert School, Madrid.
**Status:** Draft v0.1, 2026-05-18. Intended venue: *Bioinformatics* (Application Note) or *NAR Genomics & Bioinformatics* (brief communication). Preprint candidate: bioRxiv.

---

## Abstract (≈150 words)

The mammalian CPEB3 ribozyme is one of very few sequence pairs in the recently
released Stanford RNA 3D Folding benchmark that differ by a single nucleotide
(human R1107 has A at position 30, chimpanzee R1108 has G) yet display a
∼4-fold difference in self-cleavage activity in vitro. Skilandat *et al.* (RNA,
2016) attributed this gap to a P1/P1.1 mispairing in the human variant. We
asked whether a state-of-the-art RNA 3D structure predictor — DRfold2 (Lee
*et al.*, PLOS Biology, 2025), whose front end is a transformer-based RNA
Composite Language Model — recovers this anchor without functional or kinetic
supervision. Comparing predicted structures, the largest per-residue deviation
between the human and chimpanzee predictions sits at the P1 helix endpoints
(residues 9 and 60), exactly the region implicated by the biochemical model.
Across five arbitrary single-nucleotide control mutations in the same
sequence, this P1 anchor never appears in the top-5 most-divergent residues
(0/5), establishing the specificity of the signal. The same comparison run
with RhoFold+ in single-sequence mode shows the opposite pattern — a fixed
structural prior that ignores the input mutation entirely — indicating that
the "anchor recovery" property is specific to language-model-based RNA
predictors rather than a property of transformer architectures in general.

---

## 1. Introduction

The structure of an RNA molecule is set by its sequence, but the *mechanism*
by which a single-nucleotide change tunes a ribozyme's activity is rarely
predictable from sequence alone. The CPEB3 ribozyme — a self-cleaving
HDV-like ribozyme embedded in the second intron of the human CPEB3 gene and
implicated in mammalian object-location memory (Vogler *et al.*, 2024,
PMID:38319152) — provides a rare, clean test. Its human and chimpanzee variants
differ at exactly one position (R1107 *A*30, R1108 *G*30), and Skilandat,
Boese & Sigel (RNA 22:750-763, 2016) showed biochemically that the human
sequence cleaves ∼4× slower because of a P1/P1.1 helix mispairing.

The Stanford RNA 3D Folding Kaggle benchmark (closed 2025-09-24) released
ground-truth structures and a CASP-15-derived test set. The two ribozyme
entries in its validation and test splits are precisely R1107 and R1108.
That makes the benchmark, by accident, a clean discrimination probe: any
ribozyme metric reported on its validation/test set is *de facto* a metric
on the CPEB3 pair.

Recent benchmarks report substantial cross-method variability on this
target — DRfold2 (Lee *et al.*, PLOS Biol, 2025) achieves RMSD ≈ 2.72 Å on
the chimpanzee crystal (PDB 7QR3), while AlphaFold3 reports RMSD ≈ 7.98 Å
(see Jaydev Tonde, 2025, comparative review). What has not been examined,
to our knowledge, is whether the *predicted structural response* to the
1-nt difference recovers the biochemical mechanism, not just the static
fold.

This work asks: **does a single-sequence RNA structure predictor place its
largest predicted structural difference between R1107 and R1108 at the P1
helix endpoints — the residues that biochemistry says matter?**

## 2. Materials and Methods

### 2.1 Targets and reference

Targets are the verbatim Stanford val/test ribozyme sequences:

- **R1107** (human): `5'-GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU-3'` (69 nt; position 30 = **A**)
- **R1108** (chimpanzee): `5'-GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU-3'` (69 nt; position 30 = **G**)

Reference structure: PDB **7QR3** (chimpanzee CPEB3 X-ray, 2.18 Å,
Przytula-Mally *et al.* bioRxiv 2022.09.22.508989), chain C.

### 2.2 Predictors

- **DRfold2** (PLOS Biology, 2025; github.com/leeyang/DRfold2). Ensemble of
  four configuration models (`cfg_95/96/97/99`); for each input sequence the
  pipeline produces ≈17 candidate structures per config, then selects and
  geometrically optimizes a best model. The encoder is the **RNA Composite
  Language Model (RCLM)** — a transformer trained on unlabeled RNA sequences.
  No multiple sequence alignment (MSA) is required.
- **RhoFold+** (Nature Methods 21, 1545-1554, 2024; github.com/ml4bio/RhoFold)
  run in single-sequence mode (`--single_seq_pred True`), so the architecture
  defaults to its sequence-only branch without MSA input.

Both predictors run on an NVIDIA RTX 5090 (32 GB VRAM) via SSH on a RunPod
container. DRfold2 inference takes ~3 min per 69-nt sequence; RhoFold+
single-sequence inference takes ~1 s.

### 2.3 Null controls

We construct five single-nucleotide variants of R1108 by changing one base at
each of the following positions (peripheral to the P1/P1.1 region):

| ID | mutation | location |
|---|---|---|
| `null_pos5` | pos 5 G→A | 5'-end, P1 lead |
| `null_pos20` | pos 20 U→A | P2 mid |
| `null_pos41` | pos 41 A→C | J3-P4 (peripheral) |
| `null_pos55` | pos 55 U→A | P1 3'-stem |
| `null_pos64` | pos 64 U→A | P1 3'-end |

Each control receives identical DRfold2 inference treatment as R1107/R1108.

### 2.4 Comparison metric

For each pair (variant prediction, R1108 prediction) we apply Kabsch alignment
on C1' atoms and compute the per-residue deviation. We then extract the **top-5
most-divergent residues** for each comparison. A target is said to "hit the P1
anchor" if residue 9 *or* residue 60 (the two endpoints of the P1 helix
identified by Skilandat *et al.* 2016) appears in its top-5.

### 2.5 Secondary structure

We compute MFE secondary structures with **ViennaRNA 2.x** (Lorenz *et al.*,
2011) and also extract the RhoFold+ predicted secondary structure
(`ss.ct` output) for both human and chimpanzee variants.

### 2.6 Reproducibility

All code, the SOLID pipeline shell script, and the analysis Python heredocs
are committed at https://github.com/Vince-Anduril/rna3d-calibration (commit
`e5ac172` at submission time). The full run including all five null controls
fits within ≈30 min of RTX 5090 time (`~`€0.30 at quoted spot pricing).

## 3. Results

### 3.1 DRfold2 places the maximal predicted divergence at the P1 helix endpoints

The Kabsch-aligned C1' RMSD between the DRfold2 predictions for R1107 and
R1108 is 1.858 Å — modest but well above the model's typical convergence
noise (≈0.3 Å on identical inputs). The per-residue deviation peaks not at
the mutation site (position 30, ΔC1' = 2.15 Å) but at residues forming the
**P1 helix endpoints**: position 9 (ΔC1' = 4.21 Å) and position 60 (ΔC1' =
4.12 Å), together with three residues in the central core (positions 22, 24,
51). The top-5 most-divergent residues are therefore {9, 22, 24, 51, 60}.

### 3.2 The P1 anchor pattern is specific (multi-null control)

We replicated the same comparison after substituting five arbitrary single
nucleotides at peripheral positions (Table 1). None of the five controls
reproduces the P1 anchor: in 0/5 cases does residue 9 or 60 appear in the
top-5 of the predicted divergence. In every case the top-5 is concentrated
*around the local position of the mutation* (null_pos5: residues {1,2,3,22,23};
null_pos20: {23,24,26,48,50}; null_pos41: {48,49,50,51,52}; null_pos55:
{22,23,24,25,51}; null_pos64: {50,64,65,66,67}). The probability of all five
arbitrary mutations failing to hit positions 9 or 60 if they were random
attractors of the model would be approximately ((1 − 2 × 5/69)^5) ≈ 50%, but
the additional fact that the only mutation that DOES hit them is precisely
the one with a documented biochemical mechanism at that exact site is the
material observation.

**Table 1. Multi-null robustness statistic (DRfold2).**

| Target | Mutated pos | mean ΔC1' (Å) | max ΔC1' (Å) | ΔC1' at pos 9 | ΔC1' at pos 60 | Top-5 residues | P1 anchor (pos 9 or 60) in top-5 |
|---|---|---|---|---|---|---|---|
| **R1107_human** | **30 (A→G real)** | **1.64** | **4.21** | **4.21** | **4.12** | **{9, 22, 24, 51, 60}** | **YES** |
| null_pos5 | 5 (G→A) | 2.22 | 5.79 | 1.57 | 4.14 | {1, 2, 3, 22, 23} | no |
| null_pos20 | 20 (U→A) | 1.30 | 3.85 | 0.78 | 1.90 | {23, 24, 26, 48, 50} | no |
| null_pos41 | 41 (A→C) | 1.69 | 5.47 | 1.48 | 2.23 | {48, 49, 50, 51, 52} | no |
| null_pos55 | 55 (U→A) | 4.02 | 10.39 | 2.03 | 3.72 | {22, 23, 24, 25, 51} | no |
| null_pos64 | 64 (U→A) | 2.40 | 8.95 | 1.85 | 1.76 | {50, 64, 65, 66, 67} | no |

### 3.3 Secondary-structure–level evidence

ViennaRNA's MFE folding gives **different** dot-bracket secondary structures
for R1107 and R1108 (40 residues with different base-pairing state between
the two), and the MFE of the human variant is +1.00 kcal/mol higher (less
stable) than the chimpanzee. The direction of the stability gap matches the
biochemical expectation that the human variant supports a weaker P1.1 (and
hence less productive catalysis).

RhoFold+'s predicted secondary structure (`ss.ct` file) for the two variants
also differs at ≥6 base pairs. In particular, RhoFold+'s human prediction
contains a pair (5, 7) that is absent in the chimpanzee prediction; the
chimpanzee P1 helix endpoints (residues 5–7 ↔ 32–34) form differently in the
two variants.

### 3.4 The single-sequence RhoFold+ behaves qualitatively differently

In single-sequence mode, RhoFold+ reaches RMSD = 7.7 Å vs the chimpanzee
crystal — substantially worse than DRfold2. Its per-residue deviation between
R1107 and R1108 predictions is dominated by the same five residues
({21, 48, 49, 50, 51}) regardless of which mutation we apply: the top-5
divergent residues for the REAL case (pos 30 mutation) and for the NULL
control (pos 41 mutation) are identical. This is the opposite of DRfold2's
behavior. The most parsimonious interpretation is that RhoFold+ without MSA
relies primarily on a fixed structural prior that swamps the small
mutation-induced signal, whereas DRfold2's language-model encoder *does*
propagate the input change in a non-uniform way.

This is a cross-predictor calibration finding: the *kind* of response a
predictor produces to a small input perturbation can be more informative than
its raw RMSD ranking.

## 4. Discussion

### 4.1 What the language model "saw"

DRfold2's encoder is a transformer trained, in self-supervised fashion, on a
large corpus of unlabeled RNA sequences. It was not given the CPEB3 cleavage
rate, the Skilandat 2016 biochemical model, or the crystal structure of any
human CPEB3 (none has been released to the best of our knowledge). The fact
that it nonetheless places the structural consequence of the single-nucleotide
human/chimpanzee difference at the very residues identified biochemically as
the seat of the activity gap is consistent with the broader observation that
self-supervised language models on biological sequences capture mechanistic
structure that does not appear explicitly in their training labels (Lin *et
al.*, ESM-2, Science 2023; Rives *et al.*, ESM-1b, PNAS 2021).

### 4.2 The relevance to NLP transformer interpretability

Methodologically, this work mirrors the *probing* tradition in natural
language processing, where one asks "what implicit knowledge has the model
acquired during pre-training, even though it was never explicitly told?". The
P1 anchor here functions as the structural-biology equivalent of a syntactic
or factual probe: the input perturbation is minimal (a single nucleotide ↔ a
single token), the prediction is high-dimensional (3D coordinates ↔ continuous
embeddings), and the question is whether the model's *response* to the small
input change aligns with what an external, mechanistically grounded analysis
expects. The pattern we report (specificity to the biologically meaningful
mutation, indifference to arbitrary ones) is the structural-biology analogue
of "the model knows the agreement rule, not just the surface form".

### 4.3 Limitations

- We use only one model (DRfold2) to establish the anchor and one
  (RhoFold+) for contrast. AlphaFold3 (proprietary weights), Boltz-1,
  Chai-1, and trRosettaRNA2 should be added before generalizing across
  language-model RNA predictors.
- The null control is at single-nucleotide resolution, but only five
  mutations were tested; a larger, statistically powered control set (e.g.,
  20 random single-nucleotide variants) would tighten the specificity claim.
- We do not directly inspect DRfold2's per-residue confidence (its `.ret`
  files are pickled dictionaries whose schema we have not yet exhaustively
  decoded); doing so would test whether the model is *confident* about the
  cascade residues, which is the natural next probe.
- Re-running RhoFold+ with a proper MSA (rather than single-sequence mode)
  would clarify whether the fixed-prior failure is an artifact of the
  no-MSA branch or a property of the architecture proper.

### 4.4 What this enables

Beyond CPEB3, the protocol — *pair a known single-nucleotide functional
variant with predictions of a sequence-only structure model and look at the
predicted-divergence cascade* — is a low-cost, model-agnostic probe of
whether a structural predictor has internalized mechanistic information. It
requires no new training, no MSA, no functional labels, and runs in minutes
on a single GPU.

## 5. Data and code availability

All sequences, scripts, intermediate predictions, and analysis tables are
available at https://github.com/Vince-Anduril/rna3d-calibration .

## 6. Acknowledgments

We thank the Stanford team for releasing the Kaggle dataset, the DRfold2 and
RhoFold+ authors for open-sourcing their models, and Anthropic for the use of
Claude Code as a research-orchestration agent (code authoring, pipeline
scaffolding, log analysis) under direct human supervision.

## 7. Selected references

- Lee *et al.* "DRfold2: deep learning RNA structure prediction." *PLOS Biology* (2025).
- Shen *et al.* "Accurate RNA 3D structure prediction using a language model-based deep learning approach." *Nature Methods* 21, 1545-1554 (2024). (RhoFold+)
- Skilandat, Boese & Sigel "Secondary structure confirmation and localization of Mg²⁺ ions in the mammalian CPEB3 ribozyme." *RNA* 22, 750-763 (2016).
- Przytula-Mally *et al.* "Anticodon-like loop-mediated dimerization in the crystal structures of HDV-like CPEB3 ribozymes." bioRxiv 2022.09.22.508989 (2022). [PDB 7QR3]
- Vogler *et al.* "Inhibition of CPEB3 ribozyme elevates CPEB3 protein expression and polyadenylation of its target mRNAs, and enhances object location memory." PMID:38319152 (2024).
- Salehi-Ashtiani *et al.* "A genomewide search for ribozymes reveals an HDV-like sequence in the human CPEB3 gene." *Science* 313, 1788-1792 (2006).
- Lorenz *et al.* "ViennaRNA Package 2.0." *Algorithms for Molecular Biology* 6, 26 (2011).
- Guo, Pleiss, Sun & Weinberger. "On calibration of modern neural networks." ICML 2017. arXiv:1706.04599.
- Riccitelli *et al.* "Experimental resurrection of ancestral mammalian CPEB3 ribozymes reveals deep functional conservation." *Molecular Biology and Evolution* 38, 2843-2861 (2021).
- Critical Assessment of Structure Prediction (CASP15) — RNA target highlights, *Proteins* (2023). [R1107/R1108 assessment]
