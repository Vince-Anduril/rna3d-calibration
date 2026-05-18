# Paper outline (draft v0.1 — pre-data)

> **Status:** speculative outline, written before the first real GPU run. Sections
> with TBD-data are placeholders that will be filled once the eval run completes.
> Sections marked **TO VERIFY** require a literature check before the paper is sent
> out (some claims are from my training knowledge and need confirmation from
> primary sources during the actual paper-writing phase).

---

## Working title

*"Confidence calibration of a custom RNA structure predictor on ribozyme
ground-truth structures: a focused analysis of the Stanford RNA 3D Folding
benchmark."*

(Backup title: "A small predictor, a focused class: characterizing
prediction-error calibration in ribozyme structures.")

## Abstract sketch (~150 words)

The Stanford RNA 3D Folding Kaggle benchmark provides ground-truth structures
across diverse RNA classes, enabling fine-grained study of where structure
prediction is well-grounded and where it is not. We focus on ribozymes — a
class whose catalytic function makes active-site geometry the structurally
load-bearing feature — and ask a calibration question: when a predictor reports
high per-residue confidence in a ribozyme catalytic core, how reliable is that
confidence? We train a small, deliberately controllable transformer-based
predictor (~5M parameters) from scratch on the public training split, then
analyze its calibration on the validation and test splits with particular
attention to the four ribozyme structures present. We find [TBD-data]. Our
contribution is not a SOTA accuracy claim but a focused, reproducible
calibration characterization on a biologically central RNA class — a
complementary angle to recent large-model RNA structure prediction work.

## 1. Introduction (~600 words)

- Para 1: RNA structure prediction has advanced rapidly post-AlphaFold2 (AF3
  RNA, RhoFold+, RibonanzaNet, Boltz-1, Chai-1). **TO VERIFY** — confirm
  current SOTA list and a recent benchmark/survey.
- Para 2: Calibration of confidence has been extensively studied for protein
  structure (pLDDT — Mariani 2013, Tunyasuvunakool 2021 AlphaFold pLDDT
  follow-ups). Equivalent systematic analysis for RNA is less developed.
  **TO VERIFY** — search recent (2024-2025) RNA calibration literature.
- Para 3: We narrow scope to *ribozymes* because (a) catalytic function makes
  active-site geometry the load-bearing structural feature, (b) the class is
  well-bounded with established catalytic families (hammerhead, hairpin,
  group I/II intron, RNase P, glmS, HDV, twister, pistol, hatchet, **TO
  VERIFY** completeness from Lilley 2011 or more recent review).
- Para 4: We deliberately use a small (~5M-param) custom predictor rather than
  attempting to deploy SOTA models, for two reasons: (i) we control the
  confidence channel by design; (ii) we can analyze and reproduce the entire
  pipeline within a single GPU session. This trades raw accuracy for
  transparency.
- Para 5: Contributions:
  - A focused calibration analysis on ribozyme structures from the Stanford
    benchmark, stratified by class.
  - A small, fully open custom predictor whose calibration channel is
    explicit by design.
  - A finding **[TBD-data]** about where confidence and realized error align
    or diverge for this class.

## 2. Related work (~500 words)

- 2.1 RNA 3D structure prediction (AF3, RibonanzaNet-3D variants, RhoFold+,
  vfold, trRosettaRNA — **TO VERIFY** which exist + are open).
- 2.2 Calibration of structure predictors (pLDDT for proteins — well-studied;
  RNA confidence — **TO VERIFY** recent papers).
- 2.3 Conformal prediction and post-hoc calibration in scientific ML (Guo
  2017 temperature scaling; Angelopoulos 2021 conformal tutorial).
- 2.4 Ribozyme structural biology and the catalytic-core problem (Lilley
  2011, Doudna 2002 group I, **TO VERIFY** recent reviews 2020+).

## 3. Data and methods (~800 words)

### 3.1 Dataset

Stanford RNA 3D Folding benchmark (Kaggle competition closed 2025-09-24).
868 sequences across train (844) / val (12) / test (12). Sequence lengths
3-4298 nt. PDB ground-truth structures for 8,672 PDB files (some sequences
have multiple PDB entries / chains).

### 3.2 Family categorization

We categorize sequences from the free-text `description` field into 11
biological classes via priority-ordered regex rules (run 01b in our
pipeline). Ribozyme class: **[TBD-counts]** train, val, test sequences
after the v2 rescue pass.

### 3.3 Model

Small transformer encoder:
- Vocabulary: {[PAD], [MASK], [CLS], A, U, G, C, N} (size 8)
- 6-layer transformer, d_model=256, 8 heads, FF=1024 (~5M parameters)
- Two heads on top of the encoder: per-residue (x, y, z) coordinate head
  (3-layer MLP) and per-residue confidence head (3-layer MLP + sigmoid).

### 3.4 Training

Two-stage:
- **Stage 1** (~30 min, GPU): masked-LM pretraining on the 5,135 unique
  sequences in v2 train, mask probability 0.15.
- **Stage 2** (~2.3h, GPU): structure supervision on the ~800 train
  sequences with PDB ground-truth. Joint loss:
  L = L_coord + 0.5 × L_pair + 0.5 × L_conf
  where:
  - L_coord = MSE on per-residue (x, y, z) after per-batch Kabsch
    alignment.
  - L_pair = MSE on 64 random pair-distances per sample (auxiliary).
  - L_conf = BCE between predicted confidence and target =
    exp(-||pred - true|| / 5Å) computed on the same Kabsch-aligned
    coordinates.

### 3.5 Evaluation

- Per-sequence: Kabsch-aligned RMSD over all valid residues.
- Per-residue: confidence × realized error pairs → reliability diagrams +
  ECE/MCE/Brier score (Guo 2017).
- Stratified by category (ribozyme highlighted).

## 4. Results (~800 words) — TBD-data

(Fill from `.research/30_experiments/runs/GPU_training_eval/SUMMARY.md`)

- 4.1 Global accuracy and calibration (val + test, all categories).
- 4.2 Per-class breakdown (Table 1).
- 4.3 Ribozyme deep-dive: per-residue reliability diagram (Figure 2),
  catalytic-core sub-RMSD (Figure 3, if active sites are identifiable in
  ground truth).
- 4.4 Comparison to a trivial-confidence baseline (confidence = constant).
- 4.5 Sensitivity analysis: how do results change with bin granularity,
  random seed, max sequence length truncation.

## 5. Discussion (~600 words)

- What our calibration finding suggests about predictor over-confidence
  on ribozyme active-site geometry.
- Why a small custom model has methodological value even when SOTA is
  more accurate.
- Limitations:
  - Small custom model is not SOTA — we explicitly do NOT claim
    competitive accuracy.
  - Small val/test sets (n=4 ribozymes total) — we apply bootstrap CIs
    and present per-sequence data.
  - No multimer support, no RNA-protein complexes.
  - Single MSA strategy (we used the v1 MSAs; v2 deferred).
- Future work: extend to riboswitch ON/OFF conformer pairs (alternate
  focus from our Sprint-0 shortlist); apply same analysis to RhoFold+ or
  Boltz-1 predictions when those are available.

## 6. Reproducibility

All code at https://github.com/Vince-Anduril/rna3d-calibration. Single 3-hour
RTX 5090 pod session reproduces the entire trained model + eval. Random seeds
fixed. PDB extraction snapshot dates noted in `data/MANIFEST.md`.

## Figures (planned)

1. **Pipeline schema** — sequences + PDB → small transformer → coord + conf.
2. **Reliability diagram** — global vs ribozyme-only (the key figure).
3. **Active-site overlay** — predicted vs ground-truth for one representative
   ribozyme (e.g., the RNase P entry if present in val/test).
4. **Confidence-vs-error scatter** — per-residue, colored by category.
5. (Supplementary) **Loss curves** + **per-class RMSD distribution**.

## Tables (planned)

1. **Per-category eval** — n, median RMSD, mean confidence, ECE per category.
2. **Ribozyme sub-class breakdown** — per ribozyme family (hammerhead,
   hairpin, etc., if represented).

## Venue plan

Following the spec (Section 6.2 of the original design):
- Preprint: bioRxiv first
- Primary submission: **Bioinformatics** (Oxford UP) — short-paper format fits
- Backup: BMC Bioinformatics; NeurIPS LMRL / MLSB workshops if timing aligns

## Risks / counter-arguments to anticipate

- **"You trained a small model — your calibration findings are model-specific."**
  Response: yes, by design. Our claim is about characterizing calibration on a
  focused class with a controllable predictor, not about ranking SOTA. The
  same analysis pipeline can be applied to RhoFold+/Boltz-1 predictions when
  those are accessible — explicit future work.
- **"n=4 ribozyme val+test sequences is too small."**
  Response: we present per-sequence data + bootstrap CIs on aggregate. The
  paper does not over-claim. We also analyze ribozyme sequences in the train
  set (n=~60) for the structural patterns themselves (acknowledging the
  train/test distinction).
- **"This is not novel — calibration is well-studied."**
  Response: not for RNA structure prediction, and not for ribozymes
  specifically. **TO VERIFY** this claim during writing — if a 2024-2025 paper
  has done exactly this we either cite + differentiate or pivot.
