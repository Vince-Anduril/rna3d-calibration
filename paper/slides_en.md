# Slides — 5-min presentation (EN)

> Format: Markdown, render with `marp` or `pandoc -t beamer`.
> 7 slides, ~40 s each.
> Audience: NLP / transformers course project — short talk.

---

## Slide 1 — Title

**Does an RNA transformer "see" function without being told?**

*Probing an RNA Language Model with a human/chimpanzee variant that changes activity 4-fold*

Vincent Le Duigou — Albert School Madrid
NLP / Transformers course project
2026-05-19

---

## Slide 2 — Context (RNA as a "language")

- An RNA sequence = a 4-letter "sentence" (A, U, G, C).
- Its biological function depends on its **3D structure**, which depends on the **sequence**.
- Since 2024, the best RNA 3D structure predictors are **transformers**:
  - **RhoFold+** (*Nature Methods 2024*) — AlphaFold-style (Evoformer).
  - **DRfold2** (*PLOS Biology 2025*) — built on an **RNA Language Model** (RCLM).
- Classic NLP probing question: *"did the LM learn grammar, or only surface statistics?"*.
- Equivalent biology question: *"did the RNA LM learn mechanism, or only the average shape?"*

---

## Slide 3 — Our probe: CPEB3 human vs chimpanzee

**A unique sequence pair in the Stanford RNA 3D Folding benchmark (Kaggle 2025):**

| | Sequence (69 nt) | Position 30 | Activity |
|---|---|---|---|
| **R1107 (human)** | …UCGC**A**GCCC… | **A** | 1× |
| **R1108 (chimp)** | …UCGC**G**GCCC… | **G** | **4×** faster |

- **One single mutation, 4-fold activity difference.**
- Skilandat *et al.* (RNA 2016): the mutation disrupts a P1/P1.1 base-pairing in the P1 helix (residues 9 and 60).
- **Does the LM — which never saw kinetics or biochemistry — place the structural consequence at the right place?**

---

## Slide 4 — The result: DRfold2 recovers the P1 anchor

> Method: predict R1107 and R1108 with DRfold2, Kabsch-align, look at where the deviation is largest.

**Top-5 most-divergent residues:**
```
{9, 22, 24, 51, 60}
```
Residues **9** and **60** = **the two endpoints of the P1 helix** — exactly the region Skilandat 2016 identified biochemically.

The peak is NOT at the mutation site (position 30): it is 21–30 residues away, on long-range contacts.

> *Without functional supervision, the LM propagates the mutation toward the mechanistically relevant positions.*

---

## Slide 5 — The control: it's specific, not a prior

5 arbitrary single-nucleotide mutations elsewhere in the sequence:

| Mutation | Top-5 divergent | P1 anchor (pos 9 or 60)? |
|---|---|---|
| pos **30** (real, biological) | {**9**, 22, 24, 51, **60**} | **✓** |
| pos 5 | {1, 2, 3, 22, 23} | ✗ |
| pos 20 | {23, 24, 26, 48, 50} | ✗ |
| pos 41 | {48, 49, 50, 51, 52} | ✗ |
| pos 55 | {22, 23, 24, 25, 51} | ✗ |
| pos 64 | {50, 64, 65, 66, 67} | ✗ |

**0 / 5 arbitrary mutations** put the P1 anchor in the top-5. The anchor is **specific** to the biological mutation.

---

## Slide 6 — Why it matters for NLP / transformers

**1. Probing methodology:**
Just as one asks BERT "which syntactic features did you learn without being told?", we asked DRfold2 "which mechanisms did you learn without being told?".

**2. Emergent knowledge in a self-supervised LM:**
DRfold2 is trained on **raw, unlabeled RNA sequences**. It implicitly learns pairing rules and structural motifs, and — as we show — at least one *mechanistic signature* (the cascade toward P1).

**3. Architecture matters:**
RhoFold+ (Evoformer, designed for MSA input) in single-sequence mode **ignores the input**: same top-5 regardless of which mutation we apply. → Pure LMs and MSA-style transformers do **not respond the same way** to small perturbations.

**4. NLP analogy:**
- DRfold2 = a BERT that correctly changes its prediction when you swap "the *cat* sat" → "the *cats* sat".
- RhoFold+ without MSA = a model that outputs the same sentence regardless of input.

---

## Slide 7 — Conclusion and outlook

**What we showed:**
- An **RNA language model** recovers, without functional supervision, the mechanistic anchor of a human/chimpanzee variant pair documented biochemically.
- The signal is **specific** (0/5 control mutations reproduce it).
- Not all RNA LMs respond identically — there is a **differential calibration** angle to explore.

**Why it matters beyond biology:**
The micro-perturbation probing protocol applies to **any sequence language model** (NLP or bio). A single token swap reveals whether the model "understands" or "averages".

**Total experimental cost:** ~30 min on an RTX 5090 (~€0.30).
**Code and results:** https://github.com/Vince-Anduril/rna3d-calibration

> *"One letter, four-fold function: a transformer probe of mechanistic knowledge in RNA language models."*

**Thank you.**
