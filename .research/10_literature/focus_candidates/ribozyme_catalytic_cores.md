# Ribozyme catalytic cores — reference for active-site analysis

> Owner: rna-research (extended). Reference notes from training knowledge.
> Use for: identifying catalytic residues in Stanford PDB ground-truth → enables
> the "sub-RMSD on active site only" metric proposed in `H-001` / `FOCUS_DECISION.md`.
>
> **TO VERIFY** during paper writing: every numerical claim and conserved-position
> assignment must be cross-checked against the cited primary source before
> appearing in the manuscript. Treat this as a working draft, not a peer-reviewed
> consensus.

---

## Background

Ribozymes are catalytic RNAs. They can be grouped into:

- **Small self-cleaving ribozymes** — cleave a single phosphodiester bond,
  typically 30-200 nt: hammerhead, hairpin, HDV, glmS, twister, pistol, hatchet
- **Large catalytic ribozymes** — multi-domain, often >200 nt: group I intron,
  group II intron, RNase P, ribosome (rRNA)
- **Ribosomal RNA (rRNA) PTC** — peptidyl transferase center, part of large
  ribosomal subunit; technically a ribozyme

The **active site** in most small ribozymes is concentrated on a small number of
specific nucleotides whose 3D positioning is essential for catalysis. Mis-modeling
those positions by more than ~2 Å typically abolishes any meaningful claim about
catalytic geometry.

## Per-class active-site landmarks

Each entry below names the ribozyme class, gives a brief description, and lists
the residues conventionally considered to form the catalytic core. Coordinate
identifiers refer to the canonical numbering used in the cited reference
structures; in our Stanford data the actual residue indices depend on the PDB
entry.

### Hammerhead ribozyme

- **Function:** cleaves substrate RNA by 2'-OH attack on adjacent phosphate.
- **Catalytic core (minimal):** ~13 nt across two stem loops + a junction.
- **Key residues:** the central conserved CUGANGA box (positions ~3-9 in the
  classic numbering), plus the cleavage-site invariant nucleotide downstream.
  **TO VERIFY** exact position table from Martick & Scott 2006 or
  Lilley 2011 review.
- **Mechanistic stars:** G8 and G12 are widely implicated as general acid /
  general base candidates (depending on the protonation model).
- **Notes:** ~50 nt minimal active hammerheads have been crystallized; full
  natural hammerheads are larger with auxiliary tertiary contacts.

### Hairpin ribozyme

- **Function:** reversible cleavage / ligation.
- **Catalytic core:** internal-loop A and internal-loop B docked together;
  the catalytic site sits at the interface.
- **Key residues:** G8 (loop A) and A38 (loop B) are conventionally the
  catalytic stars. **TO VERIFY** numbering against Rupert et al. 2002
  (PDB 1HP6 for the loop-A loop-B complex).
- **Notes:** the docking of loop A onto loop B is the rate-limiting
  conformational event in many models — meaning the *interface geometry*
  matters as much as any single residue position.

### Hepatitis Delta Virus (HDV) ribozyme

- **Function:** self-cleavage during HDV genome replication.
- **Catalytic core:** double-pseudoknot fold (nested pseudoknots P1.1, P2).
- **Key residues:** C75 (or C76 in some numbering) is the most-discussed
  catalytic residue; M2+ ions in the active site are also essential.
- **Notes:** HDV is one of the classic test cases for pseudoknot-aware
  prediction — a model that misses the pseudoknot topology has zero chance
  of placing C75 correctly.

### glmS ribozyme

- **Function:** self-cleavage triggered by glucosamine-6-phosphate cofactor binding.
- **Catalytic core:** double pseudoknot with a defined cofactor pocket.
- **Key residues:** G33 and A-1 (relative to cleavage site) are commonly
  cited as catalytic. The bound cofactor itself is required.
- **Notes:** one of the few catalytic riboswitches; the cofactor pocket is
  geometrically demanding to predict without cofactor context (we don't
  predict cofactor positions in our pipeline — explicit limitation).

### Twister, pistol, hatchet (newer small ribozymes)

- **Function:** small self-cleaving ribozymes discovered ~2014-2017.
- **Catalytic cores:** characterized by Liu/Breaker labs; each has 3-5
  conserved residues at the cleavage site.
- **Notes:** under-represented in older training data. Likely candidates for
  our paper's "models do/don't generalize to newer classes" angle.

### Group I intron

- **Function:** self-splicing intron; uses an exogenous guanosine to initiate.
- **Catalytic core:** multi-helix active site; the catalytic core is built from
  ~6 paired domains (P1-P9) that fold into a precise geometry.
- **Key residues:** the guanosine-binding pocket residues + the 5' splice-site
  recognition region. Specific atom-level interactions are complex.
- **Notes:** these are LARGE molecules (~200-400 nt). Several have val/test
  presence in the Stanford set (TBD count).

### Group II intron

- **Function:** self-splicing intron; uses lariat intermediate (transesterification).
- **Catalytic core:** "domain V" hairpin is the most conserved catalytic feature.
- **Key residues:** the AGC triad in DV plus a bulged adenosine.

### RNase P RNA (M1)

- **Function:** processes pre-tRNA into mature tRNA (5' leader removal).
- **Catalytic core:** S-domain + C-domain joined by a junction; the C-domain
  contains the catalytic residues.
- **Key residues:** U69 region (specifics depend on the bacterial vs
  archaeal vs eukaryotic variant). **TO VERIFY** from Marvin & Engelke 2009
  or Esakova/Krasilnikov 2010 reviews.
- **Notes:** **rescued from `other` to `ribozyme` by the v2 regex** —
  RNase P entries may dominate our ribozyme set numerically.

### Ribosomal peptidyl transferase center (PTC)

- **Function:** peptide bond formation in protein synthesis.
- **Catalytic core:** within 23S rRNA (bacteria) / 28S rRNA (eukaryotes).
- **Key residues:** A2451 / A2602 (E. coli numbering) in the A-site / P-site
  vicinity. **TO VERIFY** from a current ribosome structural review.
- **Notes:** the PTC is inside large rRNA — most of our `ribosome_subunit`
  category is large rRNA structures. We **explicitly do NOT focus** on the
  PTC because (a) `ribosome_subunit` was filtered out as over-studied, and
  (b) ribosomal entries are typically too large for our 256-token model.

---

## Implication for our analysis pipeline

When we evaluate predictions on the 4 val+test ribozyme sequences (and the
~58 train ribozyme sequences), we can compute two metrics in parallel:

1. **Global RMSD** — over all residues (what `kabsch_rmsd` in `eval.py` already
   computes).
2. **Active-site sub-RMSD** — restricted to the catalytic core residues. To
   do this we need a per-sequence mapping (target_id → list of catalytic
   residue indices in the Stanford sequence numbering). This mapping is
   **not in the Stanford data**; we have to construct it from the PDB
   metadata + this reference list (or from external annotation sources like
   RNA-Puzzles or PDB ligand-binding-site annotations).

**Construction strategy (post-pod):**

1. For each ribozyme PDB ID in our val+test, look up the PDB entry on RCSB
   (https://www.rcsb.org) to get the canonical name + reference paper.
2. Map the name to one of the classes above using the conservative regex
   in `code/jobs_gpu/data.py`.
3. From the reference paper / Lilley 2011 / class-specific reviews,
   identify the catalytic residues in the canonical numbering.
4. Use sequence alignment (e.g., ViennaRNA / biotite) to map canonical
   numbering → Stanford-sequence numbering for each entry.
5. Cache the mapping in `.research/30_experiments/runs/ribozyme_active_sites.csv`.

This is **not in the first pod run** — it's a Week-2 deliverable. The first pod
run produces global RMSD + calibration; the sub-RMSD comes after.

---

## What this file is NOT

- A peer-reviewed consensus. Numbering and "key residues" assignments above
  are working notes from training knowledge — every claim must be re-verified
  against the cited primary source before any paper claim is made.
- A complete reference. Many small ribozyme classes and detailed catalytic
  mechanism debates are omitted for brevity.
- A reason to skip running the analysis. Even before we map active sites,
  global RMSD + per-residue calibration on ribozymes is publishable
  content — the sub-RMSD is the *extra* lens this file enables.
