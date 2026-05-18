# Focus pick — Ribozyme

**Decision date:** 2026-05-18 (autonomous, per user grant of direction authority)
**Decided by:** orchestrator (Claude main session)
**To be overridden by user before pod launch if desired**

## Pick: RIBOZYME

From the FOCUS_SHORTLIST.md candidates (ribozyme STRONG, riboswitch STRONG, viral_rna MEDIUM), the orchestrator selects **ribozyme**.

## Empirical signal (from run 01_description_categorization, v2 regex applied)

- Train: 58 (or higher with v2 RNase P rescue, expected ~70)
- Val: 2
- Test: 2
- Val+test: 4 — sufficient for evaluation, comparable to riboswitch

## Justification (4 criteria from 00_thesis.md)

1. **Sharpest falsifiable metric.** Catalytic ribozymes have a **small, well-defined active site** (3-7 residues per class — see Lilley 2011 catalytic core review). A sub-RMSD restricted to those residues is the cleanest small-n metric in the entire shortlist. By contrast, riboswitch ON/OFF requires having *both* conformers in the data, which we cannot guarantee.

2. **Biological centrality.** Catalysis IS the defining function of a ribozyme. A finding of the form "the models reproduce overall fold but miss catalytic geometry by X Å with high-confidence predictions" is biomechanically meaningful and immediately interpretable to the broader community.

3. **Bounded scope (3-month feasibility).** Ribozyme classes are few and well-characterized: hammerhead, hairpin, group I intron, group II intron, RNase P, glmS, twister, pistol, hatchet, HDV. We can stratify by class within the family. Easier to write a clean paper than the heterogeneous viral_rna class.

4. **Visual story.** Active-site overlay (predicted vs ground-truth) on a single ribozyme makes a single compelling Figure 1.

## Secondary lens (kept regardless)

Pseudoknot re-detection across all categories (Week-2 work) remains queued — independent of the ribozyme focus.

## Override

If the user prefers riboswitch (conformational duality angle) or viral_rna (data
volume angle), update `.research/00_thesis.md` and this file with the new pick
before launching the pod. Default pick proceeds as ribozyme until overridden.
