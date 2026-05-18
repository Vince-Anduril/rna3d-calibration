# 20 — Hypotheses (living list)

> Owner: rna-scientist. Each hypothesis is a falsifiable claim with explicit
> success and failure conditions. Picks up opportunities from
> `10_literature/novel_angles.md` or from focus-family selection in Sprint 0.

Entry format:

> ## H-NNN [STATUS] — Hypothesis title
> **From:** rna-scientist, YYYY-MM-DD, picked up from novel_angles.md entry "..."
> **Claim:** the falsifiable statement.
> **Metric:** how we measure it.
> **Pass condition:** quantitative threshold for "supported".
> **Fail condition:** quantitative threshold for "falsified".
> **Protocol:** path to `20_hypotheses/protocols/H-NNN.md`.

STATUS: PROPOSED | ACTIVE | SUPPORTED | FALSIFIED | ABANDONED

## H-001 [PROPOSED] — Focus-family inter-model error variance is concentrated, not diffuse
**From:** rna-scientist, 2026-05-18, placeholder seeded from CONTEXT.md thesis pivot (focus family F to be selected end of Week 1 from data signals).
**Claim:** Once we select focus family F at end of Week 1, the inter-model error variance on F (across RibonanzaNet, RhoFold+, and a third comparator — Vfold or Boltz-1) will be at least 2x the dataset-wide median inter-model error variance.
**Metric:** Inter-model variance of per-residue RMSD (over the 3 models, computed per residue then averaged per sequence) on sequences belonging to family F, compared to the median of the same per-sequence variance taken over all Stanford test sequences (n=62) and/or validation sequences (n=62).
**Pass condition:** median(per-seq inter-model variance | seq in F) >= 2.0 * median(per-seq inter-model variance | all seqs).
**Fail condition:** median(per-seq inter-model variance | seq in F) < 1.2 * median(per-seq inter-model variance | all seqs).
**Inconclusive band:** ratio in [1.2, 2.0) -> revisit choice of F or comparator set.
**Protocol:** `20_hypotheses/protocols/H-001.md` (to be written once F is fixed; will pin model versions, RMSD computation, family-membership definition derived from Stanford `description` clustering per data/MANIFEST.md, and bootstrap CI procedure for the small-n setting).
**Notes:** This is a placeholder hypothesis serving as a worked example of the entry format. The 2x threshold is a working-hypothesis effect size and may be revised once Sprint 0 Week 1 yields the empirical variance distribution; the revision (if any) must happen BEFORE F is chosen, to avoid post-hoc threshold tuning. Concretely, F will be one of the categories produced by the description-based categorization scheme defined in [ANSWERED] Q-2026-05-18-002 (riboswitch, ribozyme, pseudoknot, viral_rna, nmr_solution_motif, tRNA, loop_motif, ribosome_subunit) — restricted to categories that pass the >=10-sequence viability threshold across train+val+test as reported by rna-code (Q-2026-05-18-004).

