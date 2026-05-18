# Run 01 — Description categorization — SUMMARY

Scheme: 11 priority-ordered regex rules (v1, from Q-2026-05-18-002).
Splits processed: train (844), val (12), test (12).
Total rows: 868.
Viability threshold: >= 10 across train+val+test combined.

## Per-category counts

| category | train | val | test | total | viable |
|---|---:|---:|---:|---:|:--:|
| ribosome_subunit | 199 | 0 | 0 | 199 | OK |
| riboswitch | 17 | 2 | 2 | 21 | OK |
| ribozyme | 54 | 2 | 2 | 58 | OK |
| tRNA | 64 | 1 | 1 | 66 | OK |
| viral_rna | 75 | 1 | 1 | 77 | OK |
| pseudoknot | 11 | 0 | 0 | 11 | OK |
| nmr_solution_motif | 131 | 0 | 0 | 131 | OK |
| loop_motif | 50 | 1 | 1 | 52 | OK |
| complex_with_protein | 39 | 2 | 2 | 43 | OK |
| synthetic_designed | 2 | 1 | 1 | 4 | NO |
| other | 202 | 2 | 2 | 206 | OK |

**Viable categories (>= 10 total, excluding `other`):** ribosome_subunit, riboswitch, ribozyme, tRNA, viral_rna, pseudoknot, nmr_solution_motif, loop_motif, complex_with_protein.

## Confidence distribution

- high: 388
- medium: 274
- low: 206

## Ambiguity / edge cases

- 274 rows had >=2 rules match before priority resolution (medium confidence). Inspect by filtering `category_assignment.csv` on `confidence == 'medium'`.
- 206 rows fell back to `other` (low confidence). See `uncategorized.csv`.
- Sample uncategorized descriptions:
    - 'P1 HELIX NUCLEIC ACIDS (DNA/RNA) RIBONUCLEIC ACID'
    - 'RNA BACTERIOPHAGE MS2 COAT PROTEIN/RNA COMPLEX'
    - 'MS2 PROTEIN CAPSID/RNA COMPLEX'
    - 'MS2 PROTEIN CAPSID/RNA COMPLEX'
    - 'Bacteriophage Lambda N-protein-NutboxB-RNA Complex'

## Notes

- `nmr_solution_motif` requires both (`nmr` or `solution`) AND `len(sequence) < 60`; without the length gate that category absorbs most cryo-EM ribosome assemblies that mention `solution`.
- Priority order (ribosome -> riboswitch -> ribozyme -> tRNA -> viral -> pseudoknot -> nmr_solution_motif -> loop -> complex_with_protein -> synthetic_designed -> other) is taken verbatim from the scientist's Q-002 answer.
- All regexes case-insensitive. DOTALL is enabled so multi-line descriptions (validation/test rows wrap lines inside quotes) are matched correctly.
