# 01_description_categorization

Applies the 11-class priority-ordered regex categorization scheme (v1, scheme
source: [ANSWERED] Q-2026-05-18-002 in `questions_for_scientist.md`) to the
three Stanford Kaggle sequence CSVs (`train_sequences.csv`,
`validation_sequences.csv`, `test_sequences.csv`) on the `description` column,
case-insensitive, first-match-wins. The `nmr_solution_motif` rule additionally
requires `len(sequence) < 60`. Outputs a per-sequence assignment CSV, a
grouped bar chart of counts per category x split, an `uncategorized.csv` for
manual review, and a `SUMMARY.md` flagging which categories pass the >=10
viability threshold (counted across train+val+test combined).

Re-run:

```
source /workspace/rna3d/venv/bin/activate
cd /workspace/rna3d/.research/30_experiments/runs/01_description_categorization
python script.py
```

Outputs land in `results/` (committed alongside the script).
