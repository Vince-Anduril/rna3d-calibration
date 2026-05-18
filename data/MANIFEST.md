# Data Manifest

> All files in `data/kaggle_raw/` are from the Stanford RNA 3D Folding Kaggle
> competition (closed 2025-09-24), downloaded via `kaggle competitions download
> -c stanford-rna-3d-folding`. The raw zip is deleted post-extraction.

## File inventory (as of 2026-05-18 post-extraction)

| Path | Size | Rows / files | Description |
|---|---|---|---|
| `kaggle_raw/train_sequences.csv` | 3.0 MB | 19,087 rows | v1 train: target_id, sequence, temporal_cutoff, description, all_sequences |
| `kaggle_raw/train_sequences.v2.csv` | 54 MB | 313,838 rows | v2 train (newer release, larger) |
| `kaggle_raw/train_labels.csv` | 9.3 MB | 137,096 rows | v1 train labels (per-residue) |
| `kaggle_raw/train_labels.v2.csv` | 256 MB | 3,677,096 rows | v2 train labels |
| `kaggle_raw/validation_sequences.csv` | 10 KB | 63 rows | val seqs (62 + header) |
| `kaggle_raw/validation_labels.csv` | 2.4 MB | 2,516 rows | val labels |
| `kaggle_raw/test_sequences.csv` | 10 KB | 63 rows | test seqs (62 + header). **Very small.** |
| `kaggle_raw/sample_submission.csv` | 185 KB | 2,516 rows | Kaggle submission format |
| `kaggle_raw/MSA/` | 448 MB | 856 files | v1 MSAs as FASTA, one per target_id |
| `kaggle_raw/MSA_v2/` | 3.2 GB | 2,534 files | v2 MSAs (newer, larger) |
| `kaggle_raw/PDB_RNA/` | 57 GB | 8,672 files | RNA structure files (PDB format) — ground-truth source |
| `kaggle_winner_notebooks/adamlogman_randomness.ipynb` | 19 KB | — | Winner-tier notebook for reference |

## Sample rows

```
$ head -2 train_sequences.csv
target_id,sequence,temporal_cutoff,description,all_sequences
1SCL_A,GGGUGCUCAGUACGAGAGGAACCGCACCC,1995-01-26,"THE SARCIN-RICIN LOOP, A MODULAR RNA",">1SCL_1|..."
```

## Key observation for focus-candidate selection

The `description` column in `train_sequences.csv` contains free-text biological
labels (e.g. "THE SARCIN-RICIN LOOP, A MODULAR RNA"). These can be clustered or
keyword-matched into RNA-functional classes WITHOUT needing an external Rfam
mapping. Sprint 0 Week 1 should explore this:

1. Tokenize / cluster descriptions into approximate families.
2. Cross-reference with `target_id` prefixes (PDB IDs, e.g. "1SCL_A" → 1SCL chain A).
3. Identify families that are: (a) reasonably well-represented in train (≥ a
   few dozen sequences), (b) present in val or test, (c) not dominated by
   any single famous structure that's already been over-studied (e.g. avoid
   "tRNA" if it's 5 sequences all from one organism).

## Provenance and reproducibility

- Source: https://www.kaggle.com/competitions/stanford-rna-3d-folding
- Download command: `kaggle competitions download -c stanford-rna-3d-folding`
- Downloaded zip size: 13.4 GB (deleted after extraction)
- Extraction: `unzip -o stanford-rna-3d-folding.zip`
- Date downloaded on this pod: 2026-05-18
