# GPU training — evaluation (EMPTY)

Checkpoint: `/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage1_final.pt`

**No residues had valid ground-truth coordinates.** Usually means the
PDB file lookup failed (filenames/chain IDs don't match). Inspect
`per_sequence.csv` for `n_valid==0` rows and check `data.py`.

- per_sequence rows: 48
- per_residue rows: 0
- mean n_valid per sequence: 0.0
