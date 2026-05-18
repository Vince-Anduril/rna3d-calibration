# 00_pod_smoke

End-to-end pod execution smoke test. Confirms the RTX 5090 + torch 2.8 (cu128) venv works, that `train_sequences.csv` loads cleanly, and produces one tiny project-meaningful signal: the 15 most common single words (lowercase, punctuation-stripped) in the `description` column. This serves as a first, very coarse view of which biological classes dominate the training set (ribosome, NMR/cryo-EM, loop, complex...) before any clustering work.

## Re-run

```bash
cd /workspace/rna3d/.research/30_experiments/runs/00_pod_smoke
source /workspace/rna3d/venv/bin/activate && python script.py > results/output.txt
```

## Expected

- `torch.cuda.is_available = True`, device name `NVIDIA GeForce RTX 5090`.
- Train set: 844 rows, 844 unique target_ids, mean sequence length ~162 nt (min 3 / max 4298).
- Word counts dominated by structural-biology vocabulary: `structure`, `rna`, `complex`, `solution`, `cryo`, `em`, `loop`, `ribosome`, `nmr`.
