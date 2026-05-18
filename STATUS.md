# Pod run status — 2026-05-18T13:24:32Z

- Training output: `.research/30_experiments/runs/GPU_training_stage2/`
- Eval output: `.research/30_experiments/runs/GPU_training_eval/`
- Final checkpoint: `/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage1_final.pt`
- Full log: `pod_run.log`

## Last 20 lines of log

```
    raise ValueError(
ValueError: num_samples should be a positive integer value, but got num_samples=0
--- starting eval on /workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage1_final.pt ---
/usr/local/lib/python3.12/dist-packages/torch/nn/modules/transformer.py:392: UserWarning: enable_nested_tensor is True, but self.use_nested_tensor is False because encoder_layer.norm_first was True
  warnings.warn(
Device: cuda
Loaded checkpoint: step=151668, params=4.92M
Traceback (most recent call last):
  File "/workspace/rna3d/code/jobs_gpu/eval.py", line 206, in <module>
    main()
  File "/workspace/rna3d/code/jobs_gpu/eval.py", line 116, in main
    for cat in ["all"] + sorted(per_res["category"].dropna().unique().tolist()):
                                ~~~~~~~^^^^^^^^^^^^
  File "/workspace/rna3d/venv/lib/python3.12/site-packages/pandas/core/frame.py", line 4113, in __getitem__
    indexer = self.columns.get_loc(key)
              ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/workspace/rna3d/venv/lib/python3.12/site-packages/pandas/core/indexes/range.py", line 417, in get_loc
    raise KeyError(key)
KeyError: 'category'
--- final commit ---
```
