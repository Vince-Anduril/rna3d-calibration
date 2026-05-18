# Pod run status — 2026-05-18T13:32:15Z

- Training output: `.research/30_experiments/runs/GPU_training_stage2/`
- Eval output: `.research/30_experiments/runs/GPU_training_eval/`
- Final checkpoint: `/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage1_final.pt`
- Full log: `pod_run.log`

## Last 20 lines of log

```
    main()
  File "/workspace/rna3d/code/jobs_gpu/train.py", line 299, in main
    stage2_structure(model, device, ds, out_dir, log_fp, deadline, batch_size=args.batch_struct)
  File "/workspace/rna3d/code/jobs_gpu/train.py", line 181, in stage2_structure
    struct = structure_loss(out["coords"].float(), coords, valid)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/workspace/rna3d/code/jobs_gpu/loss.py", line 64, in structure_loss
    aligned = _kabsch_align(pred_coords, true_coords, valid_mask)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/workspace/rna3d/code/jobs_gpu/loss.py", line 41, in _kabsch_align
    U, S, Vt = torch.linalg.svd(H)
               ^^^^^^^^^^^^^^^^^^^
NotImplementedError: "svd_cuda_gesvdjBatched" not implemented for 'BFloat16'
--- starting eval on /workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage1_final.pt ---
/usr/local/lib/python3.12/dist-packages/torch/nn/modules/transformer.py:392: UserWarning: enable_nested_tensor is True, but self.use_nested_tensor is False because encoder_layer.norm_first was True
  warnings.warn(
Device: cuda
Loaded checkpoint: step=151668, params=4.92M
Eval wrote empty SUMMARY (no per-residue ground truth). Check /workspace/rna3d/.research/30_experiments/runs/GPU_training_eval/SUMMARY.md
--- final commit ---
```
