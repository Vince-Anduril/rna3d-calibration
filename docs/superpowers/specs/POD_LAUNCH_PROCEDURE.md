# Pod launch procedure (when pod restarts)

Concise checklist. The pod will land at a NEW IP/port (RunPod re-assigns on every restart).

## Pre-flight (you, 30 seconds)

1. Start the pod from the RunPod console.
2. Copy the "SSH over exposed TCP" line. It will look like:
   `ssh root@<NEW_IP> -p <NEW_PORT> -i ~/.ssh/id_ed25519`
3. Paste that line in chat — Claude does the rest.

## What Claude does (no input needed)

1. Verify TCP reachable: `nc -z <NEW_IP> <NEW_PORT>`.
2. If `/workspace/rna3d/` survived the stop (likely — it's persistent network mount):
   ```
   ssh -p <NEW_PORT> root@<NEW_IP> 'cd /workspace/rna3d && git pull --ff-only origin main'
   ```
   This pulls `code/jobs_gpu/` (the GPU training pipeline that was authored on Mac).
3. If `/workspace/rna3d/` was wiped (rare but possible on hardware change), re-bootstrap:
   - Clone repo: `git clone git@github.com:Vince-Anduril/rna3d-calibration.git /workspace/rna3d`
   - Re-download dataset (~10 min): `kaggle competitions download -c stanford-rna-3d-folding`
   - Re-extract (~2 min): `unzip ...`
   - Re-create venv: `python3 -m venv venv --system-site-packages && pip install -r requirements.txt`
4. Launch the autonomous training:
   ```
   ssh -p <NEW_PORT> root@<NEW_IP> '
     cd /workspace/rna3d
     nohup bash code/jobs_gpu/pod_launch.sh > pod_run.log 2>&1 &
     disown
     echo "launched, PID=$!"
   '
   ```
5. Verify it's running and GPU is starting:
   ```
   ssh -p <NEW_PORT> root@<NEW_IP> 'pgrep -af pod_launch; nvidia-smi'
   ```

## Expected timeline

| Phase | Duration | GPU util |
|---|---|---|
| Bootstrap recovery (only if /workspace was wiped) | ~15 min | 0% |
| Stage 1: MLM pretraining | ~30 min | 70-90% |
| Stage 2: structure supervision | ~2.3 h | 80-95% |
| Eval | ~5 min | 30-50% |
| **Total** | **~3 h** (or +15 min if recovery needed) | |

## Auto-commit cadence

The pipeline pushes to GitHub:
- Every 15 min of training (loss curves, log)
- At the end of Stage 1 (final.pt name)
- At the end of Stage 2 (checkpoint names + final losses)
- At the end of eval (per_category.csv, calibration.csv, SUMMARY.md, reliability_diagrams.png)

You can monitor live at: https://github.com/Vince-Anduril/rna3d-calibration/commits/main

## When it finishes

Claude (next session) reads:
- `.research/30_experiments/runs/GPU_training_eval/SUMMARY.md` — the headline numbers
- `calibration_per_category.csv` — per-class ECE with bootstrap CIs
- `reliability_diagrams.png` — the key figure
- Determines whether H-001 (ribozyme/global ECE ratio ≥ 2x) is SUPPORTED, FALSIFIED, or INCONCLUSIVE
- If SUPPORTED → drafts Findings entry + paper figure
- If FALSIFIED → pivots to riboswitch backup or refines model
- If INCONCLUSIVE → adjusts threshold / runs longer / adds bootstrap iterations

## Cost estimate

RTX 5090 spot pricing on RunPod is typically around 0.50-0.80€/h. A single 3h run = ~1.50-2.50€.

## Failure recovery

- **SSH fails:** check `chmod 600 ~/.ssh/authorized_keys` on the pod (sshd's StrictModes rejects 0700, which is the RunPod default).
- **Job dies mid-training:** the auto-commit cadence means GitHub has the last checkpoint name + losses up to 15 min before the failure. The pod's `/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2/stage2_step_NNNNNN.pt` files are still on the pod's disk. Resume training from the latest checkpoint by adapting `train.py` (not currently wired; add `--resume` flag if needed).
- **Out of memory:** the training script logs OOMs and skips problem sequences. If too many sequences fail, reduce `--batch-struct` from 8 to 4 or `--max-len` from 256 to 192.
- **GitHub push fails:** the auto-push is best-effort; if the deploy key was rotated, fix and re-push from the next commit.
