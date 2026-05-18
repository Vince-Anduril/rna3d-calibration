# Ready for Pod — GPU Training Run Briefing

**Date:** 2026-05-18
**Author:** orchestrator (Claude main session)
**For:** Vincent (user)

---

## What you'll get out of this run

A trained custom RNA structure model + a stratified calibration evaluation,
focused on the **ribozyme** family (auto-picked from the Sprint 0 Week 1
shortlist — see `_authoring/docs/FOCUS_DECISION.md`).

Concrete artifacts committed to GitHub after the run:

- `.research/30_experiments/runs/GPU_training_stage2/`
  - `stage1_final.pt` — model after MLM pretraining
  - `stage2_step_NNNNNN.pt` — incremental structure-supervision checkpoints (every 10 min)
  - `stage2_final.pt` — final weights
  - `stage2_losses.csv` — full loss curves
  - `train.log` — verbose training log
- `.research/30_experiments/runs/GPU_training_eval/`
  - `per_sequence.csv` — Kabsch-aligned RMSD + mean confidence per val/test sequence
  - `per_residue.csv` — per-residue (confidence, error) — the calibration raw data
  - `per_category.csv` — aggregated by category (ribozyme highlighted)
  - `reliability_diagrams.png` — predicted-confidence vs realized-error, global and ribozyme
  - `SUMMARY.md` — paper-ready narrative

## What's running

1. **Stage 1 (~30 min): Masked-LM pretraining** on the ≈5k unique v2-train sequences. Learns RNA token semantics. Saturates GPU.
2. **Stage 2 (~2.3 h): Structure supervision** — supervised on the ≈800 train sequences that have PDB ground-truth. Joint loss: Kabsch-aligned coord MSE + pair-distance auxiliary + per-residue calibration loss (confidence → exp(-error/5Å)).
3. **Eval (~5 min)** on val + test (62+62 sequences). Computes per-sequence RMSD + reliability diagrams per category, with ribozyme deep-dive.

**Total runtime:** ~3 hours. **GPU utilization:** expected 80–95% throughout (transformer with batch=8 sequences × 256 tokens on RTX 5090).

## Model details

- 6-layer transformer encoder, d_model=256, 8 heads, FF=1024
- Coord head: 3-layer MLP → (x, y, z)
- Confidence head: 3-layer MLP → sigmoid scalar
- **~5M parameters** (deliberately small — this is a CONTROLLED predictor for calibration analysis, not a SOTA attempt)

## Why this design

We can't reliably install RhoFold / RibonanzaNet-3D / Boltz-1 in scoped pod time.
Instead, we train our own small predictor we fully understand — this gives us:

1. **Explicit calibration channel.** Confidence head is trained from scratch with a
   well-defined target (exp(-error/5Å)). No ambiguity about "what does this score mean".
2. **Ribozyme calibration analysis** with a model whose biases we can fully characterize.
3. **A defensible paper claim.** "We trained a small, controllable predictor to study
   prediction-error calibration on ribozyme structures from the Stanford 3D Folding
   benchmark." Humble, novel, focused.
4. **Reproducible.** ≈5M params + 868 sequences fits in a single 3-hour pod session.

## To launch (your action)

```bash
# 1) Restart the pod — get the new IP and TCP port from the RunPod console.
# 2) Update memory with the new connection details (or tell Claude the new IP/port).
# 3) From the Mac, SCP the authoring code up to the pod:
scp -P <NEW_PORT> -r "/Users/leduigouvincent/Documents/cours madrid/Data/Kaggle Stanford ADN /_authoring/code/" \
    root@<NEW_IP>:/workspace/rna3d/code/jobs_gpu/
# 4) SSH in and launch:
ssh -p <NEW_PORT> root@<NEW_IP> '
  cd /workspace/rna3d
  git pull --ff-only origin main
  nohup bash code/jobs_gpu/pod_launch.sh > pod_run_nohup.log 2>&1 &
  disown
  echo "launched, PID=$!"
'
```

…or just tell Claude "pod is on at <NEW_IP>:<NEW_PORT>, launch" and it'll handle the SCP + ssh + nohup chain itself.

## What gets pushed to GitHub during the run

Auto-commit cadence (every 15 minutes of training) pushes:
- The loss curves (`stage2_losses.csv`)
- The latest checkpoint NAME (the .pt files themselves stay on pod — they're too big)
- The training log

Final commit (at end of run) pushes everything in
`GPU_training_stage2/` and `GPU_training_eval/`.

You can watch progress at:
https://github.com/Vince-Anduril/rna3d-calibration/commits/main

## After the run — analysis cadence

When you next sit with Claude:
1. `rna-scientist` reads `eval/SUMMARY.md` + `per_category.csv` + `per_residue.csv`
2. Decides whether the ribozyme calibration signal is publishable (e.g. ECE on ribozyme
   ≥2× global ECE → meets H-001 spirit threshold)
3. If yes → drafts the paper's Figure 1 + Table 1 from these results
4. If no → revises focus (riboswitch fallback) or model architecture

## Risks I've accepted

- **Custom model won't match SOTA accuracy.** That's fine — the paper claim is about
  calibration characterization, not raw accuracy ranking.
- **Some PDB files may fail to parse** (different formats). Job logs these; non-fatal.
- **Some sequences may be too long for max_len=256.** They get truncated. Logged.
- **The pod may die mid-training.** Checkpoints every 10 min + auto-push every 15 min
  means we recover at most 15 minutes of work.

## What's NOT in this run (deferred)

- RhoFold / Boltz-1 / other established 3D RNA predictors — would need a separate,
  scoped install session. Worth a Week-2 pod run.
- Pseudoknot detection across categories (Week-2 queued item).
- Multi-model ensemble for inter-model variance (needs at least 2 working 3D models).

## Files staged on Mac (will SCP at pod launch)

```
_authoring/code/
  model.py          # transformer + dual head (~5M params)
  loss.py           # MLM + Kabsch-aligned structure + confidence
  data.py           # Stanford CSV + PDB loader (+ SyntheticRNADataset for smoke)
  train.py          # 2-stage trainer with incremental checkpoint + commit
  eval.py           # per-category metrics + reliability diagrams
  pod_launch.sh     # autonomous launcher
_authoring/docs/
  FOCUS_DECISION.md # ribozyme pick + justification
  READY_FOR_POD.md  # this file
```

All smoke-tested on Mac CPU with synthetic data — forward, loss, and a minimal
training loop all run cleanly (4.92M params, MLM/structure/confidence losses all OK).
