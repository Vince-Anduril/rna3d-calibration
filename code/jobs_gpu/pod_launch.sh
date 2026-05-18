#!/usr/bin/env bash
# Autonomous pod launcher — runs the full GPU training + eval pipeline.
# Designed to be invoked once on the pod, then run unattended for ~3 hours.
# Auto-commits + auto-pushes incrementally so the GitHub mirror always has the latest state.
set +e  # don't abort on individual failures
cd /workspace/rna3d

# Make sure git knows our identity (may have been wiped on a fresh pod)
git config user.email "vleduigou@c2tclinic.com" 2>/dev/null
git config user.name "Vincent Le Duigou" 2>/dev/null

LOG=/workspace/rna3d/pod_run.log
echo "===== POD AUTONOMOUS RUN started $(date -u +%FT%TZ) =====" | tee -a $LOG
echo "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)" | tee -a $LOG

# 1) Pull latest code (in case the user pushed updates from Mac)
echo "--- git pull ---" | tee -a $LOG
git pull --ff-only origin main 2>&1 | tail -3 | tee -a $LOG

# 2) Ensure venv is sound (idempotent)
echo "--- venv check ---" | tee -a $LOG
if [ ! -d /workspace/rna3d/venv ]; then
  python3 -m venv /workspace/rna3d/venv --system-site-packages
fi
source /workspace/rna3d/venv/bin/activate
pip install -q --upgrade pip 2>&1 | tail -1 | tee -a $LOG
# Install any missing deps in the venv (idempotent)
pip install -q -r /workspace/rna3d/requirements.txt 2>&1 | tail -3 | tee -a $LOG

# 3) Verify GPU
echo "--- torch GPU check ---" | tee -a $LOG
python -c "import torch; assert torch.cuda.is_available(); print('CUDA OK:', torch.cuda.get_device_name(0))" 2>&1 | tee -a $LOG

# 4) Place the authored code under code/ (it's the Mac authoring tree; the pod-side
#    expectation is that `code/jobs_gpu/` already has these scripts, or the user
#    SCPed them up before launch; if not present, abort).
if [ ! -f /workspace/rna3d/code/jobs_gpu/train.py ]; then
  echo "ERROR: /workspace/rna3d/code/jobs_gpu/train.py missing. Did you scp _authoring/code/* to /workspace/rna3d/code/jobs_gpu/?" | tee -a $LOG
  exit 1
fi

# 5) Run training (default 2.8h)
TRAIN_OUT=/workspace/rna3d/.research/30_experiments/runs/GPU_training_stage2
mkdir -p $TRAIN_OUT
echo "--- starting training (target ${HOURS:-2.8}h) ---" | tee -a $LOG
python /workspace/rna3d/code/jobs_gpu/train.py \
  --hours "${HOURS:-2.8}" \
  --out-dir "$TRAIN_OUT" \
  --repo-root /workspace/rna3d \
  2>&1 | tee -a $LOG

# 6) Run eval on the final checkpoint
LATEST_CKPT=$(ls -t $TRAIN_OUT/stage2_*.pt 2>/dev/null | head -1)
if [ -z "$LATEST_CKPT" ]; then
  LATEST_CKPT=$(ls -t $TRAIN_OUT/stage1_*.pt 2>/dev/null | head -1)
fi
EVAL_OUT=/workspace/rna3d/.research/30_experiments/runs/GPU_training_eval
mkdir -p $EVAL_OUT
echo "--- starting eval on $LATEST_CKPT ---" | tee -a $LOG
if [ -n "$LATEST_CKPT" ]; then
  python /workspace/rna3d/code/jobs_gpu/eval.py \
    --ckpt "$LATEST_CKPT" \
    --out-dir "$EVAL_OUT" \
    2>&1 | tee -a $LOG
else
  echo "no checkpoint found; skipping eval" | tee -a $LOG
fi

# 7) Final commit + push (training script also commits incrementally)
echo "--- final commit ---" | tee -a $LOG
git add .research/30_experiments/runs/GPU_training_stage2/ \
        .research/30_experiments/runs/GPU_training_eval/ \
        STATUS.md 2>/dev/null
cat > /workspace/rna3d/STATUS.md <<EOF
# Pod run status — $(date -u +%FT%TZ)

- Training output: \`.research/30_experiments/runs/GPU_training_stage2/\`
- Eval output: \`.research/30_experiments/runs/GPU_training_eval/\`
- Final checkpoint: \`$LATEST_CKPT\`
- Full log: \`pod_run.log\`

## Last 20 lines of log

\`\`\`
$(tail -20 $LOG)
\`\`\`
EOF
git add STATUS.md
git commit -m "feat(GPU_run): training + eval complete on RTX 5090 ($(date -u +%FT%TZ))" 2>&1 | tail -3 | tee -a $LOG
git push origin main 2>&1 | tail -3 | tee -a $LOG

echo "===== POD AUTONOMOUS RUN finished $(date -u +%FT%TZ) =====" | tee -a $LOG

# 8) AUTO-STOP THE POD to avoid burning money.
# Tries multiple paths in order of preference:
#   (a) runpodctl with cached config + RUNPOD_POD_ID env var
#   (b) runpodctl with hostname as pod id
#   (c) graceful container shutdown
echo "--- auto-stop attempt ---" | tee -a $LOG
POD_ID="${RUNPOD_POD_ID:-$(hostname)}"
if command -v runpodctl >/dev/null 2>&1; then
  echo "trying: runpodctl stop pod $POD_ID" | tee -a $LOG
  runpodctl stop pod "$POD_ID" 2>&1 | tee -a $LOG
  # Give the API ~30s; if we are still alive, fall through.
  sleep 30
fi
echo "fallback: poweroff" | tee -a $LOG
poweroff 2>&1 | tee -a $LOG || shutdown -h now 2>&1 | tee -a $LOG || true
