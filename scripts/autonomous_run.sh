#!/usr/bin/env bash
# Autonomous 4h pipeline: chains 3 jobs, commits/pushes after each.
# Failures are logged and don't block the next job.
set +e
cd /workspace/rna3d
source venv/bin/activate

GLOBAL_LOG=/workspace/rna3d/autonomous_run.log
echo "===== Autonomous run started $(date -u +%FT%TZ) =====" | tee -a $GLOBAL_LOG

commit_and_push() {
  local msg="$1"
  git add .research/ models/ code/ scripts/ STATUS.md 2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$msg" 2>&1 | tee -a $GLOBAL_LOG
    git push origin main 2>&1 | tail -3 | tee -a $GLOBAL_LOG
  else
    echo "(nothing to commit for $msg)" | tee -a $GLOBAL_LOG
  fi
}

write_status() {
  local phase="$1"
  cat > /workspace/rna3d/STATUS.md <<EOF
# Autonomous Run Status

Last updated: $(date -u +%FT%TZ)
Current phase: $phase

See \`autonomous_run.log\` for verbose details.

## Phase history
$(grep -E "^---" $GLOBAL_LOG 2>/dev/null | tail -20)
EOF
}

# --- Job 01b: v2 regex categorization ---
echo "--- $(date -u +%FT%TZ) starting Job 01b: v2 regex ---" | tee -a $GLOBAL_LOG
write_status "Job 01b: v2 regex categorization"
python /workspace/rna3d/code/jobs/job_01b_v2_regex.py 2>&1 | tee -a $GLOBAL_LOG
commit_and_push "feat(01b): v2 regex categorization (RNase P->ribozyme, phages->viral_rna)"

# --- Job 02: install models ---
echo "--- $(date -u +%FT%TZ) starting Job 02: install models ---" | tee -a $GLOBAL_LOG
write_status "Job 02: installing RibonanzaNet + RhoFold"
bash /workspace/rna3d/code/jobs/job_02_install_models.sh
commit_and_push "chore(02): RibonanzaNet/RhoFold install attempt + INSTALL.md"

# --- Job 03: inference on smallest 50 sequences (only if RibonanzaNet imports) ---
echo "--- $(date -u +%FT%TZ) starting Job 03: inference small ---" | tee -a $GLOBAL_LOG
write_status "Job 03: RibonanzaNet inference on smallest 50 sequences"
python /workspace/rna3d/code/jobs/job_03_inference_small.py 2>&1 | tee -a $GLOBAL_LOG
commit_and_push "feat(03): RibonanzaNet inference on smallest 50 Stanford sequences"

# --- Final status ---
echo "===== Autonomous run finished $(date -u +%FT%TZ) =====" | tee -a $GLOBAL_LOG
write_status "DONE — see commits since launch"
commit_and_push "chore: autonomous run complete — final STATUS.md update"
