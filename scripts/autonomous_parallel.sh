#!/usr/bin/env bash
# Launches 3 jobs in parallel (2 CPU-heavy on different data, 1 GPU-heavy).
# Each job logs to its own file. A watcher commits/pushes whenever a SUMMARY.md appears.
set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/autonomous_parallel.log
echo "===== Autonomous PARALLEL run started $(date -u +%FT%TZ) =====" | tee -a $LOG

# Launch all 3 in background
echo "--- Launching Job A (PDB parse, CPU heavy, 64 procs) ---" | tee -a $LOG
nohup python /workspace/rna3d/code/jobs/job_A_pdb_parse.py > /workspace/rna3d/.research/30_experiments/runs/A_pdb_parse/stdout.log 2>&1 &
A_PID=$!
echo "JOB_A_PID=$A_PID" | tee -a $LOG

echo "--- Launching Job B (RNA-FM embeddings, GPU heavy) ---" | tee -a $LOG
nohup python /workspace/rna3d/code/jobs/job_B_rnafm_embed.py > /workspace/rna3d/.research/30_experiments/runs/B_rnafm_embeddings/stdout.log 2>&1 &
B_PID=$!
echo "JOB_B_PID=$B_PID" | tee -a $LOG

echo "--- Launching Job C (MSA stats, CPU heavy, 32 procs) ---" | tee -a $LOG
nohup python /workspace/rna3d/code/jobs/job_C_msa_stats.py > /workspace/rna3d/.research/30_experiments/runs/C_msa_stats/stdout.log 2>&1 &
C_PID=$!
echo "JOB_C_PID=$C_PID" | tee -a $LOG

# Watcher: commits + pushes every 90s while any job is running; final commit when all done
commit_push() {
  local msg="$1"
  cd /workspace/rna3d
  git add .research/ 2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$msg" 2>&1 | tee -a $LOG
    git push origin main 2>&1 | tail -3 | tee -a $LOG
  fi
}

write_status() {
  cat > /workspace/rna3d/STATUS.md <<EOF
# Autonomous Parallel Run Status

Last updated: $(date -u +%FT%TZ)

## Jobs

| Job | PID | Status | Output |
|---|---|---|---|
| A — PDB parse (CPU) | $A_PID | $(kill -0 $A_PID 2>/dev/null && echo RUNNING || echo DONE) | .research/30_experiments/runs/A_pdb_parse/ |
| B — RNA-FM embed (GPU) | $B_PID | $(kill -0 $B_PID 2>/dev/null && echo RUNNING || echo DONE) | .research/30_experiments/runs/B_rnafm_embeddings/ |
| C — MSA stats (CPU) | $C_PID | $(kill -0 $C_PID 2>/dev/null && echo RUNNING || echo DONE) | .research/30_experiments/runs/C_msa_stats/ |

## GPU snapshot

\`\`\`
$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.free --format=csv 2>/dev/null)
\`\`\`

## Last 5 log lines

\`\`\`
$(tail -5 $LOG)
\`\`\`
EOF
}

ITER=0
while kill -0 $A_PID 2>/dev/null || kill -0 $B_PID 2>/dev/null || kill -0 $C_PID 2>/dev/null; do
  ITER=$((ITER+1))
  write_status
  commit_push "chore(autonomous): incremental snapshot (iter $ITER)"
  sleep 90
done

# Final commit
echo "===== All jobs done $(date -u +%FT%TZ) =====" | tee -a $LOG
write_status
commit_push "feat(autonomous): final commit — all 3 parallel jobs complete"
