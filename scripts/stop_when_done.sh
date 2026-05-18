#!/usr/bin/env bash
# Watchdog: waits for pod_launch.sh to finish, then stops the pod via runpodctl.
# Idempotent — safe to launch multiple times.
POD_ID="6sp7kzqjg5ypb3"
LOG=/workspace/rna3d/stop_watcher.log

echo "[$(date -u +%FT%TZ)] watchdog starting, will stop pod $POD_ID when training finishes" >> $LOG

# Wait while training is still running
while pgrep -af "pod_launch.sh|train.py" | grep -v grep > /dev/null 2>&1; do
  sleep 60
done

# Grace period for final commits/push
echo "[$(date -u +%FT%TZ)] training done, sleeping 120s grace then stopping pod" >> $LOG
sleep 120

echo "[$(date -u +%FT%TZ)] stopping pod $POD_ID" >> $LOG
runpodctl stop pod "$POD_ID" 2>&1 >> $LOG

# Belt-and-suspenders: also try poweroff
sleep 30
echo "[$(date -u +%FT%TZ)] fallback poweroff" >> $LOG
poweroff 2>&1 >> $LOG
