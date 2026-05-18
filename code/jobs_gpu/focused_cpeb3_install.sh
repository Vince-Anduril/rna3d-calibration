#!/usr/bin/env bash
# Install candidate 3D RNA predictors on the pod for the focused CPEB3 study.
# Best-effort: each model is independent, failure of one does not block the others.
# Plan: get at least TWO working out of {DRfold2, RhoFold+, trRosettaRNA}.
#
# Run from /workspace/rna3d/ on the pod. Activates the venv first.

set +e  # do not abort the whole script on individual install failures
cd /workspace/rna3d
source venv/bin/activate

INSTALL_LOG=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/INSTALL.md
mkdir -p "$(dirname "$INSTALL_LOG")"
exec > >(tee "$INSTALL_LOG") 2>&1

echo "# CPEB3 focused study — install report"
echo
echo "Started: $(date -u +%FT%TZ)"
echo

# ---------------------------------------------------------------------------
echo "## 1. DRfold2"
# Published PLOS Biology 2025. Best reported perf on CPEB3 chimp (RMSD 2.72 Å).
# Repo: https://github.com/leeyang/DRfold2 (verify URL at install time).
echo
DRFOLD2_REPO_URL="https://github.com/leeyang/DRfold2.git"  # PROBABLE — verify before relying on this
mkdir -p /workspace/rna3d/models
cd /workspace/rna3d/models
if [ ! -d "DRfold2" ]; then
  git clone "$DRFOLD2_REPO_URL" DRfold2 2>&1 | tail -5
fi
if [ -d "DRfold2" ]; then
  cd DRfold2
  echo "DRfold2 repo present. Files:"
  ls -la | head -20
  # Try standard install steps; the README in the repo is authoritative.
  pip install -q -e . 2>&1 | tail -3 || echo "DRfold2 editable install failed — read README"
  python -c "import drfold2; print('DRfold2 importable')" 2>&1 | head -2 || \
    echo "TODO: import path may be different — check repo's example/inference.py"
  cd /workspace/rna3d
else
  echo "DRfold2 clone failed — try alternative URL"
fi
echo

# ---------------------------------------------------------------------------
echo "## 2. RhoFold+"
# Open-source, Nature Methods 2024.
# Real repo URL is uncertain; try a few:
echo
cd /workspace/rna3d/models
for url in \
  "https://github.com/biomap-research/rhofold.git" \
  "https://github.com/RFOLD/RhoFold-PyTorch.git" \
  "https://github.com/RFOLD/RhoFold.git"
do
  echo "Trying: $url"
  git clone "$url" rhofold_test 2>&1 | tail -2
  if [ -d "rhofold_test" ] && [ "$(ls rhofold_test)" ]; then
    mv rhofold_test rhofold
    echo "  → cloned from $url"
    break
  fi
  rm -rf rhofold_test
done
if [ -d "rhofold" ]; then
  cd rhofold
  echo "RhoFold contents:"
  ls -la | head -20
  pip install -q -r requirements.txt 2>&1 | tail -3 || echo "(no requirements.txt or pip failed)"
  pip install -q -e . 2>&1 | tail -3 || echo "(no setup.py or pip failed)"
  cd /workspace/rna3d
else
  echo "RhoFold not found via any URL — manual install needed (search github)"
fi
echo

# ---------------------------------------------------------------------------
echo "## 3. trRosettaRNA (fallback)"
echo
cd /workspace/rna3d/models
for url in \
  "https://github.com/yangzhanglab/trRosettaRNA.git" \
  "https://github.com/zhng-lab/trRosettaRNA.git"
do
  echo "Trying: $url"
  git clone "$url" trrosetta_test 2>&1 | tail -2
  if [ -d "trrosetta_test" ] && [ "$(ls trrosetta_test)" ]; then
    mv trrosetta_test trrosettaRNA
    break
  fi
  rm -rf trrosetta_test
done
[ -d "trrosettaRNA" ] && echo "trRosettaRNA cloned" || echo "trRosettaRNA not found"
echo

# ---------------------------------------------------------------------------
echo "## 4. AF3 server (no install; manual submission)"
echo
cat <<'EOF'
For AlphaFold Server (AF3), no install. Process:
  1. Go to https://alphafoldserver.com/ (free Google account, ~20 jobs/day).
  2. Submit each sequence (R1107, R1108) as a single-RNA job.
  3. Download the JSON result.
  4. Place files at:
     /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/R1107_human/af3_server.json
     /workspace/rna3d/.research/30_experiments/runs/cpeb3_focused/R1108_chimp/af3_server.json
EOF
echo

# ---------------------------------------------------------------------------
echo "## 5. GPU sanity"
nvidia-smi --query-gpu=name,memory.used,memory.free --format=csv,noheader
echo
echo "Finished: $(date -u +%FT%TZ)"
