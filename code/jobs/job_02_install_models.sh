#!/usr/bin/env bash
set +e
source /workspace/rna3d/venv/bin/activate
cd /workspace/rna3d
OUT=/workspace/rna3d/.research/30_experiments/runs/02_model_installation
mkdir -p "$OUT/results"
exec > >(tee "$OUT/INSTALL.md") 2>&1

echo "# Model installation report — $(date -u +%FT%TZ)"
echo ""
echo "## RibonanzaNet"
pip install -q huggingface_hub einops 2>&1 | tail -3
pip install -q git+https://github.com/Shujun-He/RibonanzaNet.git 2>&1 | tail -5 || echo "git+pip failed"
python -c "from ribonanzanet import RibonanzaNet; print('RibonanzaNet import OK')" 2>&1 | tail -2
if ! python -c "import ribonanzanet" 2>/dev/null; then
  cd /workspace/rna3d/models
  git clone https://github.com/Shujun-He/RibonanzaNet.git ribonanzanet_repo 2>&1 | tail -2
  cd ribonanzanet_repo 2>/dev/null && pip install -q -e . 2>&1 | tail -3
  cd /workspace/rna3d
fi
echo ""
echo "## RhoFold+"
cd /workspace/rna3d/models
git clone https://github.com/RFOLD/RhoFold.git rhofold_repo 2>&1 | tail -3
echo ""
echo "## GPU memory"
nvidia-smi --query-gpu=memory.used,memory.free --format=csv 2>&1
echo "DONE_02 $(date -u +%FT%TZ)"
