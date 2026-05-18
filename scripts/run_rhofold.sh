#!/usr/bin/env bash
# RhoFold+ second-predictor reproduction for the CPEB3 study.
#
# Install (from real repo github.com/ml4bio/RhoFold), download weights from
# HuggingFace, run single-sequence prediction on R1107 + R1108 + null_pos41,
# compare with DRfold2 results.

set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/rhofold_run.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
mkdir -p "$OUT"
exec > >(tee -a $LOG) 2>&1
echo "===== RhoFold+ pipeline started $(date -u +%FT%TZ) ====="

commit_push() {
  git add .research/30_experiments/runs/cpeb3_focused/SUMMARY.md \
          .research/30_experiments/runs/cpeb3_focused/*/rhofold/*.pdb \
          .research/40_findings/FINDINGS.md \
          2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -1
    git push origin main 2>&1 | tail -2
  fi
}

# ========== Step 1: install ==========
echo
echo "--- [1] install RhoFold+ ---"
cd /workspace/rna3d/models
if [ ! -d RhoFold ]; then
  git clone https://github.com/ml4bio/RhoFold.git 2>&1 | tail -3
fi
cd RhoFold

# Install via pip the deps from environment_linux.yaml (skip conda — our venv handles most)
# Inspect env file then pip-install best-effort
if [ -f envs/environment_linux.yaml ]; then
  echo "--- env file ---"
  head -40 envs/environment_linux.yaml
fi

# Install the package itself (it has setup.py)
pip install -q -e . 2>&1 | tail -5 || echo "(pip install -e . issue — will continue and import-test)"

# Common deps that may be missing
pip install -q einops dm-tree fair-esm omegaconf 2>&1 | tail -3

# RhoFold inference.py top-level imports simtk.openmm via relax module.
# We do NOT need relaxation (running with --single_seq_pred True for fast inference);
# patch the import to be optional so missing simtk does not crash.
INFER=/workspace/rna3d/models/RhoFold/inference.py
if grep -q "^from rhofold.relax.relax import AmberRelaxation$" $INFER; then
  echo "patching inference.py to make AmberRelaxation import optional"
  python3 -c "
p = '$INFER'
s = open(p).read()
s = s.replace(
    'from rhofold.relax.relax import AmberRelaxation',
    'try:\n    from rhofold.relax.relax import AmberRelaxation\nexcept Exception as _e:\n    AmberRelaxation = None\n    print(f\"[patch] AmberRelaxation unavailable: {_e}\")'
)
open(p, 'w').write(s)
print('patched')
"
fi
# Also patch any usage of AmberRelaxation(...) to be conditional
# (we add a guard in the main loop later if needed; for now just protect the import)

# ========== Step 2: download weights from HF ==========
echo
echo "--- [2] download RhoFold weights ---"
mkdir -p pretrained
if [ ! -f pretrained/RhoFold_pretrained.pt ]; then
  wget -q https://huggingface.co/cuhkaih/rhofold/resolve/main/rhofold_pretrained_params.pt -O pretrained/RhoFold_pretrained.pt
fi
ls -la pretrained/RhoFold_pretrained.pt 2>&1 | awk '{print $5, $9}'

# ========== Step 3: inference ==========
echo
echo "--- [3] inference on R1107 + R1108 + null_pos41 ---"
for tup in \
  "R1107_human:GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU" \
  "R1108_chimp:GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU" \
  "null_pos41:GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCCGCCAUUGCACUCCGGCUGCGAAUUCUGCU"
do
  tid="${tup%%:*}"; seq="${tup#*:}"
  td="$OUT/$tid/rhofold"
  mkdir -p "$td"
  if [ -f "$td/unrelaxed_model.pdb" ] || [ -f "$td/relaxed_1000_model.pdb" ]; then
    echo ">>> [skip] $tid already has rhofold output"
    continue
  fi
  fasta="$td/input.fasta"
  cat > "$fasta" <<EOF
>${tid}
${seq}
EOF
  echo ">>> running RhoFold on $tid <<<"
  cd /workspace/rna3d/models/RhoFold
  # --single_seq_pred True → no MSA, no relaxation
  python inference.py \
    --input_fas "$fasta" \
    --single_seq_pred True \
    --output_dir "$td" \
    --ckpt pretrained/RhoFold_pretrained.pt 2>&1 | tail -15
  ls "$td" 2>&1 | head -5
  # Verify we got a PDB before moving on
  if ! ls "$td"/*.pdb >/dev/null 2>&1; then
    echo "[WARN] no PDB produced for $tid — see error above"
  fi
done

# ========== Step 4: analyze and compare with DRfold2 ==========
echo
echo "--- [4] analysis ---"
python - <<'PYEND'
import numpy as np, datetime
from pathlib import Path
import biotite.structure.io as bsio

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
REF = Path("/workspace/rna3d/data/kaggle_raw/PDB_RNA/7qr3.cif")

def c1(p, chain=None):
    if not p.exists(): return None
    try:
        s = bsio.load_structure(str(p))
        if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
        if chain is not None and hasattr(s, "chain_id"):
            m = s.chain_id == chain
            if m.any(): s = s[m]
        m = s.atom_name == "C1'"
        return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None
    except Exception:
        return None

def kabsch(p, q):
    if p is None or q is None or p.shape != q.shape or p.shape[0] < 3: return float("nan")
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return float(np.sqrt(((al-q)**2).sum(-1).mean()))

def per_res(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al-q, axis=-1)

def rhofold_pdb(t):
    for name in ("relaxed_1000_model.pdb", "unrelaxed_model.pdb"):
        p = OUT/t/"rhofold"/name
        if p.exists(): return p
    return None

ref_c = c1(REF, chain="C")

# RhoFold predictions
rh = {t: c1(rhofold_pdb(t)) for t in ("R1107_human", "R1108_chimp", "null_pos41") if rhofold_pdb(t) is not None}
print(f"RhoFold loaded: {list(rh.keys())}")
for t, x in rh.items():
    print(f"  {t}: shape {None if x is None else x.shape}")

# Append to SUMMARY.md
with open(OUT/"SUMMARY.md", "a") as f:
    f.write("\n\n---\n\n# RhoFold+ second-predictor results (added " + datetime.datetime.now(datetime.UTC).isoformat() + ")\n\n")
    f.write("## RhoFold+ RMSD vs 7QR3 chain C\n\n| target | RMSD (A) |\n|---|---|\n")
    for t, x in rh.items():
        if x is not None and ref_c is not None and x.shape == ref_c.shape:
            f.write(f"| {t} | {kabsch(x, ref_c):.3f} |\n")
    if "R1107_human" in rh and "R1108_chimp" in rh:
        real_delta = per_res(rh["R1107_human"], rh["R1108_chimp"])
        top5_r = np.argsort(real_delta)[-5:][::-1].tolist()
        f.write(f"\n## RhoFold+ REAL case (R1107 vs R1108, pos 30 A->G)\n\n")
        f.write(f"- Kabsch RMSD: {kabsch(rh['R1107_human'], rh['R1108_chimp']):.3f} A\n")
        f.write(f"- Position 30 delta: {real_delta[29]:.3f} A\n")
        f.write(f"- mean/max delta: {real_delta.mean():.2f} / {real_delta.max():.2f} A\n")
        f.write(f"- Top-5 divergent residues: " + ", ".join(f"pos {i+1} ({real_delta[i]:.2f}A)" for i in top5_r) + "\n")
    if "null_pos41" in rh and "R1108_chimp" in rh:
        null_delta = per_res(rh["null_pos41"], rh["R1108_chimp"])
        top5_n = np.argsort(null_delta)[-5:][::-1].tolist()
        f.write(f"\n## RhoFold+ NULL case (pos 41 A->C)\n\n")
        f.write(f"- Kabsch RMSD: {kabsch(rh['null_pos41'], rh['R1108_chimp']):.3f} A\n")
        f.write(f"- Position 41 delta: {null_delta[40]:.3f} A\n")
        f.write(f"- Top-5 divergent residues: " + ", ".join(f"pos {i+1} ({null_delta[i]:.2f}A)" for i in top5_n) + "\n")
        if "R1107_human" in rh:
            top5_real_set = set(int(i)+1 for i in np.argsort(real_delta)[-5:])
            top5_null_set = set(int(i)+1 for i in np.argsort(null_delta)[-5:])
            overlap = sorted(top5_real_set & top5_null_set)
            f.write(f"\n## RhoFold+ cross-validation verdict\n\n- Top-5 overlap (REAL & NULL): **{overlap}** ({len(overlap)} positions)\n")
            f.write(f"- DRfold2 had overlap=1. If RhoFold+ overlap is also small (<=2) AND its REAL top-5 includes residues near 9 or 60 → **F-001 is cross-model reproduced**.\n")
PYEND
commit_push "feat(cpeb3-rhofold): RhoFold+ predictions + cross-model comparison with DRfold2"

echo
echo "===== RhoFold+ pipeline finished $(date -u +%FT%TZ) ====="
