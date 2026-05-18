#!/usr/bin/env bash
# Full CPEB3 focused pipeline — single script combining inference + sanity checks.
# Idempotent: if folds/opt_0_*.pdb already exist for a target, skips re-inference.
# Auto-commits + auto-pushes after each phase. Watchdog should be launched separately.

set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/cpeb3_full.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
mkdir -p "$OUT"
exec > >(tee -a $LOG) 2>&1
echo "===== CPEB3 FULL pipeline started $(date -u +%FT%TZ) ====="

commit_push() {
  git add .gitignore .research/30_experiments/runs/cpeb3_focused/SUMMARY.md \
          .research/30_experiments/runs/cpeb3_focused/*/input.fasta \
          .research/30_experiments/runs/cpeb3_focused/*/drfold2/folds/opt_0_*.pdb \
          .research/40_findings/FINDINGS.md \
          2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -1
    git push origin main 2>&1 | tail -2
  fi
}

# Inference helper: only runs DRfold2 if folds/opt_0_*.pdb is missing
infer_if_needed() {
  local target_id="$1"
  local sequence="$2"
  local td="$OUT/$target_id"
  mkdir -p "$td"
  if ls "$td/drfold2/folds/opt_0_"*.pdb >/dev/null 2>&1; then
    echo "  [skip] $target_id already has folds/opt_0_*.pdb"
    return
  fi
  echo "  [run] $target_id"
  cat > "$td/input.fasta" <<EOF
>$target_id
$sequence
EOF
  mkdir -p "$td/drfold2"
  (cd /workspace/rna3d/models/DRfold2 && python DRfold_infer.py "$td/input.fasta" "$td/drfold2" 2>&1 | tail -5)
  ls "$td/drfold2/folds/" 2>&1 | head
}

# ========== Phase 1: ensure originals exist ==========
echo
echo "--- [1] R1107 + R1108 inference (skip if folds exist) ---"
infer_if_needed "R1107_human" "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"
infer_if_needed "R1108_chimp" "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"

# ========== Phase 2: NULL CONTROL — mutation at pos 41 (peripheral) ==========
echo
echo "--- [2] NULL CONTROL: R1108 with A->C at position 41 ---"
# pos 41 0-indexed 40: original is 'A' (in UCAGCC**A**UUG); mutate to C
infer_if_needed "null_pos41" "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCCGCCAUUGCACUCCGGCUGCGAAUUCUGCU"

# ========== Phase 3: full analysis ==========
echo
echo "--- [3] analysis: RMSD vs 7QR3 + REAL vs NULL discrimination ---"
python - <<'PYEND'
import sys, datetime
import numpy as np
from pathlib import Path
import biotite.structure.io as bsio

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
REF = Path("/workspace/rna3d/data/kaggle_raw/PDB_RNA/7qr3.cif")

def c1(p, chain=None):
    try:
        s = bsio.load_structure(str(p))
        if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
        if chain is not None and hasattr(s, "chain_id"):
            m = s.chain_id == chain
            if m.any(): s = s[m]
        m = s.atom_name == "C1'"
        return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None
    except Exception as e:
        print(f"  c1 failed on {p.name}: {e}", file=sys.stderr); return None

def kabsch_rmsd(p, q):
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

def best_pdb(t):
    g = sorted((OUT/t/"drfold2/folds").glob("opt_0_*.pdb"))
    if not g: return None, None
    return c1(g[0]), g[0].name

ref_c = c1(REF, chain="C")
ref_d = c1(REF, chain="D")
print(f"REF 7QR3 chain C: {None if ref_c is None else ref_c.shape}")
print(f"REF 7QR3 chain D: {None if ref_d is None else ref_d.shape}")

h, h_name = best_pdb("R1107_human")
c, c_name = best_pdb("R1108_chimp")
n, n_name = best_pdb("null_pos41")

print(f"R1107_human best: {h_name}, shape {None if h is None else h.shape}")
print(f"R1108_chimp best: {c_name}, shape {None if c is None else c.shape}")
print(f"null_pos41 best:  {n_name}, shape {None if n is None else n.shape}")

# Metrics
real_delta = per_res(h, c) if (h is not None and c is not None) else None
null_delta = per_res(n, c) if (n is not None and c is not None) else None

with open(OUT/"SUMMARY.md", "w") as f:
    f.write("# CPEB3 focused study — DRfold2 results + NULL CONTROL\n\n")
    f.write(f"Run completed: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
    f.write("## Setup\n\n")
    f.write("- Predictor: DRfold2 (PLOS Biology 2025), 4-config ensemble cfg_95/96/97/99\n")
    f.write("- Targets:\n")
    f.write("  - R1107 human CPEB3 (pos 30 = A)\n")
    f.write("  - R1108 chimpanzee CPEB3 (pos 30 = G)\n")
    f.write("  - null_pos41: R1108 with A->C at position 41 (peripheral, J3-P4 region)\n")
    f.write("- Reference: PDB 7QR3 chain C\n\n")
    f.write("## RMSD vs 7QR3\n\n")
    f.write("| target | RMSD vs chain C (A) | RMSD vs chain D (A) |\n|---|---|---|\n")
    if h is not None:
        f.write(f"| R1107 human | {kabsch_rmsd(h, ref_c):.3f} | {kabsch_rmsd(h, ref_d):.3f} |\n")
    if c is not None:
        f.write(f"| R1108 chimp | {kabsch_rmsd(c, ref_c):.3f} | {kabsch_rmsd(c, ref_d):.3f} |\n")
    if n is not None:
        f.write(f"| null_pos41 | {kabsch_rmsd(n, ref_c):.3f} | {kabsch_rmsd(n, ref_d):.3f} |\n")
    f.write("\n## Discrimination & NULL CONTROL\n\n")
    if real_delta is not None:
        f.write(f"### REAL case (R1107 vs R1108, pos 30 A->G, biologically functional)\n\n")
        f.write(f"- Kabsch RMSD between predictions: **{kabsch_rmsd(h, c):.3f} A**\n")
        f.write(f"- Position 30 (mutation site) per-residue delta: **{real_delta[29]:.3f} A**\n")
        f.write(f"- mean/max delta: {real_delta.mean():.2f} / {real_delta.max():.2f} A\n")
        top5_real = np.argsort(real_delta)[-5:][::-1].tolist()
        f.write(f"- Top-5 divergent residues (1-indexed): " +
                ", ".join(f"pos {i+1} ({real_delta[i]:.2f}A)" for i in top5_real) + "\n\n")
    if null_delta is not None:
        f.write(f"### NULL case (null_pos41 vs R1108, pos 41 A->C, peripheral)\n\n")
        f.write(f"- Kabsch RMSD between predictions: **{kabsch_rmsd(n, c):.3f} A**\n")
        f.write(f"- Position 41 (null mutation site) per-residue delta: **{null_delta[40]:.3f} A**\n")
        f.write(f"- mean/max delta: {null_delta.mean():.2f} / {null_delta.max():.2f} A\n")
        top5_null = np.argsort(null_delta)[-5:][::-1].tolist()
        f.write(f"- Top-5 divergent residues (1-indexed): " +
                ", ".join(f"pos {i+1} ({null_delta[i]:.2f}A)" for i in top5_null) + "\n\n")
        if real_delta is not None:
            top5_real_set = set(int(i)+1 for i in np.argsort(real_delta)[-5:])
            top5_null_set = set(int(i)+1 for i in np.argsort(null_delta)[-5:])
            overlap = sorted(top5_real_set & top5_null_set)
            f.write(f"### VERDICT\n\n")
            f.write(f"- Overlap of top-5 (REAL & NULL): **{overlap}** ({len(overlap)} positions)\n\n")
            if len(overlap) >= 4:
                f.write("**F-001 IS LIKELY AN ARTIFACT.** The top-5 divergent residues are essentially the same regardless of where the mutation actually sits. The 'pos 30 -> P1 cascade' is a DRfold2 structural prior, not a response to the position-30 mutation specifically.\n")
            elif len(overlap) >= 3:
                f.write("**F-001 IS PARTIALLY ARTIFACTUAL.** Strong overlap suggests DRfold2 has a structural prior that biases divergence to certain residues. F-001 should be re-framed as 'DRfold2 reliably predicts deltas in [residues 9, 60, 51, 22, 24] regardless of where small mutations occur'.\n")
            else:
                f.write("**F-001 SURVIVES THE NULL CONTROL.** Top-5 changes meaningfully between mutations. The position-30 -> P1 cascade appears specific to the position-30 mutation, supporting a real biology/model alignment.\n")
print("SUMMARY.md written")
PYEND
commit_push "feat(cpeb3): full results (R1107+R1108) + NULL CONTROL at pos 41 — verdict on F-001"

echo
echo "===== CPEB3 FULL pipeline finished $(date -u +%FT%TZ) ====="
