#!/usr/bin/env bash
# CPEB3 sanity-check pipeline — second pod run.
# Runs:
#   (a) Null control: re-predict R1108 with ONE arbitrary mutation at a residue
#       OUTSIDE the P1/P1.1 region (we use position 41, in J3-P4 — peripheral).
#       If DRfold2's top-5-divergent-residues stays at 9, 60, 51, 22, 24 →
#       F-001 is a model prior (BAD). If top-5 shifts to surround pos 41 →
#       F-001 is real (GOOD).
#   (b) Inspect the previous run's .ret files for any per-residue confidence
#       score we can extract.
#   (c) Try to install + run RhoFold+ as a second predictor (best-effort).
#
# Auto-commits + auto-pushes incrementally. Watchdog stops the pod at end.
set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/cpeb3_sanity.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
mkdir -p "$OUT/null_pos41"
exec > >(tee -a $LOG) 2>&1
echo "===== CPEB3 sanity pipeline started $(date -u +%FT%TZ) ====="

commit_push() {
  git add .research/ 2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -2
    git push origin main 2>&1 | tail -2
  fi
}

# ============== (b) Confidence channel — inspect .ret files first (fast, CPU) ==============
echo
echo "--- (b) inspect .ret files for per-residue confidence ---"
python - <<'PYEND'
from pathlib import Path
import numpy as np
OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
for target in ("R1107_human", "R1108_chimp"):
    d = OUT / target / "drfold2" / "rets_dir"
    if not d.exists(): continue
    rets = sorted(d.glob("*.ret"))[:3]
    print(f"\n{target}: {len(list(d.glob('*.ret')))} .ret files; inspecting {len(rets)}")
    for r in rets:
        size = r.stat().st_size
        print(f"  {r.name} ({size} bytes)")
        try:
            # Try common formats
            arr = np.load(str(r), allow_pickle=True)
            print(f"    np.load OK: shape={getattr(arr, 'shape', '?')} dtype={getattr(arr, 'dtype', '?')}")
            if isinstance(arr, np.ndarray) and arr.dtype == object:
                print(f"    object array; first item type: {type(arr.item()).__name__ if arr.size == 1 else type(arr.flat[0]).__name__}")
        except Exception:
            try:
                import pickle
                with open(r, "rb") as f:
                    obj = pickle.load(f)
                print(f"    pickle OK: type={type(obj).__name__}")
                if hasattr(obj, 'keys'):
                    print(f"    keys: {list(obj.keys())[:10]}")
                elif isinstance(obj, (list, tuple)):
                    print(f"    len={len(obj)}, first type={type(obj[0]).__name__}")
            except Exception as e:
                # Last resort: peek as text
                with open(r, "rb") as f:
                    head = f.read(80)
                print(f"    not np/pickle; head bytes: {head[:40]}")
PYEND
commit_push "feat(cpeb3-sanity): inspected .ret files for confidence-channel format"

# ============== (a) Null control: mutate R1108 at position 41 ==============
echo
echo "--- (a) NULL CONTROL: mutate R1108 at position 41 (G->C, peripheral residue) ---"
# R1108: GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU
# position 41 (1-indexed) → 0-indexed 40 → currently 'A' (the 'A' in UCAGCC**A**UUG)
# Let's mutate to C as a control single-nt change in a non-P1 region.
NULL_SEQ="GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCCGCCAUUGCACUCCGGCUGCGAAUUCUGCU"
# Sanity: confirm the change is at pos 41 only
python -c "
ref = 'GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU'
new = '$NULL_SEQ'
assert len(ref) == len(new) == 69
diffs = [i for i in range(len(ref)) if ref[i] != new[i]]
print(f'null mutation diff at 1-indexed: {[i+1 for i in diffs]} ({[ref[i] for i in diffs]} -> {[new[i] for i in diffs]})')
assert diffs == [40], f'expected only position 41 to change, got {[i+1 for i in diffs]}'
"
cat > "$OUT/null_pos41/input.fasta" <<EOF
>R1108_chimp_NULL_pos41_A2C
$NULL_SEQ
EOF
mkdir -p "$OUT/null_pos41/drfold2"
cd /workspace/rna3d/models/DRfold2
python DRfold_infer.py "$OUT/null_pos41/input.fasta" "$OUT/null_pos41/drfold2" 2>&1 | tail -10
ls -la "$OUT/null_pos41/drfold2/folds/" 2>&1 | head
commit_push "feat(cpeb3-sanity): DRfold2 null-control mutation @ pos 41 (R1108 with A2C, peripheral)"

# ============== (a-cont'd) Compute null vs R1108 deltas + compare to A2G case ==============
echo
echo "--- (a-cont'd) compute null delta + compare to position-30 case ---"
python - <<'PYEND'
import numpy as np
from pathlib import Path
import biotite.structure.io as bsio

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
def c1(p):
    s = bsio.load_structure(str(p))
    if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
    m = s.atom_name == "C1'"
    return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None
def kabsch_per_res(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al - q, axis=-1)

def best(t):
    g = list((OUT / t / "drfold2" / "folds").glob("opt_0_*.pdb"))
    return c1(g[0]) if g else None

chimp = best("R1108_chimp")          # original R1108 (chimp, has G at pos 30)
human = best("R1107_human")          # original R1107 (human, has A at pos 30 — the real biological variant)
null_p41 = best("null_pos41")        # R1108 with A->C at pos 41 (null peripheral mutation)
assert all(x is not None for x in (chimp, human, null_p41)), "missing prediction(s)"

# real A2G case (pos 30): human vs chimp
real = kabsch_per_res(human, chimp)
# null A2C case (pos 41): null vs chimp
null = kabsch_per_res(null_p41, chimp)

print("=== POSITION-WISE DELTA COMPARISON ===")
print(f"REAL (R1107 vs R1108, single-nt diff at pos 30): mean={real.mean():.2f}, max={real.max():.2f}")
print(f"NULL (R1108_pos41_A2C vs R1108, single-nt diff at pos 41): mean={null.mean():.2f}, max={null.max():.2f}")
print()
print("Top-5 divergent residues, REAL case (pos 30 A->G):")
top5_real = np.argsort(real)[-5:][::-1]
for i in top5_real:
    print(f"  pos {i+1}: {real[i]:.2f} A")
print()
print("Top-5 divergent residues, NULL case (pos 41 A->C):")
top5_null = np.argsort(null)[-5:][::-1]
for i in top5_null:
    print(f"  pos {i+1}: {null[i]:.2f} A")
print()
overlap = set([int(i)+1 for i in top5_real]) & set([int(i)+1 for i in top5_null])
print(f"Overlap of top-5 between REAL and NULL: {sorted(overlap)}")
print()
print("VERDICT:")
print("  if overlap is large (e.g. 4-5 positions identical) → top-5 is a MODEL PRIOR; F-001 is artifact.")
print("  if overlap is small (e.g. 0-2 positions identical) AND the NULL top-5 surrounds pos 41 → F-001 is REAL.")

# Append result to SUMMARY.md
with open(OUT/"SUMMARY.md", "a") as f:
    f.write("\n\n---\n\n## Sanity check (a): NULL CONTROL — mutation at peripheral position 41 (A->C)\n\n")
    f.write(f"REAL case (R1107 vs R1108 — pos 30 A->G, biologically functional):\n")
    f.write(f"- mean delta = {real.mean():.2f} A, max = {real.max():.2f} A\n")
    f.write(f"- top-5 residues: {[int(i)+1 for i in top5_real]}\n\n")
    f.write(f"NULL case (R1108_p41_A2C vs R1108 — pos 41 A->C, peripheral):\n")
    f.write(f"- mean delta = {null.mean():.2f} A, max = {null.max():.2f} A\n")
    f.write(f"- top-5 residues: {[int(i)+1 for i in top5_null]}\n\n")
    f.write(f"Overlap of top-5 (REAL ∩ NULL): {sorted(overlap)}\n\n")
    if len(overlap) >= 3:
        f.write("→ **F-001 IS LIKELY AN ARTIFACT**: top-5 divergent residues are largely the same regardless of where the mutation actually sits. The 'cascade to P1' is a DRfold2 prior, not a response to the position-30 mutation specifically.\n")
    else:
        f.write("→ **F-001 SURVIVES THE NULL CONTROL**: top-5 changes meaningfully between mutations. The position-30 → P1 cascade appears to be specific to the position-30 mutation, supporting a real biology↔model alignment.\n")
PYEND
commit_push "feat(cpeb3-sanity): null-control analysis — top-5 overlap between REAL (pos30) and NULL (pos41) cases"

# ============== (c) Best-effort RhoFold+ install ==============
echo
echo "--- (c) try RhoFold+ install via SSH URL ---"
cd /workspace/rna3d/models
for url in \
  "git@github.com:RFOLD/RhoFold.git" \
  "git@github.com:biomap-research/rhofold.git" \
  "https://github.com/sokrypton/rhofold.git"
do
  echo "trying $url"
  git clone "$url" rhofold_test 2>&1 | tail -3
  if [ -d "rhofold_test" ] && [ -n "$(ls rhofold_test 2>/dev/null)" ]; then
    mv rhofold_test rhofold
    echo "→ cloned from $url"
    break
  fi
  rm -rf rhofold_test
done
if [ -d "rhofold" ]; then
  cd rhofold && ls && pip install -q -r requirements.txt 2>&1 | tail -3
  cd /workspace/rna3d
  echo "RhoFold present; wiring inference is deferred to next run"
else
  echo "RhoFold not obtainable via tried URLs; defer to AF3 server (manual)"
fi
commit_push "chore(cpeb3-sanity): rhofold install attempt + final"

echo
echo "===== CPEB3 sanity pipeline finished $(date -u +%FT%TZ) ====="
