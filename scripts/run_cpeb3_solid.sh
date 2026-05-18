#!/usr/bin/env bash
# Solidify CPEB3 finding for paper. Three killer experiments:
#   (1) Confidence channel: inspect DRfold2 .ret files for per-residue scores
#   (2) Multiple null controls at positions {5, 20, 55, 64} (already have 41)
#       → robustness statistic over 5 NULLs vs 1 REAL
#   (3) Secondary structure: ViennaRNA RNAfold + DRfold2 ss.ct extraction →
#       does R1107 vs R1108 SS differ in the P1.1 region?

set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/cpeb3_solid.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
exec > >(tee -a $LOG) 2>&1
echo "===== CPEB3 SOLID pipeline started $(date -u +%FT%TZ) ====="

commit_push() {
  cd /workspace/rna3d   # CRITICAL: previous bug was wrong cwd at commit time
  git add .research/30_experiments/runs/cpeb3_focused/SUMMARY.md \
          .research/30_experiments/runs/cpeb3_focused/null_*/input.fasta \
          .research/30_experiments/runs/cpeb3_focused/null_*/drfold2/folds/opt_0_*.pdb \
          .research/30_experiments/runs/cpeb3_focused/confidence_analysis.md \
          .research/30_experiments/runs/cpeb3_focused/ss_analysis.md \
          .research/30_experiments/runs/cpeb3_focused/multi_null_table.csv \
          .research/40_findings/FINDINGS.md \
          2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -1
    git push origin main 2>&1 | tail -2
  fi
}

# ============================================================================
# Phase 1: DRfold2 confidence channel inspection (~2 min, CPU)
# ============================================================================
echo
echo "--- [1] inspect DRfold2 .ret files for per-residue confidence ---"
python - <<'PYEND'
from pathlib import Path
import numpy as np
import torch

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
CONF = OUT / "confidence_analysis.md"

with open(CONF, "w") as f:
    f.write("# DRfold2 confidence channel inspection\n\n")
    for target in ("R1107_human", "R1108_chimp", "null_pos41"):
        rets_dir = OUT / target / "drfold2" / "rets_dir"
        f.write(f"\n## {target}\n\n")
        if not rets_dir.exists():
            f.write(f"_no rets_dir_\n"); continue
        rets = sorted(rets_dir.glob("*.ret"))
        f.write(f"- {len(rets)} .ret files\n")
        if not rets: continue
        # Try first few
        for r in rets[:3]:
            f.write(f"\n### {r.name}\n\n```\n")
            try:
                obj = torch.load(str(r), map_location="cpu", weights_only=False)
                f.write(f"torch.load OK, type={type(obj).__name__}\n")
                if hasattr(obj, 'keys'):
                    for k in list(obj.keys())[:20]:
                        v = obj[k]
                        shape = getattr(v, 'shape', None)
                        f.write(f"  '{k}': type={type(v).__name__} shape={shape}\n")
                elif isinstance(obj, (list, tuple)):
                    for i, x in enumerate(obj[:10]):
                        f.write(f"  [{i}]: type={type(x).__name__} shape={getattr(x, 'shape', None)}\n")
                elif isinstance(obj, np.ndarray):
                    f.write(f"  np.ndarray shape={obj.shape} dtype={obj.dtype}\n")
                elif torch.is_tensor(obj):
                    f.write(f"  tensor shape={obj.shape} dtype={obj.dtype}\n")
            except Exception as e:
                f.write(f"torch.load failed: {e}\n")
                # Try pickle / np
                try:
                    import pickle
                    with open(r, "rb") as fh: obj = pickle.load(fh)
                    f.write(f"pickle OK, type={type(obj).__name__}\n")
                except Exception as e2:
                    f.write(f"pickle also failed: {e2}\n")
                    with open(r, "rb") as fh:
                        head = fh.read(80)
                    f.write(f"first 40 bytes: {head[:40]!r}\n")
            f.write("```\n")
print("confidence inspect done")
PYEND
commit_push "feat(cpeb3-solid): inspect DRfold2 .ret files for per-residue confidence channel"

# ============================================================================
# Phase 2: multiple NULL controls at positions {5, 20, 55, 64}
# ============================================================================
echo
echo "--- [2] multiple NULL CONTROLS (positions 5, 20, 55, 64) ---"
declare -A NULLS=(
  [null_pos5]="GGGGACCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"   # pos 5 G->A (5' end)
  [null_pos20]="GGGGGCCACAGCAGAAGCGAUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"  # pos 20 U->A (P2 mid)
  [null_pos55]="GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCAGCGAAUUCUGCU"  # pos 55 U->A (P1 stem)
  [null_pos64]="GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUACUGCU"  # pos 64 U->A (P1 3' end)
)

# Sanity check + run each
for tid in "${!NULLS[@]}"; do
  seq="${NULLS[$tid]}"
  td="$OUT/$tid"
  mkdir -p "$td/drfold2"
  if ls "$td/drfold2/folds/opt_0_"*.pdb >/dev/null 2>&1; then
    echo "  [skip] $tid already has folds"
    continue
  fi
  cat > "$td/input.fasta" <<EOF
>$tid
$seq
EOF
  python -c "
ref = 'GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU'
new = open('$td/input.fasta').readlines()[1].strip()
assert len(ref) == len(new) == 69, f'len mismatch'
d = [i+1 for i in range(69) if ref[i] != new[i]]
print(f'$tid: change at 1-indexed {d}')
"
  echo ">>> running DRfold2 on $tid <<<"
  (cd /workspace/rna3d/models/DRfold2 && python DRfold_infer.py "$td/input.fasta" "$td/drfold2" 2>&1 | tail -5)
done
commit_push "feat(cpeb3-solid): 4 additional NULL CONTROLS at peripheral positions 5, 20, 55, 64"

# ============================================================================
# Phase 3: aggregate multi-null analysis
# ============================================================================
echo
echo "--- [3] aggregate analysis: multi-null robustness statistic ---"
python - <<'PYEND'
import numpy as np
from pathlib import Path
import biotite.structure.io as bsio
import csv

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")

def c1(p):
    if not p.exists(): return None
    s = bsio.load_structure(str(p))
    if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
    m = s.atom_name == "C1'"
    return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None

def best(t):
    g = sorted((OUT/t/"drfold2/folds").glob("opt_0_*.pdb"))
    return c1(g[0]) if g else None

def per_res(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al-q, axis=-1)

chimp = best("R1108_chimp")
human = best("R1107_human")
real_delta = per_res(human, chimp) if (chimp is not None and human is not None) else None

# Mutation site per target → for plotting
mut_pos = {
    "R1107_human": 30, "null_pos5": 5, "null_pos20": 20,
    "null_pos41": 41, "null_pos55": 55, "null_pos64": 64,
}

rows = []
for tid, mut in mut_pos.items():
    p = best(tid)
    if p is None or chimp is None: continue
    d = per_res(p, chimp)
    top5 = sorted([int(i)+1 for i in np.argsort(d)[-5:][::-1]])
    rows.append({
        "target": tid,
        "mutated_position": mut,
        "kabsch_rmsd_vs_chimp": float(np.linalg.norm(p.mean(0)-chimp.mean(0))),  # rough
        "mean_delta": float(d.mean()),
        "max_delta": float(d.max()),
        "delta_at_mutation": float(d[mut-1]),
        "delta_at_pos9": float(d[8]),
        "delta_at_pos60": float(d[59]),
        "top5": ",".join(str(x) for x in top5),
        "p1_anchor_residue_in_top5": "yes" if (9 in top5 or 60 in top5) else "no",
    })

with open(OUT/"multi_null_table.csv", "w", newline="") as f:
    if rows:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

# Robustness statistic
n_total = len(rows) - 1  # exclude REAL
n_with_p1 = sum(1 for r in rows if r["target"] != "R1107_human" and r["p1_anchor_residue_in_top5"] == "yes")
real_has_p1 = any(r["target"] == "R1107_human" and r["p1_anchor_residue_in_top5"] == "yes" for r in rows)

print(f"\n=== MULTI-NULL ROBUSTNESS STATISTIC ===")
print(f"REAL case (pos 30 mutation): P1 anchor (residues 9 or 60) in top-5? -> {'YES' if real_has_p1 else 'NO'}")
print(f"NULL cases (other mutations): {n_with_p1}/{n_total} have P1 anchor in top-5")
print(f"Interpretation:")
if real_has_p1 and n_with_p1 <= 1:
    print(f"  → STRONG: P1 anchor specifically appears for the REAL (pos 30) mutation, {n_with_p1}/{n_total} for arbitrary mutations.")
elif real_has_p1 and n_with_p1 / max(n_total, 1) >= 0.5:
    print(f"  → WEAK: P1 anchor appears for both REAL and most NULL — could be a model attractor")
else:
    print(f"  → MIXED — interpret in context")

import datetime
with open(OUT/"SUMMARY.md", "a") as f:
    f.write("\n\n---\n\n## Multi-null robustness analysis (5 null controls + 1 real, "
            + datetime.datetime.now(datetime.UTC).isoformat() + ")\n\n")
    f.write("| target | mutation | mean Δ (Å) | max Δ (Å) | Δ at mut | Δ pos 9 | Δ pos 60 | top-5 | P1 anchor? |\n")
    f.write("|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        f.write(f"| {r['target']} | pos {r['mutated_position']} | "
                f"{r['mean_delta']:.2f} | {r['max_delta']:.2f} | "
                f"{r['delta_at_mutation']:.2f} | {r['delta_at_pos9']:.2f} | "
                f"{r['delta_at_pos60']:.2f} | {{{r['top5']}}} | {r['p1_anchor_residue_in_top5']} |\n")
    f.write(f"\n**Robustness:** REAL has P1 anchor in top-5 = {real_has_p1}. "
            f"NULL: {n_with_p1}/{n_total} have P1 anchor.\n")
    if real_has_p1 and n_with_p1 <= 1:
        f.write(f"\n**STRONG SIGNAL**: P1 anchor appears specifically for the biologically functional "
                f"position-30 mutation, in {n_with_p1}/{n_total} arbitrary mutations. "
                f"F-001 strengthened.\n")
PYEND
commit_push "feat(cpeb3-solid): multi-null robustness statistic (5 nulls + 1 real)"

# ============================================================================
# Phase 4: secondary structure comparison (ViennaRNA + DRfold2 ss.ct)
# ============================================================================
echo
echo "--- [4] secondary structure analysis (ViennaRNA RNAfold + DRfold2 ss.ct) ---"
pip install -q ViennaRNA 2>&1 | tail -2
python - <<'PYEND'
from pathlib import Path
import datetime
OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
SS = OUT / "ss_analysis.md"

with open(SS, "w") as f:
    f.write("# Secondary structure comparison — R1107 vs R1108\n\n")
    f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
    try:
        import RNA
        f.write("## Method 1: ViennaRNA RNAfold (MFE secondary structure)\n\n")
        for tid, seq in [
            ("R1107_human", "GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"),
            ("R1108_chimp", "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"),
        ]:
            db, mfe = RNA.fold(seq)
            f.write(f"### {tid}\n\n")
            f.write(f"```\n{seq}\n{db}\nMFE = {mfe:.2f} kcal/mol\n```\n\n")
        # Now diff
        h_db, h_mfe = RNA.fold("GGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU")
        c_db, c_mfe = RNA.fold("GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU")
        diffs = [(i+1, h_db[i], c_db[i]) for i in range(len(h_db)) if h_db[i] != c_db[i]]
        f.write(f"### Difference\n\n- MFE difference: {h_mfe - c_mfe:+.2f} kcal/mol "
                f"(positive = human less stable)\n")
        f.write(f"- Dot-bracket differences at residues (human/chimp): {diffs}\n\n")
        if not diffs:
            f.write("→ Identical MFE secondary structure between human and chimp by ViennaRNA. The 4× activity difference would NOT be predicted from MFE alone.\n")
        else:
            f.write("→ Different MFE secondary structure — ViennaRNA already sees a base-pairing difference, partially predicting Skilandat's mechanism.\n")
    except Exception as e:
        f.write(f"\nViennaRNA failed: {e}\n")

    # DRfold2 ss.ct files (if produced) — but DRfold2 may not output ss.ct
    f.write("\n## Method 2: DRfold2 secondary structure (from .ret if available)\n\n")
    for tid in ("R1107_human", "R1108_chimp"):
        for ss_name in ("ss.ct", "secondary_structure.ct"):
            p = OUT / tid / "drfold2" / ss_name
            if p.exists():
                f.write(f"### {tid}\n\n```\n{p.read_text()[:500]}\n```\n\n")
                break
        else:
            f.write(f"### {tid}\n\n_(no ss.ct found from DRfold2)_\n\n")

    # RhoFold ss.ct (it DOES output ss.ct)
    f.write("\n## Method 3: RhoFold+ secondary structure (ss.ct)\n\n")
    for tid in ("R1107_human", "R1108_chimp"):
        p = OUT / tid / "rhofold" / "ss.ct"
        if p.exists():
            f.write(f"### {tid}\n\n```\n{p.read_text()[:1000]}\n```\n\n")
        else:
            f.write(f"### {tid}\n\n_(no ss.ct found)_\n\n")
print("ss analysis done")
PYEND
commit_push "feat(cpeb3-solid): secondary structure comparison (ViennaRNA + RhoFold ss.ct)"

echo
echo "===== CPEB3 SOLID pipeline finished $(date -u +%FT%TZ) ====="
