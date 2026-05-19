#!/usr/bin/env bash
# Extension B — full 69-position mutational scan of R1108.
#
# For EACH of the 69 positions, generate one single-nt mutation (random alt
# base != original), run DRfold2 inference, then aggregate into a heatmap
# of "P1 anchor activation" per position.
#
# Runtime: 69 × ~3 min = ~3.5 hours on RTX 5090. Idempotent (skips done
# positions). Auto-commits/pushes per inference + after final aggregate.

set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/cpeb3_fullscan.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
exec > >(tee -a $LOG) 2>&1
echo "===== CPEB3 full 69-position scan started $(date -u +%FT%TZ) ====="

commit_push() {
  cd /workspace/rna3d
  git add .research/30_experiments/runs/cpeb3_focused/scan_*/input.fasta \
          .research/30_experiments/runs/cpeb3_focused/scan_*/drfold2/folds/opt_0_*.pdb \
          .research/30_experiments/runs/cpeb3_focused/full_scan_results.csv \
          .research/30_experiments/runs/cpeb3_focused/full_scan_summary.md 2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -1
    git push origin main 2>&1 | tail -2
  fi
}

# Phase A: generate 69 single-nt variants (one per position, deterministic alt)
echo
echo "--- [A] generate 69 single-nt variants (deterministic alt = next base AUGCA) ---"
python <<'PYEND'
from pathlib import Path
R1108 = "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"
# Cyclic-next-base rule: A→U, U→G, G→C, C→A. Deterministic, always changes the base.
NEXT = {"A": "U", "U": "G", "G": "C", "C": "A"}
OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")
for pos in range(1, 70):
    tid = f"scan_pos{pos:02d}"
    td = OUT / tid; td.mkdir(parents=True, exist_ok=True)
    orig = R1108[pos-1]
    alt = NEXT[orig]
    seq = R1108[:pos-1] + alt + R1108[pos:]
    assert len(seq) == 69 and seq != R1108
    (td / "input.fasta").write_text(f">{tid}\n{seq}\n")
print("69 input.fasta files written")
PYEND
commit_push "feat(fullscan): generate 69 deterministic single-nt scan variants"

# Phase B: inference for each (skip if folds present)
echo
echo "--- [B] DRfold2 inference on 69 positions ---"
for pos in $(seq 1 69); do
  tid="scan_pos$(printf '%02d' $pos)"
  td="$OUT/$tid"
  mkdir -p "$td/drfold2"
  if ls "$td/drfold2/folds/opt_0_"*.pdb >/dev/null 2>&1; then
    continue
  fi
  echo ">>> $tid <<<"
  (cd /workspace/rna3d/models/DRfold2 && python DRfold_infer.py "$td/input.fasta" "$td/drfold2" 2>&1 | tail -3)
  # Commit after every 5 inferences to keep GitHub fresh
  if [ $((pos % 5)) -eq 0 ]; then
    commit_push "feat(fullscan): inference progress through pos $pos"
  fi
done
commit_push "feat(fullscan): all 69 inferences complete"

# Phase C: aggregate
echo
echo "--- [C] aggregate full scan results ---"
python <<'PYEND'
import numpy as np, csv, datetime
from pathlib import Path
import biotite.structure.io as bsio

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")

def c1(p):
    if not p.exists(): return None
    s = bsio.load_structure(str(p))
    if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
    m = s.atom_name == "C1'"
    return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None

def best(tid):
    g = sorted((OUT / tid / "drfold2" / "folds").glob("opt_0_*.pdb"))
    return c1(g[0]) if g else None

def kabsch_per_res(p, q):
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al-q, axis=-1)

chimp = best("R1108_chimp")
human = best("R1107_human")
if chimp is None: raise SystemExit("missing R1108_chimp")

# Real case as reference
real_d = kabsch_per_res(human, chimp) if human is not None else None

rows = []
for pos in range(1, 70):
    tid = f"scan_pos{pos:02d}"
    p = best(tid)
    if p is None:
        continue
    d = kabsch_per_res(p, chimp)
    top5 = set(int(i)+1 for i in np.argsort(d)[-5:])
    rows.append({
        "position": pos, "target_id": tid,
        "mean_delta": float(d.mean()),
        "max_delta": float(d.max()),
        "delta_at_mut": float(d[pos-1]),
        "delta_at_9": float(d[8]),
        "delta_at_60": float(d[59]),
        "pos9_in_top5": 9 in top5,
        "pos60_in_top5": 60 in top5,
        "both_anchors_in_top5": (9 in top5) and (60 in top5),
        "either_anchor_in_top5": (9 in top5) or (60 in top5),
        "top5_1based": ",".join(str(i) for i in sorted(top5)),
    })

# CSV
with open(OUT / "full_scan_results.csv", "w", newline="") as f:
    if rows:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

# Summary
n_total = len(rows)
n_both = sum(1 for r in rows if r["both_anchors_in_top5"])
n_either = sum(1 for r in rows if r["either_anchor_in_top5"])
n_pos9 = sum(1 for r in rows if r["pos9_in_top5"])
n_pos60 = sum(1 for r in rows if r["pos60_in_top5"])
positions_both = [r["position"] for r in rows if r["both_anchors_in_top5"]]
positions_either = [r["position"] for r in rows if r["either_anchor_in_top5"]]

with open(OUT / "full_scan_summary.md", "w") as f:
    f.write(f"# Full 69-position mutational scan — DRfold2\n\n")
    f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
    f.write(f"## Headline\n\n")
    f.write(f"- Positions scanned: {n_total} / 69\n")
    f.write(f"- Positions whose mutation puts BOTH anchors (9 and 60) in top-5: **{n_both} / {n_total}**\n")
    f.write(f"- Positions whose mutation puts EITHER anchor in top-5: **{n_either} / {n_total}**\n")
    f.write(f"- Positions whose mutation puts pos 9 in top-5: **{n_pos9} / {n_total}**\n")
    f.write(f"- Positions whose mutation puts pos 60 in top-5: **{n_pos60} / {n_total}**\n\n")
    if positions_both:
        f.write(f"## Positions that activate BOTH anchors\n\n{positions_both}\n\n")
    f.write(f"## REAL biological mutation (pos 30 A→G, R1107)\n\n")
    if real_d is not None:
        top5_real = set(int(i)+1 for i in np.argsort(real_d)[-5:])
        f.write(f"- top-5: {{{','.join(str(i) for i in sorted(top5_real))}}}\n")
        f.write(f"- both anchors in top-5: {9 in top5_real and 60 in top5_real}\n\n")
    f.write(f"## Interpretation\n\n")
    if n_both == 1:
        f.write(f"**EXCELLENT SPECIFICITY.** Of 69 positions, only ONE mutation activates both anchors. If that position is 30 (the real biological mutation), the F-001 cascade is statistically distinguished from background at p < 1/69 ≈ 0.014.\n")
    elif n_both <= 5:
        f.write(f"**STRONG SPECIFICITY.** {n_both}/69 positions activate both anchors. Report the empirical p-value as {n_both}/69 ≈ {n_both/69:.3f}.\n")
    else:
        f.write(f"**WEAK SPECIFICITY.** {n_both}/69 positions activate both anchors. The 'both anchor' pattern is a model attractor more than a mutation-specific response.\n")
    f.write(f"\n## Per-position detail (sorted by delta_at_9 + delta_at_60)\n\n")
    f.write("| pos | mut Δ | mean Δ | max Δ | Δ@9 | Δ@60 | both anchors? | top-5 |\n")
    f.write("|---|---|---|---|---|---|---|---|\n")
    for r in sorted(rows, key=lambda r: -(r["delta_at_9"] + r["delta_at_60"])):
        mark = "**YES**" if r["both_anchors_in_top5"] else "no"
        f.write(f"| {r['position']} | {r['delta_at_mut']:.2f} | {r['mean_delta']:.2f} | {r['max_delta']:.2f} | {r['delta_at_9']:.2f} | {r['delta_at_60']:.2f} | {mark} | {{{r['top5_1based']}}} |\n")
PYEND
commit_push "feat(fullscan): aggregate 69-position scan — full statistical map"

echo
echo "===== full scan finished $(date -u +%FT%TZ) ====="
