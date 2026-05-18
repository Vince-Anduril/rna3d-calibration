#!/usr/bin/env bash
# Experiment #2 — n=20 random NULL controls + permutation test.
# Generates 20 single-nt mutations of R1108 at uniformly random positions
# (excluding pos 30 which is the real biological mutation), runs DRfold2 on
# each, then computes proper frequency statistics.
#
# Runtime: ~60 min on RTX 5090 (~3 min per DRfold2 inference × 20). The
# script is idempotent: it skips any null_rand_NN that already has folds.

set +e
cd /workspace/rna3d
source venv/bin/activate

LOG=/workspace/rna3d/cpeb3_n20.log
OUT=/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused
exec > >(tee -a $LOG) 2>&1
echo "===== CPEB3 n=20 random NULLs started $(date -u +%FT%TZ) ====="

commit_push() {
  cd /workspace/rna3d
  git add .research/30_experiments/runs/cpeb3_focused/null_rand_*/input.fasta \
          .research/30_experiments/runs/cpeb3_focused/null_rand_*/drfold2/folds/opt_0_*.pdb \
          .research/30_experiments/runs/cpeb3_focused/SUMMARY.md \
          .research/30_experiments/runs/cpeb3_focused/n20_results.csv \
          .research/30_experiments/runs/cpeb3_focused/n20_results.md \
          .research/40_findings/FINDINGS.md \
          2>/dev/null
  if ! git diff --staged --quiet; then
    git commit -q -m "$1" 2>&1 | tail -1
    git push origin main 2>&1 | tail -2
  fi
}

# ============ Phase A: generate 20 random single-nt variants ============
echo
echo "--- [A] generate 20 random single-nt variants of R1108 ---"
python <<'PYEND'
import random, json
from pathlib import Path

R1108 = "GGGGGCCACAGCAGAAGCGUUCACGUCGCGGCCCCUGUCAGCCAUUGCACUCCGGCUGCGAAUUCUGCU"
assert len(R1108) == 69
BASES = ["A", "U", "G", "C"]
random.seed(42)

OUT = Path("/workspace/rna3d/.research/30_experiments/runs/cpeb3_focused")

# Exclude pos 30 (the real biological mutation) from random sampling.
positions = list(range(1, 70))
positions.remove(30)
# Sample 20 unique positions
sampled = random.sample(positions, 20)
sampled.sort()

manifest = []
for pos in sampled:
    orig = R1108[pos - 1]
    alt_choices = [b for b in BASES if b != orig]
    alt = random.choice(alt_choices)
    mutated = R1108[:pos-1] + alt + R1108[pos:]
    assert len(mutated) == 69
    tid = f"null_rand_{pos:02d}"
    td = OUT / tid
    td.mkdir(parents=True, exist_ok=True)
    td.joinpath("input.fasta").write_text(f">{tid}\n{mutated}\n")
    manifest.append({"target_id": tid, "position": pos, "original": orig, "alt": alt, "sequence": mutated})

with open(OUT / "n20_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)

print(f"Generated 20 random NULL targets at positions: {sampled}")
PYEND
commit_push "feat(cpeb3-n20): generate 20 random NULL controls (seed=42)"

# ============ Phase B: run DRfold2 on each (skip if folds exist) ============
echo
echo "--- [B] DRfold2 inference on 20 NULLs (skip if folds exist) ---"
for tid_dir in "$OUT"/null_rand_*; do
  tid=$(basename "$tid_dir")
  mkdir -p "$tid_dir/drfold2"
  if ls "$tid_dir/drfold2/folds/opt_0_"*.pdb >/dev/null 2>&1; then
    echo "  [skip] $tid (folds present)"
    continue
  fi
  echo ">>> running DRfold2 on $tid <<<"
  (cd /workspace/rna3d/models/DRfold2 && python DRfold_infer.py "$tid_dir/input.fasta" "$tid_dir/drfold2" 2>&1 | tail -4)
  # Commit + push after each (so failure mid-loop doesn't lose progress)
  commit_push "feat(cpeb3-n20): DRfold2 inference for $tid"
done

# ============ Phase C: aggregate stats + permutation test ============
echo
echo "--- [C] aggregate stats + permutation test ---"
python <<'PYEND'
import numpy as np, json, csv, datetime
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
if chimp is None or human is None:
    print("Missing R1107_human or R1108_chimp"); raise SystemExit(1)

# Real case
real_delta = kabsch_per_res(human, chimp)
real_top5 = set(np.argsort(real_delta)[-5:][::-1].tolist())  # 0-indexed
real_top5_1based = sorted(i+1 for i in real_top5)
real_has_anchor = (8 in real_top5) and (59 in real_top5)  # 0-indexed for residues 9 and 60

# All 20 random NULLs
rows = [{"target_id": "R1107_human", "position": 30, "kind": "REAL",
         "mean_delta": float(real_delta.mean()),
         "max_delta": float(real_delta.max()),
         "delta_at_9": float(real_delta[8]),
         "delta_at_60": float(real_delta[59]),
         "top5_1based": ",".join(str(i+1) for i in sorted(real_top5)),
         "pos9_in_top5": 8 in real_top5,
         "pos60_in_top5": 59 in real_top5,
         "both_anchors_in_top5": real_has_anchor,
         "either_anchor_in_top5": (8 in real_top5) or (59 in real_top5)}]

with open(OUT / "n20_manifest.json") as f:
    manifest = json.load(f)

for entry in manifest:
    tid = entry["target_id"]; pos = entry["position"]
    p = best(tid)
    if p is None:
        rows.append({"target_id": tid, "position": pos, "kind": "NULL",
                     "mean_delta": None, "max_delta": None,
                     "delta_at_9": None, "delta_at_60": None,
                     "top5_1based": "MISSING", "pos9_in_top5": False,
                     "pos60_in_top5": False, "both_anchors_in_top5": False,
                     "either_anchor_in_top5": False})
        continue
    d = kabsch_per_res(p, chimp)
    top5 = set(np.argsort(d)[-5:][::-1].tolist())
    rows.append({"target_id": tid, "position": pos, "kind": "NULL",
                 "mean_delta": float(d.mean()), "max_delta": float(d.max()),
                 "delta_at_9": float(d[8]), "delta_at_60": float(d[59]),
                 "top5_1based": ",".join(str(i+1) for i in sorted(top5)),
                 "pos9_in_top5": 8 in top5, "pos60_in_top5": 59 in top5,
                 "both_anchors_in_top5": (8 in top5) and (59 in top5),
                 "either_anchor_in_top5": (8 in top5) or (59 in top5)})

# CSV
with open(OUT / "n20_results.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

nulls = [r for r in rows if r["kind"] == "NULL" and r["top5_1based"] != "MISSING"]
n_nulls = len(nulls)
n_either = sum(1 for r in nulls if r["either_anchor_in_top5"])
n_both = sum(1 for r in nulls if r["both_anchors_in_top5"])

# Empirical p-value: under H0 (random structural prior), what fraction of NULL
# mutations produce the "both 9 AND 60 in top-5" pattern?
p_emp = (n_both + 1) / (n_nulls + 1)  # add-1 smoothed

# Theoretical baseline: probability of any 2 specific residues both being in
# top-5 of 69 if top-5 were uniformly random = C(67, 3) / C(69, 5) ≈ 4 / 4641 ≈ 0.0086
import math
theoretical = math.comb(67, 3) / math.comb(69, 5)

print(f"\n=== n={n_nulls} random NULL statistics ===")
print(f"REAL case has both anchors in top-5: {rows[0]['both_anchors_in_top5']}")
print(f"NULL cases with pos 9 in top-5:      {sum(1 for r in nulls if r['pos9_in_top5'])} / {n_nulls}")
print(f"NULL cases with pos 60 in top-5:     {sum(1 for r in nulls if r['pos60_in_top5'])} / {n_nulls}")
print(f"NULL cases with BOTH 9 AND 60:       {n_both} / {n_nulls}")
print(f"NULL cases with EITHER 9 OR 60:      {n_either} / {n_nulls}")
print(f"Empirical add-1-smoothed p-value (both anchors): {p_emp:.4f}")
print(f"Theoretical baseline (uniform top-5): {theoretical:.4f}")

with open(OUT / "n20_results.md", "w") as f:
    f.write(f"# Experiment #2 — n={n_nulls} random NULL controls\n\n")
    f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
    f.write(f"## Headline\n\n")
    f.write(f"- **REAL case** (pos 30, biological): both anchors (9, 60) in top-5 = `{rows[0]['both_anchors_in_top5']}`\n")
    f.write(f"- **NULL controls** (n={n_nulls}, random single-nt at uniformly sampled positions):\n")
    f.write(f"  - with pos 9 in top-5:  **{sum(1 for r in nulls if r['pos9_in_top5'])}/{n_nulls}**\n")
    f.write(f"  - with pos 60 in top-5: **{sum(1 for r in nulls if r['pos60_in_top5'])}/{n_nulls}**\n")
    f.write(f"  - with **BOTH** 9 AND 60 in top-5: **{n_both}/{n_nulls}**\n")
    f.write(f"  - with EITHER 9 OR 60 in top-5:  **{n_either}/{n_nulls}**\n\n")
    f.write(f"## Statistical interpretation\n\n")
    f.write(f"- Add-1-smoothed empirical p-value for ``both anchors in top-5'' under H0 (random prior): **{p_emp:.4f}**\n")
    f.write(f"- Theoretical baseline (uniform random top-5 of 5 out of 69): {theoretical:.4f}\n\n")
    if n_both == 0:
        f.write(f"**Outcome:** No random NULL mutation reproduces the ``both anchors'' pattern. The empirical p-value upper bound is {p_emp:.3f} ({1}/(n+1) with add-1 smoothing). The biological mutation is statistically distinguishable from the random-mutation null at this n.\n")
    else:
        f.write(f"**Outcome:** {n_both}/{n_nulls} random NULLs DO reproduce the both-anchors pattern. Empirical p ≈ {p_emp:.3f}. The signal is NOT specific at the resolution of this test.\n")
    f.write(f"\n## Detailed table\n\n")
    f.write("| target | pos | kind | mean Δ | max Δ | Δ@9 | Δ@60 | top-5 | both anchors? |\n")
    f.write("|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        if r["mean_delta"] is None:
            f.write(f"| {r['target_id']} | {r['position']} | {r['kind']} | - | - | - | - | MISSING | - |\n")
        else:
            mark = "**YES**" if r["both_anchors_in_top5"] else "no"
            f.write(f"| {r['target_id']} | {r['position']} | {r['kind']} | {r['mean_delta']:.2f} | {r['max_delta']:.2f} | {r['delta_at_9']:.2f} | {r['delta_at_60']:.2f} | {{{r['top5_1based']}}} | {mark} |\n")
PYEND
commit_push "feat(cpeb3-n20): final aggregation + permutation-style p-value"

echo
echo "===== CPEB3 n=20 pipeline finished $(date -u +%FT%TZ) ====="
