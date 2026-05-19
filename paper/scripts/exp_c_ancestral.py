"""Extension C — Ancestral CPEB3 resurrection probe.

WHAT THIS DOES (run on pod after manual data fetch):
    1. Read ancestral sequences from data/ancestral/sequences.fasta
       (user must fetch them first — see exp_c_ancestral_fetch.md).
    2. Run DRfold2 inference on each ancestor (skip if folds exist).
    3. Compute Kabsch C1' RMSD pair-wise between every ancestor and
       (a) R1107 human, (b) R1108 chimpanzee, (c) the consensus ancestral.
    4. KEY ANALYSIS: does predicted structural divergence track the
       measured experimental cleavage activity reported in Bendixsen 2021?
    5. Output: ancestral_results.csv + ancestral_results.md.

USAGE (on pod):
    python paper/scripts/exp_c_ancestral.py

REQUIRED DATA:
    /workspace/rna3d/data/ancestral/sequences.fasta — multi-FASTA with
    one record per ancestor; record IDs should match the Bendixsen 2021
    node labels (e.g. "ancestral_mammalian", "ancestral_primate", ...).

    Optionally /workspace/rna3d/data/ancestral/activity.csv with columns
    (ancestor_id, fraction_cleaved) for the experimental activity rates.
"""
from __future__ import annotations
import sys, csv, datetime
from pathlib import Path
import numpy as np

ROOT = Path("/workspace/rna3d")
DATA = ROOT / "data/ancestral"
OUT = ROOT / ".research/30_experiments/runs/cpeb3_focused/ancestral"
OUT.mkdir(parents=True, exist_ok=True)


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    if not path.exists(): return seqs
    cur = None; buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur and buf:
                seqs[cur] = "".join(buf).replace("T", "U").upper()
            cur = line[1:].strip().split()[0]; buf = []
        elif line.strip():
            buf.append(line.strip())
    if cur and buf:
        seqs[cur] = "".join(buf).replace("T", "U").upper()
    return seqs


def run_drfold2(target_id: str, seq: str):
    td = OUT / target_id
    td.mkdir(parents=True, exist_ok=True)
    fasta = td / "input.fasta"
    fasta.write_text(f">{target_id}\n{seq}\n")
    out_dir = td / "drfold2"
    out_dir.mkdir(exist_ok=True)
    folds = list(out_dir.glob("folds/opt_0_*.pdb"))
    if folds:
        print(f"  [skip] {target_id}")
        return folds[0]
    import subprocess
    print(f"  [run] {target_id} ({len(seq)} nt)")
    subprocess.run(
        ["python", "/workspace/rna3d/models/DRfold2/DRfold_infer.py",
         str(fasta), str(out_dir)],
        cwd="/workspace/rna3d/models/DRfold2",
        check=False,
    )
    folds = list(out_dir.glob("folds/opt_0_*.pdb"))
    return folds[0] if folds else None


def c1(p):
    if p is None or not p.exists(): return None
    import biotite.structure.io as bsio
    s = bsio.load_structure(str(p))
    if hasattr(s, "stack_depth") and s.stack_depth() > 1: s = s[0]
    m = s.atom_name == "C1'"
    return np.asarray(s.coord[m], dtype=np.float32) if m.any() else None


def kabsch_rmsd(p, q):
    if p is None or q is None or p.shape != q.shape or p.shape[0] < 3: return float("nan")
    pc, qc = p.mean(0), q.mean(0); pp, qq = p-pc, q-qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0,1.0,d]) @ U.T
    al = (R @ pp.T).T + qc
    return float(np.sqrt(((al-q)**2).sum(-1).mean()))


def main():
    fasta_path = DATA / "sequences.fasta"
    if not fasta_path.exists():
        print(f"NO DATA at {fasta_path}.")
        print("Please follow exp_c_ancestral_fetch.md to obtain the sequences from")
        print("https://gitlab.com/devinbendixsen/cpeb3_phylo_fls (Bendixsen 2021 supp).")
        return

    ancestors = parse_fasta(fasta_path)
    print(f"Loaded {len(ancestors)} ancestral sequences:")
    for name, seq in ancestors.items():
        print(f"  {name}: {len(seq)} nt")

    # Activity table (optional)
    activity: dict[str, float] = {}
    act_path = DATA / "activity.csv"
    if act_path.exists():
        with open(act_path) as f:
            for row in csv.DictReader(f):
                activity[row["ancestor_id"]] = float(row.get("fraction_cleaved", "nan"))
        print(f"Loaded {len(activity)} activity values")

    # Run DRfold2 on each
    pdb_paths: dict[str, Path | None] = {}
    for name, seq in ancestors.items():
        pdb_paths[name] = run_drfold2(name, seq)

    # Also load R1107 and R1108 references (already on pod)
    main_run = ROOT / ".research/30_experiments/runs/cpeb3_focused"
    def best(tid):
        g = sorted((main_run / tid / "drfold2" / "folds").glob("opt_0_*.pdb"))
        return c1(g[0]) if g else None
    r1107 = best("R1107_human"); r1108 = best("R1108_chimp")

    # Pair-wise comparison
    coords = {name: c1(p) for name, p in pdb_paths.items() if p is not None}
    coords["R1107_human"] = r1107
    coords["R1108_chimp"] = r1108

    rows = []
    for name, c in coords.items():
        rows.append({
            "id": name,
            "n_residues": int(c.shape[0]) if c is not None else 0,
            "rmsd_vs_R1107": kabsch_rmsd(c, r1107) if c is not None else float("nan"),
            "rmsd_vs_R1108": kabsch_rmsd(c, r1108) if c is not None else float("nan"),
            "activity_fraction_cleaved": activity.get(name, float("nan")),
        })

    csv_path = OUT.parent / "ancestral_results.csv"
    with open(csv_path, "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    md_path = OUT.parent / "ancestral_results.md"
    with open(md_path, "w") as f:
        f.write(f"# Ancestral CPEB3 resurrection — Extension C\n\n")
        f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
        f.write("| id | n_res | RMSD vs R1107 | RMSD vs R1108 | activity |\n")
        f.write("|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['id']} | {r['n_residues']} | {r['rmsd_vs_R1107']:.3f} | {r['rmsd_vs_R1108']:.3f} | {r['activity_fraction_cleaved']:.3f} |\n")
        f.write("\n## Hypothesis\n\nDoes the predicted structural divergence from an ancestor correlate with the measured activity gap? If yes, structural predictions track function across evolutionary time.\n")

    print(f"\nWritten: {csv_path}\nWritten: {md_path}")


if __name__ == "__main__":
    main()
