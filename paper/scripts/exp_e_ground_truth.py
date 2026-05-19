"""Experiment E — ground-truth comparison vs PDB 7QR3 chain C.

7QR3 is the X-ray crystal of the chimpanzee CPEB3 ribozyme (Przytula-Mally
et al. 2022), 2.18 Å. Chains C and D are the two RNA monomers in the
crystallographic dimer; both are 69 nt. We use chain C as the ground
truth for R1108.

We compute Kabsch C1' RMSD of:
  - DRfold2 R1107_human vs 7QR3:C  (cross-species, 1 nt mismatch — best
    available proxy for human ground truth)
  - DRfold2 R1108_chimp vs 7QR3:C  (matched sequence)
  - AF3 R1107 vs 7QR3:C
  - AF3 R1108 vs 7QR3:C

And we look at where (per-residue) each model diverges most from the
crystal.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
CRYSTAL = ROOT / "data/crystal/7QR3.cif"
AF3_CIF = ROOT / "data/af3/raw/fold_2026_05_19_09_18_model_0.cif"
DRF_RUN = ROOT / ".research/30_experiments/runs/cpeb3_focused"
OUT = DRF_RUN


def parse_cif_c1(path: Path, chain_filter: str | None = None):
    """Return dict {chain: ndarray[L,3]} of C1' coords, ordered by residue."""
    lines = path.read_text().splitlines()
    cols, data_start = [], None
    in_header = False
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s); in_header = True
        elif in_header and not s.startswith("_atom_site."):
            data_start = i; break
    by_chain: dict[str, list[tuple[int, np.ndarray]]] = {}
    for ln in lines[data_start:]:
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        parts = ln.split()
        if len(parts) < len(cols): continue
        d = dict(zip(cols, parts))
        aid = d["_atom_site.label_atom_id"].strip('"')
        if aid != "C1'": continue
        ch = d["_atom_site.label_asym_id"]
        if chain_filter and ch != chain_filter: continue
        try:
            ri = int(d["_atom_site.label_seq_id"])
            xyz = np.array([float(d["_atom_site.Cartn_x"]),
                            float(d["_atom_site.Cartn_y"]),
                            float(d["_atom_site.Cartn_z"])], dtype=np.float32)
        except (KeyError, ValueError):
            continue
        by_chain.setdefault(ch, []).append((ri, xyz))
    out = {}
    for c, lst in by_chain.items():
        lst.sort()
        out[c] = np.stack([x for _, x in lst], axis=0)
    return out


def load_drfold_c1(target):
    folds = sorted((DRF_RUN / target / "drfold2" / "folds").glob("opt_0_*.pdb"))
    if not folds:
        return None
    coords = []
    for ln in folds[0].read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "C1'":
            try:
                coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
    return np.array(coords, dtype=np.float32) if coords else None


def kabsch_per_res(p, q):
    if p is None or q is None or p.shape != q.shape:
        return None
    pc, qc = p.mean(0), q.mean(0)
    pp, qq = p - pc, q - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al - q, axis=-1)


def rmsd(p, q):
    dpr = kabsch_per_res(p, q)
    return float(np.sqrt((dpr * dpr).mean())) if dpr is not None else float("nan")


def main():
    print("Loading 7QR3 crystal (chain C, R1108 chimp ground truth)...")
    cry = parse_cif_c1(CRYSTAL, chain_filter=None)
    crystal_C = cry.get("C")
    crystal_D = cry.get("D")
    print(f"  chain C: {crystal_C.shape if crystal_C is not None else None}")
    print(f"  chain D: {crystal_D.shape if crystal_D is not None else None}")
    # Crystal C ↔ D sanity (should be near 0 if same conformation)
    if crystal_C is not None and crystal_D is not None:
        print(f"  crystal C ↔ D RMSD = {rmsd(crystal_C, crystal_D):.2f} Å (sanity: same molecule in dimer)")

    print("\nLoading AF3 multi-chain predictions...")
    af3 = parse_cif_c1(AF3_CIF)
    af3_R1107 = af3.get("A"); af3_R1108 = af3.get("B")
    print(f"  AF3 R1107 (chain A): {af3_R1107.shape}")
    print(f"  AF3 R1108 (chain B): {af3_R1108.shape}")

    print("\nLoading DRfold2 predictions...")
    drf_R1107 = load_drfold_c1("R1107_human")
    drf_R1108 = load_drfold_c1("R1108_chimp")
    print(f"  DRfold2 R1107: {drf_R1107.shape if drf_R1107 is not None else None}")
    print(f"  DRfold2 R1108: {drf_R1108.shape if drf_R1108 is not None else None}")

    print("\n" + "=" * 60)
    print("GROUND-TRUTH RMSD vs 7QR3:C (chimpanzee CPEB3 X-ray crystal)")
    print("=" * 60)
    # All 5 AF3 seeds
    af3_per_seed = {}
    for seed in range(5):
        cif_s = AF3_CIF.parent / f"fold_2026_05_19_09_18_model_{seed}.cif"
        d = parse_cif_c1(cif_s)
        af3_per_seed[seed] = {"A": d.get("A"), "B": d.get("B")}
    print(f"\nAF3 R1108 vs crystal (all 5 seeds):")
    af3_R1108_rmsds = []
    for s in range(5):
        r = rmsd(af3_per_seed[s]["B"], crystal_C)
        af3_R1108_rmsds.append(r)
        print(f"  seed {s}: {r:.2f} Å")
    best_af3 = min(af3_R1108_rmsds)
    worst_af3 = max(af3_R1108_rmsds)
    print(f"  best-of-5: {best_af3:.2f} Å, worst: {worst_af3:.2f} Å")

    results = []
    for name, pred in [
        ("DRfold2 R1108 (matched sequence)", drf_R1108),
        ("AF3 R1108 default-seed (matched)", af3_R1108),
        ("AF3 R1108 best-of-5    (matched)", af3_per_seed[int(np.argmin(af3_R1108_rmsds))]["B"]),
        ("DRfold2 R1107 (cross-species)",    drf_R1107),
        ("AF3 R1107 default-seed (cross-sp)", af3_R1107),
    ]:
        r = rmsd(pred, crystal_C)
        results.append((name, r))
        print(f"  {name}: {r:.2f} Å")

    # Per-residue divergence DRfold2 R1108 vs crystal — to spot where it agrees/disagrees
    pr_drf = kabsch_per_res(drf_R1108, crystal_C)
    pr_af3 = kabsch_per_res(af3_R1108, crystal_C)
    print("\n  Per-residue C1' deviation R1108 vs 7QR3:C")
    print(f"  DRfold2: mean={pr_drf.mean():.2f}, max={pr_drf.max():.2f} Å (at res {int(np.argmax(pr_drf))+1})")
    print(f"  AF3:     mean={pr_af3.mean():.2f}, max={pr_af3.max():.2f} Å (at res {int(np.argmax(pr_af3))+1})")
    print()
    print(f"  Cascade-residue distances (R1108 vs 7QR3:C):")
    for r in [9, 22, 24, 30, 51, 60]:
        print(f"    res {r}: DRfold2={pr_drf[r-1]:.2f} Å,  AF3={pr_af3[r-1]:.2f} Å")

    # Write report
    import datetime
    md = OUT / "ground_truth_results.md"
    with open(md, "w") as f:
        f.write("# Experiment E — predictions vs 7QR3 crystal (ground truth)\n\n")
        f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
        f.write("Crystal: **PDB 7QR3 chain C** — chimpanzee CPEB3 ribozyme, X-ray 2.18 Å (Przytula-Mally 2022).\n\n")
        f.write("Crystal contains 2 chains (C, D) of the same RNA in a crystallographic dimer; ")
        f.write(f"chain C↔chain D Kabsch RMSD = {rmsd(crystal_C, crystal_D):.2f} Å.\n\n")
        f.write("## Headline: Kabsch C1' RMSD vs 7QR3:C\n\n")
        f.write("| Predictor | RMSD vs crystal (Å) |\n|---|---|\n")
        for name, r in results:
            f.write(f"| {name} | **{r:.2f}** |\n")
        f.write("\n### AF3 per-seed (all 5)\n\n")
        f.write("| seed | RMSD (Å) | ranking_score |\n|---|---|---|\n")
        import json as _j
        for s in range(5):
            js = _j.loads((AF3_CIF.parent / f"fold_2026_05_19_09_18_summary_confidences_{s}.json").read_text())
            f.write(f"| {s} | {af3_R1108_rmsds[s]:.2f} | {js.get('ranking_score', 0)} |\n")
        f.write("\n## Per-residue deviation vs crystal (R1108, matched sequence)\n\n")
        f.write("| residue | DRfold2 Δ (Å) | AF3 Δ (Å) |\n|---|---|---|\n")
        for r in range(1, 70):
            f.write(f"| {r} | {pr_drf[r-1]:.2f} | {pr_af3[r-1]:.2f} |\n")
        f.write("\n## Interpretation\n\n")
        drf_b = next(r for n, r in results if "DRfold2 R1108" in n)
        af3_b = next(r for n, r in results if "AF3     R1108" in n)
        gap = af3_b - drf_b
        if gap > 5:
            f.write(f"**On the matched (chimpanzee) sequence, DRfold2 is {gap:.1f} Å closer to the experimental crystal than AlphaFold 3** (DRfold2 {drf_b:.2f} Å vs AF3 {af3_b:.2f} Å). On this short orphan-like ribozyme, the single-sequence transformer beats AF3 against ground truth.\n")
        elif gap < -5:
            f.write(f"**AF3 is {-gap:.1f} Å closer to the crystal than DRfold2** on the matched sequence, despite AF3's low self-reported pLDDT (~25). pLDDT is miscalibrated for this RNA class.\n")
        else:
            f.write(f"DRfold2 and AF3 are within {abs(gap):.1f} Å of each other vs ground truth. Both reasonable / both struggling.\n")
    print(f"\nWritten: {md}")


if __name__ == "__main__":
    main()
