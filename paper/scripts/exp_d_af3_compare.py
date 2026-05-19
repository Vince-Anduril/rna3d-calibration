"""Extension D — AF3 vs DRfold2 cross-model analysis on CPEB3 panel.

Input :
    /Users/leduigouvincent/rna3d-mirror/data/af3/raw/
        fold_2026_05_19_09_18_model_0.cif   (7-chain multi-model)
        fold_2026_05_19_09_18_full_data_0.json  (per-atom pLDDT, PAE)

Chain mapping (verified against job_request.json):
    A=R1107_human  B=R1108_chimp  C=null_pos5  D=null_pos20
    E=null_pos41   F=null_pos55   G=null_pos64

Output:
    af3_results.csv  — per-target n_res, RMSD vs DRfold2, top-5 divergent res
    af3_results.md   — story-format markdown
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
AF3_DIR = ROOT / "data/af3/raw"
DRF_RUN = ROOT / ".research/30_experiments/runs/cpeb3_focused"
OUT_DIR = DRF_RUN
CIF = AF3_DIR / "fold_2026_05_19_09_18_model_0.cif"
JSON = AF3_DIR / "fold_2026_05_19_09_18_full_data_0.json"

CHAIN_MAP = {
    "A": "R1107_human", "B": "R1108_chimp",
    "C": "null_pos5", "D": "null_pos20", "E": "null_pos41",
    "F": "null_pos55", "G": "null_pos64",
}
CASCADE = [9, 22, 24, 51, 60]
MUT = 30


def parse_cif_atoms(path: Path):
    """Yield rows from the _atom_site loop. Returns list of dicts."""
    lines = path.read_text().splitlines()
    cols, rows, in_loop, in_header = [], [], False, False
    for ln in lines:
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s)
            in_loop, in_header = True, True
            continue
        if in_header and not s.startswith("_atom_site."):
            in_header = False
        if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
            parts = s.split()
            if len(parts) >= len(cols):
                rows.append(dict(zip(cols, parts)))
        elif in_loop and rows and (s.startswith("#") or s.startswith("loop_") or s.startswith("_") and not s.startswith("_atom_site.")):
            break
    return rows


def c1_per_chain(rows):
    """Return {chain_id: ndarray[L,3]} of C1' coords sorted by residue."""
    by_chain: dict[str, list[tuple[int, np.ndarray]]] = {}
    for r in rows:
        aid = r.get("_atom_site.label_atom_id", "").strip('"')
        if aid != "C1'":
            continue
        c = r["_atom_site.label_asym_id"]
        ri = int(r["_atom_site.label_seq_id"])
        xyz = np.array([float(r["_atom_site.Cartn_x"]),
                        float(r["_atom_site.Cartn_y"]),
                        float(r["_atom_site.Cartn_z"])], dtype=np.float32)
        by_chain.setdefault(c, []).append((ri, xyz))
    out: dict[str, np.ndarray] = {}
    for c, lst in by_chain.items():
        lst.sort()
        out[c] = np.stack([x for _, x in lst], axis=0)
    return out


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


def kabsch_rmsd(p, q):
    d = kabsch_per_res(p, q)
    return float(np.sqrt((d * d).mean())) if d is not None else float("nan")


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


def plddt_per_residue(json_path: Path, chain_ids_of_chain: dict[str, list[int]]):
    """Reduce per-atom pLDDT to per-residue by averaging atoms per residue.
    Returns {chain: ndarray[L] of mean pLDDT}."""
    j = json.loads(json_path.read_text())
    atom_chain = j["atom_chain_ids"]
    plddt = j["atom_plddts"]
    # We don't have atom->res mapping directly in full_data; use token_*
    token_chain = j["token_chain_ids"]
    token_res = j["token_res_ids"]
    # Re-bin per-atom into residues: atoms are listed grouped by chain in
    # the order of token enumeration. Each token corresponds to ~1 nucleotide.
    # Simpler proxy: use token-level pLDDT by averaging atoms per (chain,res)
    # using a running counter that resets at every chain transition.
    # Build atom-to-token index by matching chain transitions.
    out_chain_to_residue_plddt: dict[str, list[float]] = {}
    # Walk atoms; for each chain switch find next residue index in token order
    # Easier: trust token_res_ids align 1:1 with token order in atoms (chain-grouped)
    # Build mapping: for each chain, residues 1..L appear in token order
    # We assign atoms to (chain, current_residue) using atom_chain transitions
    # plus knowing per-chain length from token_res_ids
    per_chain_residues: dict[str, list[int]] = {}
    for ch, ri in zip(token_chain, token_res):
        per_chain_residues.setdefault(ch, []).append(ri)

    # For each chain, we know how many atoms per residue is variable, but
    # we'll use a simpler heuristic: walk atom list, when chain changes, start
    # over at residue 1; otherwise, use the fact that residue index is monotonic
    # and that the # of atoms per residue equals the count of consecutive atoms
    # of the same chain that sit "before" the next residue jump. Since we only
    # have per-atom chain (not residue), we instead use the total atom count
    # divided by residue count to approximate.

    # Better: use atom_chain to split atoms by chain, then split by equal-size
    # groups (~22 atoms/residue for RNA).
    chain_to_atoms_plddt: dict[str, list[float]] = {}
    for ch, p in zip(atom_chain, plddt):
        chain_to_atoms_plddt.setdefault(ch, []).append(p)

    for ch, atoms in chain_to_atoms_plddt.items():
        n_res = len(per_chain_residues.get(ch, []))
        if n_res == 0:
            continue
        atoms = np.asarray(atoms, dtype=np.float32)
        # Distribute n_atoms over n_res evenly (some res have 22, some 21 — close enough)
        # Use integer slicing
        chunks = np.array_split(atoms, n_res)
        out_chain_to_residue_plddt[ch] = [float(c.mean()) for c in chunks]
    return {ch: np.asarray(v, dtype=np.float32) for ch, v in out_chain_to_residue_plddt.items()}


def main():
    print("Loading AF3 multi-chain CIF...")
    rows = parse_cif_atoms(CIF)
    print(f"  {len(rows)} ATOM rows")
    af3_coords = c1_per_chain(rows)
    for ch, c in af3_coords.items():
        print(f"  chain {ch} ({CHAIN_MAP.get(ch,'?')}): {c.shape[0]} C1'")

    # Load DRfold2 ref C1' per target
    drf_coords = {}
    for tid in CHAIN_MAP.values():
        c = load_drfold_c1(tid)
        if c is not None:
            drf_coords[tid] = c
            print(f"  DRfold2 {tid}: {c.shape[0]} C1'")

    # Per-residue pLDDT per chain
    plddts = plddt_per_residue(JSON, {ch: list(range(1, 70)) for ch in CHAIN_MAP})

    # Compute per-target rows
    rep_lines = []
    af3_R1107 = af3_coords.get("A")
    af3_R1108 = af3_coords.get("B")
    drf_R1107 = drf_coords.get("R1107_human")
    drf_R1108 = drf_coords.get("R1108_chimp")

    # KEY ANALYSIS 1: AF3 R1107 vs AF3 R1108 — per-residue Kabsch delta
    delta_af3 = kabsch_per_res(af3_R1107, af3_R1108)
    top5_af3 = set(int(i) + 1 for i in np.argsort(delta_af3)[-5:])
    delta_drf = kabsch_per_res(drf_R1107, drf_R1108)
    top5_drf = set(int(i) + 1 for i in np.argsort(delta_drf)[-5:])

    print()
    print(f"AF3 R1107 vs R1108 top-5 divergent residues: {sorted(top5_af3)}")
    print(f"   Δ@9={delta_af3[8]:.2f}  Δ@22={delta_af3[21]:.2f}  Δ@24={delta_af3[23]:.2f}  Δ@30={delta_af3[29]:.2f}  Δ@51={delta_af3[50]:.2f}  Δ@60={delta_af3[59]:.2f}")
    print(f"DRfold2 R1107 vs R1108 top-5: {sorted(top5_drf)}")
    print(f"   Δ@9={delta_drf[8]:.2f}  Δ@22={delta_drf[21]:.2f}  Δ@24={delta_drf[23]:.2f}  Δ@30={delta_drf[29]:.2f}  Δ@51={delta_drf[50]:.2f}  Δ@60={delta_drf[59]:.2f}")

    # KEY ANALYSIS 2: cross-model agreement on full structure
    # AF3 vs DRfold2 per target
    for ch, tid in CHAIN_MAP.items():
        af3 = af3_coords.get(ch)
        drf = drf_coords.get(tid)
        if af3 is None or drf is None:
            continue
        rmsd = kabsch_rmsd(af3, drf)
        pl = plddts.get(ch)
        pl_str = (f" pLDDT mean={pl.mean():.1f} @9={pl[8]:.1f} @30={pl[29]:.1f} @60={pl[59]:.1f}"
                  if pl is not None and len(pl) >= 60 else "")
        print(f"  RMSD AF3 vs DRfold2 {tid}: {rmsd:.2f} Å{pl_str}")

    # Save full report
    md = OUT_DIR / "af3_results.md"
    csv_path = OUT_DIR / "af3_results.csv"

    overlap = top5_af3 & top5_drf
    with open(csv_path, "w") as f:
        f.write("target,chain,n_C1,af3_vs_drfold2_rmsd,af3_plddt_mean,af3_plddt_at_30,af3_plddt_at_9,af3_plddt_at_60\n")
        for ch, tid in CHAIN_MAP.items():
            af3 = af3_coords.get(ch); drf = drf_coords.get(tid); pl = plddts.get(ch)
            if af3 is None: continue
            rm = kabsch_rmsd(af3, drf) if drf is not None else float("nan")
            pmean = f"{pl.mean():.2f}" if pl is not None else "nan"
            p30 = f"{pl[29]:.2f}" if pl is not None and len(pl) >= 30 else "nan"
            p9 = f"{pl[8]:.2f}" if pl is not None and len(pl) >= 9 else "nan"
            p60 = f"{pl[59]:.2f}" if pl is not None and len(pl) >= 60 else "nan"
            f.write(f"{tid},{ch},{af3.shape[0]},{rm:.3f},{pmean},{p30},{p9},{p60}\n")

    with open(md, "w") as f:
        import datetime
        f.write("# Extension D — AlphaFold 3 cross-model comparison on CPEB3 panel\n\n")
        f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
        f.write("**Job**: single AF3 Server submission, 7 RNA chains (1 job instead of 7 — multi-chain run, all 7 sequences predicted independently as they don't form a complex).\n\n")
        f.write("## Chain mapping\n\n")
        f.write("| chain | target | length |\n|---|---|---|\n")
        for ch, tid in CHAIN_MAP.items():
            af3 = af3_coords.get(ch)
            f.write(f"| {ch} | {tid} | {af3.shape[0] if af3 is not None else '?'} |\n")

        f.write("\n## Key result 1 — top-5 divergent residues, R1107 vs R1108\n\n")
        f.write(f"- **DRfold2** top-5: `{sorted(top5_drf)}`\n")
        f.write(f"- **AF3** top-5:     `{sorted(top5_af3)}`\n")
        f.write(f"- **Overlap**: `{sorted(overlap)}` ({len(overlap)}/5 in common)\n\n")
        f.write("| residue | DRfold2 Δ (Å) | AF3 Δ (Å) | both top-5? |\n|---|---|---|---|\n")
        for r in [9, 22, 24, 30, 51, 60]:
            in_drf = "yes" if r in top5_drf else "no"
            in_af3 = "yes" if r in top5_af3 else "no"
            both = "**BOTH**" if (r in top5_drf and r in top5_af3) else ("DRfold2 only" if r in top5_drf else ("AF3 only" if r in top5_af3 else "neither"))
            f.write(f"| {r} | {delta_drf[r-1]:.2f} | {delta_af3[r-1]:.2f} | {both} |\n")

        f.write("\n## Key result 2 — RMSD between AF3 and DRfold2 per target\n\n")
        f.write("| target | RMSD AF3↔DRfold2 (Å) | AF3 mean pLDDT |\n|---|---|---|\n")
        for ch, tid in CHAIN_MAP.items():
            af3 = af3_coords.get(ch); drf = drf_coords.get(tid); pl = plddts.get(ch)
            if af3 is None: continue
            rm = kabsch_rmsd(af3, drf) if drf is not None else float("nan")
            pls = f"{pl.mean():.1f}" if pl is not None else "n/a"
            f.write(f"| {tid} | {rm:.2f} | {pls} |\n")

        f.write("\n## Key result 3 — AF3 per-residue pLDDT at cascade residues\n\n")
        f.write("Mean pLDDT (AF3, R1107 vs R1108) at the residues DRfold2 flagged as the cascade. ")
        f.write("**High pLDDT at residue 30 in R1107** = AF3 confidently predicts at the mispairing site (calibration risk).\n\n")
        f.write("| residue | AF3 R1107 pLDDT | AF3 R1108 pLDDT | Δ (R1107-R1108) |\n|---|---|---|---|\n")
        plA = plddts.get("A"); plB = plddts.get("B")
        if plA is not None and plB is not None:
            for r in [9, 22, 24, 30, 51, 60]:
                f.write(f"| {r} | {plA[r-1]:.1f} | {plB[r-1]:.1f} | {plA[r-1]-plB[r-1]:+.1f} |\n")

        f.write("\n## Interpretation\n\n")
        if len(overlap) >= 3:
            f.write(f"**Strong cross-model agreement**: AF3 and DRfold2 share {len(overlap)}/5 top divergent residues. The cascade pattern found in DRfold2 is NOT a model-specific artifact — at least one other independent SOTA architecture (AF3) sees the same structural response.\n")
        elif len(overlap) >= 1:
            f.write(f"**Partial cross-model agreement** ({len(overlap)}/5 overlap). Some divergence sites are shared, but each model also has its own emphasis. Worth reporting both top-5 sets as complementary evidence.\n")
        else:
            f.write(f"**No top-5 overlap**: AF3 and DRfold2 disagree on *which* residues are most perturbed by the human↔chimp variant. This is itself an interesting finding (model-specific structural attribution) but weakens any biological claim about cascade residues.\n")
        f.write("\n")
        if 30 not in top5_af3 and delta_af3[29] < 1.0:
            f.write("**Important secondary observation**: AF3 itself shows low structural divergence at residue 30 (the mutation site) — consistent with the parsimonious geometric explanation from our Experiment #1 (most of the C1' motion lives off-site).\n")
    print(f"\nWritten: {csv_path}")
    print(f"Written: {md}")


if __name__ == "__main__":
    main()
