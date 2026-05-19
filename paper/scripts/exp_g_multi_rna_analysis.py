"""Experiment G — cross-family DRfold2 vs AF3 vs crystal analysis.

Targets: 5 RNAs (CPEB3 / 7QR3 + 4 new from data/multi_rna/).
For each target, compute:
  - DRfold2 vs crystal (C1' Kabsch RMSD)
  - AF3 single-chain (5 seeds), default & best
  - DRfold2 ensemble std (per-residue) → calibration vs error
  - AF3 per-residue pLDDT → calibration vs error

Output: multi_rna_results.{csv,md} + per-target detail
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")

# --- Target manifest ---
TARGETS = [
    {"id": "7QR3",  "crystal": ROOT / "data/crystal/7QR3.cif", "crystal_chain": "C",
     "af3_dir": ROOT / "data/af3_multi/raw/fold_r1108_chimp_single",
     "drf_dir": ROOT / ".research/30_experiments/runs/cpeb3_focused/R1108_chimp/drfold2",
     "name": "CPEB3 ribozyme (R1108)", "family": "HDV ribozyme", "deposit": "2021-12", "length": 69,
     "drfold_seen": "likely yes (Lee 2025 trained post-2022)",
     "af3_seen": "unlikely (cutoff Sep 2021)"},
    {"id": "9LJN",  "crystal": ROOT / "data/multi_rna/9LJN.cif", "crystal_chain": "A",
     "af3_dir": ROOT / "data/af3_multi/raw/fold_9ljn_guanineii",
     "drf_dir": ROOT / ".research/30_experiments/runs/multi_rna/9LJN/drfold2",
     "name": "Guanine-II riboswitch", "family": "riboswitch", "deposit": "2025-01", "length": 71,
     "drfold_seen": "unlikely", "af3_seen": "no"},
    {"id": "9UW0",  "crystal": ROOT / "data/multi_rna/9UW0.cif", "crystal_chain": "A",
     "af3_dir": ROOT / "data/af3_multi/raw/fold_9uw0_2dg",
     "drf_dir": ROOT / ".research/30_experiments/runs/multi_rna/9UW0/drfold2",
     "name": "2'-dG-III riboswitch", "family": "riboswitch", "deposit": "2025-05", "length": 63,
     "drfold_seen": "unlikely", "af3_seen": "no"},
    {"id": "9HRD",  "crystal": ROOT / "data/multi_rna/9HRD.cif", "crystal_chain": "A",
     "af3_dir": ROOT / "data/af3_multi/raw/fold_9hrd_gtpapt",
     "drf_dir": ROOT / ".research/30_experiments/runs/multi_rna/9HRD/drfold2",
     "name": "Class V GTP aptamer", "family": "aptamer", "deposit": "2024-12", "length": 67,
     "drfold_seen": "unlikely", "af3_seen": "no"},
    {"id": "12CI",  "crystal": ROOT / "data/multi_rna/12CI.cif", "crystal_chain": "A",
     "af3_dir": ROOT / "data/af3_multi/raw/fold_12ci_dopamine",
     "drf_dir": ROOT / ".research/30_experiments/runs/multi_rna/12CI/drfold2",
     "name": "Dopamine aptamer DGR-1A", "family": "aptamer", "deposit": "2026-03", "length": 82,
     "drfold_seen": "no", "af3_seen": "no"},
]

CASCADE = [9, 22, 24, 51, 60]


def parse_cif_c1(path: Path, chain_filter: str | None = None, standard_only: bool = True):
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
    STANDARD = {"A", "U", "G", "C"}
    for ln in lines[data_start:]:
        if not (ln.startswith("ATOM") or (not standard_only and ln.startswith("HETATM"))):
            continue
        parts = ln.split()
        if len(parts) < len(cols): continue
        d = dict(zip(cols, parts))
        aid = d["_atom_site.label_atom_id"].strip('"')
        if aid != "C1'": continue
        comp = d["_atom_site.label_comp_id"]
        if standard_only and comp not in STANDARD: continue
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


def pdb_c1(path: Path):
    coords = []
    for ln in path.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "C1'":
            try:
                coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError: continue
    return np.array(coords, dtype=np.float32) if coords else None


def kabsch_per_res(p, q):
    if p is None or q is None or p.shape != q.shape: return None
    pc, qc = p.mean(0), q.mean(0); pp, qq = p - pc, q - qc
    H = pp.T @ qq
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    al = (R @ pp.T).T + qc
    return np.linalg.norm(al - q, axis=-1)


def rmsd(p, q):
    pr = kabsch_per_res(p, q)
    return float(np.sqrt((pr*pr).mean())) if pr is not None else float("nan")


def af3_plddt_per_residue(json_path: Path, chain: str | None = None):
    j = json.loads(json_path.read_text())
    atom_chain = j["atom_chain_ids"]
    plddt = j["atom_plddts"]
    token_chain = j["token_chain_ids"]
    if chain is None:
        chain = atom_chain[0]
    n_res = sum(1 for c in token_chain if c == chain)
    chain_atoms = [p for c, p in zip(atom_chain, plddt) if c == chain]
    chunks = np.array_split(np.asarray(chain_atoms, dtype=np.float32), n_res)
    return np.array([c.mean() for c in chunks], dtype=np.float32)


def analyze_target(t):
    print(f"\n{'='*60}\n{t['id']} — {t['name']}\n{'='*60}")
    out = {"id": t["id"], "name": t["name"], "family": t["family"],
           "length": t["length"], "deposit": t["deposit"],
           "drfold_seen": t["drfold_seen"], "af3_seen": t["af3_seen"]}

    # Crystal
    chains = parse_cif_c1(t["crystal"], chain_filter=t["crystal_chain"])
    cry = chains.get(t["crystal_chain"])
    if cry is None:
        print(f"  [SKIP] no crystal chain {t['crystal_chain']}")
        return out
    out["crystal_n_C1"] = int(cry.shape[0])
    print(f"  Crystal chain {t['crystal_chain']}: {cry.shape[0]} C1' residues")

    # DRfold2 selected
    drf_pdbs = sorted((t["drf_dir"] / "folds").glob("opt_0_*.pdb"))
    if not drf_pdbs:
        print(f"  [WAIT] no DRfold2 fold yet at {t['drf_dir']}")
        out["drfold2_rmsd_vs_crystal"] = None
        return out
    drf_c1 = pdb_c1(drf_pdbs[0])
    if drf_c1.shape[0] != cry.shape[0]:
        print(f"  [WARN] length mismatch: DRfold2 {drf_c1.shape[0]} vs crystal {cry.shape[0]}")
        out["drfold2_rmsd_vs_crystal"] = None
    else:
        out["drfold2_rmsd_vs_crystal"] = rmsd(drf_c1, cry)
        print(f"  DRfold2 vs crystal: {out['drfold2_rmsd_vs_crystal']:.2f} Å")

    # AF3 single-chain — 5 seeds
    af3_cifs = sorted(t["af3_dir"].glob("*model_*.cif"))
    af3_rmsds = []
    af3_plddts = []
    for cif in af3_cifs:
        af3_chains = parse_cif_c1(cif)
        if not af3_chains: continue
        # take the single chain (or chain A)
        pred = af3_chains["A"] if "A" in af3_chains else next(iter(af3_chains.values()))
        if pred.shape[0] != cry.shape[0]:
            af3_rmsds.append(float("nan")); af3_plddts.append(float("nan")); continue
        af3_rmsds.append(rmsd(pred, cry))
        # find matching summary json
        seed_id = cif.stem.split("_model_")[1]
        full_data = cif.parent / cif.name.replace(f"model_{seed_id}.cif", f"full_data_{seed_id}.json")
        if full_data.exists():
            pl = af3_plddt_per_residue(full_data, chain=list(af3_chains.keys())[0])
            af3_plddts.append(float(pl.mean()))
    out["af3_rmsds_per_seed"] = af3_rmsds
    out["af3_mean_plddts_per_seed"] = af3_plddts
    if af3_rmsds:
        out["af3_rmsd_default"] = af3_rmsds[0]
        out["af3_rmsd_best"] = min(r for r in af3_rmsds if not np.isnan(r))
        print(f"  AF3 single-chain: default seed={out['af3_rmsd_default']:.2f} Å, best of {len(af3_rmsds)}={out['af3_rmsd_best']:.2f} Å")
        print(f"  AF3 mean pLDDTs per seed: {[f'{p:.1f}' for p in af3_plddts]}")
        # Calibration: per-residue plddt vs per-residue err for default seed
        default_cif = af3_cifs[0]
        default_chains = parse_cif_c1(default_cif)
        default_pred = default_chains["A"] if "A" in default_chains else next(iter(default_chains.values()))
        full_data = default_cif.parent / default_cif.name.replace("model_0.cif", "full_data_0.json")
        if full_data.exists() and default_pred.shape[0] == cry.shape[0]:
            per_res_err = kabsch_per_res(default_pred, cry)
            per_res_plddt = af3_plddt_per_residue(full_data, chain=list(default_chains.keys())[0])
            try:
                from scipy.stats import spearmanr
                rho, p = spearmanr(per_res_plddt, per_res_err)
                out["af3_calib_spearman"] = float(rho)
                out["af3_calib_p"] = float(p)
                print(f"  AF3 pLDDT vs err: Spearman = {rho:+.3f}, p={p:.3g}")
            except Exception as e:
                print(f"  [calib failed] {e}")
    return out


if __name__ == "__main__":
    results = [analyze_target(t) for t in TARGETS]

    # Write CSV
    import csv
    csv_path = ROOT / ".research/30_experiments/runs/multi_rna/multi_rna_results.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["id", "name", "family", "length", "deposit", "drfold_seen", "af3_seen",
              "drfold2_rmsd_vs_crystal",
              "af3_rmsd_default", "af3_rmsd_best",
              "af3_calib_spearman", "af3_calib_p"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k) for k in fields})
    print(f"\nCSV: {csv_path}")

    # MD summary
    md = ROOT / ".research/30_experiments/runs/multi_rna/multi_rna_results.md"
    with open(md, "w") as f:
        import datetime
        f.write("# Cross-family analysis — Experiment G\n\n")
        f.write(f"Generated: {datetime.datetime.now(datetime.UTC).isoformat()}\n\n")
        f.write("All AF3 runs in this table are **single-chain** (one job per RNA).\n\n")
        f.write("## Ground-truth RMSD per target\n\n")
        f.write("| ID | RNA | family | L | deposit | DRfold2 seen? | AF3 seen? | DRfold2 (Å) | AF3 default (Å) | AF3 best (Å) |\n")
        f.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for r in results:
            drf = f"{r['drfold2_rmsd_vs_crystal']:.2f}" if r.get('drfold2_rmsd_vs_crystal') else "—"
            af3d = f"{r['af3_rmsd_default']:.2f}" if r.get('af3_rmsd_default') is not None else "—"
            af3b = f"{r['af3_rmsd_best']:.2f}" if r.get('af3_rmsd_best') is not None else "—"
            f.write(f"| {r['id']} | {r['name']} | {r['family']} | {r['length']} | {r['deposit']} | {r['drfold_seen']} | {r['af3_seen']} | {drf} | {af3d} | {af3b} |\n")

        f.write("\n## AF3 per-residue pLDDT calibration vs crystal\n\n")
        f.write("| ID | Spearman ρ | p |\n|---|---|---|\n")
        for r in results:
            rho = f"{r['af3_calib_spearman']:+.3f}" if r.get('af3_calib_spearman') is not None else "—"
            p = f"{r['af3_calib_p']:.3g}" if r.get('af3_calib_p') is not None else "—"
            f.write(f"| {r['id']} | {rho} | {p} |\n")

        # Aggregate stats
        drf_vals = [r['drfold2_rmsd_vs_crystal'] for r in results if r.get('drfold2_rmsd_vs_crystal') is not None]
        af3_def = [r['af3_rmsd_default'] for r in results if r.get('af3_rmsd_default') is not None]
        af3_best = [r['af3_rmsd_best'] for r in results if r.get('af3_rmsd_best') is not None]
        f.write(f"\n## Aggregate (N = {len(drf_vals)} targets with full data)\n\n")
        if drf_vals:
            f.write(f"- DRfold2 mean RMSD: **{np.mean(drf_vals):.2f} Å** (median {np.median(drf_vals):.2f})\n")
        if af3_def:
            f.write(f"- AF3 default-seed mean: **{np.mean(af3_def):.2f} Å** (median {np.median(af3_def):.2f})\n")
        if af3_best:
            f.write(f"- AF3 best-of-5 mean: **{np.mean(af3_best):.2f} Å** (median {np.median(af3_best):.2f})\n")
        if len(drf_vals) == len(af3_def) and len(drf_vals) >= 3:
            from scipy.stats import wilcoxon
            try:
                w, p = wilcoxon(drf_vals, af3_def)
                f.write(f"- Wilcoxon signed-rank (DRfold2 vs AF3 default): W={w:.1f}, p={p:.3g} (N={len(drf_vals)})\n")
            except Exception as e:
                f.write(f"- Wilcoxon failed: {e}\n")
    print(f"MD:  {md}")
