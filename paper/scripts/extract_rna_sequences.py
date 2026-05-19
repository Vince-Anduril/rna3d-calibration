"""Extract RNA sequence from each multi-RNA target CIF.

For each target, find the RNA chain(s), extract per-residue type (A/U/G/C),
and write input.fasta for DRfold2 + a manifest for AF3 submission.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path("/Users/leduigouvincent/rna3d-mirror")
DATA = ROOT / "data/multi_rna"

# CASE-MAPPED 3-letter to 1-letter for RNA
RNA_MAP = {"A": "A", "U": "U", "G": "G", "C": "C",
           "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
           "RA": "A", "RU": "U", "RG": "G", "RC": "C"}


def parse_cif_chains(path: Path):
    """For each chain, list (residue_index, comp_id) of RNA residues only."""
    lines = path.read_text().splitlines()
    cols, data_start = [], None
    in_header = False
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s); in_header = True
        elif in_header and not s.startswith("_atom_site."):
            data_start = i; break
    chains: dict[str, dict[int, str]] = {}
    for ln in lines[data_start:]:
        if not ln.startswith("ATOM"): continue
        parts = ln.split()
        if len(parts) < len(cols): continue
        d = dict(zip(cols, parts))
        aid = d["_atom_site.label_atom_id"].strip('"')
        if aid != "C1'": continue
        comp = d["_atom_site.label_comp_id"]
        if comp not in RNA_MAP: continue
        ch = d["_atom_site.label_asym_id"]
        try:
            ri = int(d["_atom_site.label_seq_id"])
        except (KeyError, ValueError):
            continue
        chains.setdefault(ch, {})[ri] = RNA_MAP[comp]
    return chains


targets = ["9LJN", "9UW0", "9HRD", "12CI"]

manifest = []
for tid in targets:
    cif = DATA / f"{tid}.cif"
    chains = parse_cif_chains(cif)
    print(f"\n=== {tid} ===")
    for ch, ridx in chains.items():
        residues = sorted(ridx.items())
        seq = "".join(r for _, r in residues)
        nres = len(seq)
        print(f"  chain {ch}: {nres} nt | {seq}")
    # Pick best chain: first chain by length within 50-150
    candidates = [(ch, ridx) for ch, ridx in chains.items() if 50 <= len(ridx) <= 150]
    if not candidates:
        print(f"  WARNING: no suitable chain in {tid}")
        continue
    # Sort by length, then alphabetically
    candidates.sort(key=lambda x: (-len(x[1]), x[0]))
    pick_ch, pick_ridx = candidates[0]
    seq = "".join(r for _, r in sorted(pick_ridx.items()))
    print(f"  → pick chain {pick_ch}, {len(seq)} nt: {seq}")
    # Write FASTA
    fa = DATA / f"{tid}.fasta"
    fa.write_text(f">{tid}_{pick_ch}\n{seq}\n")
    manifest.append({"target": tid, "chain": pick_ch, "length": len(seq), "sequence": seq})

(DATA / "manifest.json").write_text(json.dumps(manifest, indent=2))
print(f"\nManifest: {DATA / 'manifest.json'}")
for m in manifest:
    print(f"  {m['target']}_{m['chain']} ({m['length']} nt) → {DATA / (m['target']+'.fasta')}")
