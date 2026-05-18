"""Data loading: read Stanford CSVs + parse PDB ground-truth into per-residue coords.

Designed to run on the pod. For Mac smoke tests, use `SyntheticRNADataset`.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import random
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

VOCAB = {"[PAD]": 0, "[MASK]": 1, "[CLS]": 2, "A": 3, "U": 4, "G": 5, "C": 6, "N": 7}
PAD_ID = VOCAB["[PAD]"]
MASK_ID = VOCAB["[MASK]"]
CLS_ID = VOCAB["[CLS]"]


def encode(seq: str, max_len: int) -> list[int]:
    ids = [CLS_ID] + [VOCAB.get(c.upper(), VOCAB["N"]) for c in seq][: max_len - 1]
    return ids[:max_len]


@dataclass
class StructureSample:
    target_id: str
    split: str
    sequence: str
    coords: np.ndarray | None  # [L, 3] or None if no ground-truth
    valid: np.ndarray | None   # [L] bool


def load_pdb_coords(pdb_path: Path) -> np.ndarray | None:
    """Extract per-residue C1' or P coordinates. Returns [L, 3] or None on failure.

    Uses biotite if available; falls back to a simple parser.
    """
    try:
        import biotite.structure.io as bsio
        s = bsio.load_structure(str(pdb_path))
        # Prefer C1' for RNA backbone; fall back to P
        for name in ("C1'", "P"):
            mask = s.atom_name == name
            if mask.any():
                coords = s.coord[mask]
                return np.asarray(coords, dtype=np.float32)
    except Exception:
        pass
    # naive fallback: parse PDB ATOM lines
    coords = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() in ("C1'", "P"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        if not coords:
            return None
        return np.asarray(coords, dtype=np.float32)
    except Exception:
        return None


class StanfordRNADataset(Dataset):
    """Sequence + PDB ground-truth dataset, with category labels for stratified eval."""

    def __init__(
        self,
        data_root: Path,
        splits: list[str] = ("train", "val", "test"),
        max_len: int = 256,
        category_csv: Path | None = None,  # from run 01b
    ):
        self.data_root = Path(data_root)
        self.max_len = max_len
        files = {
            "train": "train_sequences.csv",
            "val": "validation_sequences.csv",
            "test": "test_sequences.csv",
        }
        dfs = []
        for sp in splits:
            df = pd.read_csv(self.data_root / "kaggle_raw" / files[sp])
            df["split"] = sp
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True)
        df["length"] = df["sequence"].str.len()

        # Optional category labels (from run 01b)
        if category_csv is not None and Path(category_csv).exists():
            cats = pd.read_csv(category_csv)[["target_id", "category"]]
            df = df.merge(cats, on="target_id", how="left")
            df["category"] = df["category"].fillna("other")
        else:
            df["category"] = "unknown"

        self.df = df
        self.pdb_dir = self.data_root / "kaggle_raw" / "PDB_RNA"

    def __len__(self) -> int:
        return len(self.df)

    def __getitem__(self, idx: int) -> dict:
        row = self.df.iloc[idx]
        seq = row["sequence"]
        ids = encode(seq, self.max_len)
        ids = ids + [PAD_ID] * (self.max_len - len(ids))
        token_ids = torch.tensor(ids, dtype=torch.long)

        # Try to load coords for this target_id; if none, return zeros with valid=0
        coords = np.zeros((self.max_len, 3), dtype=np.float32)
        valid = np.zeros(self.max_len, dtype=bool)
        pdb_candidates = [
            self.pdb_dir / f"{row['target_id']}.pdb",
            self.pdb_dir / f"{row['target_id']}.cif",
        ]
        for p in pdb_candidates:
            if p.exists():
                cc = load_pdb_coords(p)
                if cc is not None and len(cc) > 0:
                    n = min(len(cc), self.max_len - 1)  # leave room for CLS
                    coords[1 : 1 + n] = cc[:n]
                    valid[1 : 1 + n] = True
                break

        return {
            "token_ids": token_ids,
            "coords": torch.from_numpy(coords),
            "valid": torch.from_numpy(valid),
            "target_id": row["target_id"],
            "split": row["split"],
            "category": row["category"],
            "length": int(row["length"]),
        }


def make_mlm_batch(
    token_ids: torch.Tensor, mlm_prob: float = 0.15
) -> tuple[torch.Tensor, torch.Tensor]:
    """Returns (masked_input, labels) where labels are -100 for unmasked positions."""
    labels = token_ids.clone()
    mask = (torch.rand_like(token_ids, dtype=torch.float) < mlm_prob) & (token_ids != PAD_ID) & (token_ids != CLS_ID)
    labels[~mask] = -100
    masked = token_ids.clone()
    masked[mask] = MASK_ID
    return masked, labels


class SyntheticRNADataset(Dataset):
    """Tiny dataset for Mac CPU smoke tests. 32 fake sequences with fake coords."""

    def __init__(self, n: int = 32, max_len: int = 64):
        self.max_len = max_len
        rng = np.random.RandomState(42)
        self.samples = []
        for i in range(n):
            L = rng.randint(20, max_len - 4)
            seq = "".join(rng.choice(list("AUGC"), size=L))
            ids = encode(seq, max_len)
            ids = ids + [PAD_ID] * (max_len - len(ids))
            coords = np.zeros((max_len, 3), dtype=np.float32)
            valid = np.zeros(max_len, dtype=bool)
            coords[1 : 1 + L] = rng.randn(L, 3).astype(np.float32) * 5
            valid[1 : 1 + L] = True
            self.samples.append({
                "token_ids": torch.tensor(ids, dtype=torch.long),
                "coords": torch.from_numpy(coords),
                "valid": torch.from_numpy(valid),
                "target_id": f"FAKE_{i:03d}",
                "split": "train",
                "category": "synthetic",
                "length": L,
            })

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, i):
        return self.samples[i]
