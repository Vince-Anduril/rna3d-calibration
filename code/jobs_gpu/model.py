"""RNA structure prediction model.

A small transformer encoder + per-residue 3D coordinate head + confidence head.
Designed to be trained from scratch on the Stanford RNA 3D Folding dataset to
serve as a controllable predictor for calibration analysis on the ribozyme focus
class.

Architecture (default config, ~7M params):
    - Embedding: vocab=8, d_model=256
    - Sinusoidal positional encoding
    - 6 transformer encoder layers, 8 heads, FF=1024
    - Coord head: 3-layer MLP -> (x, y, z) per residue
    - Confidence head: 3-layer MLP -> sigmoid scalar per residue
"""
from __future__ import annotations
import math
from dataclasses import dataclass
import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class ModelConfig:
    vocab_size: int = 8       # [PAD] [MASK] [CLS] A U G C N
    d_model: int = 256
    n_heads: int = 8
    n_layers: int = 6
    d_ff: int = 1024
    max_len: int = 512
    dropout: float = 0.1
    pad_id: int = 0


class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        pos = torch.arange(max_len).unsqueeze(1).float()
        div = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(pos * div)
        pe[:, 1::2] = torch.cos(pos * div)
        self.register_buffer("pe", pe.unsqueeze(0))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.pe[:, : x.size(1)]


class RNAStructureModel(nn.Module):
    """Encoder + dual head (3D coords + confidence).

    Inputs:
        token_ids: LongTensor [B, L]
    Returns dict with:
        mlm_logits: [B, L, vocab] — masked-LM logits (for stage-1 pretraining)
        coords:     [B, L, 3]     — predicted xyz per residue
        confidence: [B, L]        — sigmoid confidence per residue
    """

    def __init__(self, cfg: ModelConfig):
        super().__init__()
        self.cfg = cfg
        self.embed = nn.Embedding(cfg.vocab_size, cfg.d_model, padding_idx=cfg.pad_id)
        self.pos = PositionalEncoding(cfg.d_model, cfg.max_len)
        layer = nn.TransformerEncoderLayer(
            d_model=cfg.d_model,
            nhead=cfg.n_heads,
            dim_feedforward=cfg.d_ff,
            dropout=cfg.dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=cfg.n_layers)

        self.mlm_head = nn.Linear(cfg.d_model, cfg.vocab_size)

        self.coord_head = nn.Sequential(
            nn.Linear(cfg.d_model, cfg.d_model),
            nn.GELU(),
            nn.Linear(cfg.d_model, cfg.d_model),
            nn.GELU(),
            nn.Linear(cfg.d_model, 3),
        )

        self.confidence_head = nn.Sequential(
            nn.Linear(cfg.d_model, cfg.d_model // 2),
            nn.GELU(),
            nn.Linear(cfg.d_model // 2, cfg.d_model // 4),
            nn.GELU(),
            nn.Linear(cfg.d_model // 4, 1),
        )

    def forward(self, token_ids: torch.Tensor) -> dict:
        pad_mask = token_ids == self.cfg.pad_id
        h = self.embed(token_ids)
        h = self.pos(h)
        h = self.encoder(h, src_key_padding_mask=pad_mask)
        return {
            "mlm_logits": self.mlm_head(h),
            "coords": self.coord_head(h),
            "confidence": torch.sigmoid(self.confidence_head(h).squeeze(-1)),
            "pad_mask": pad_mask,
        }

    @property
    def num_params(self) -> int:
        return sum(p.numel() for p in self.parameters())
