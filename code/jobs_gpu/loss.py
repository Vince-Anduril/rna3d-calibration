"""Losses for RNA structure prediction.

- mlm_loss: cross-entropy for masked-LM pretraining (stage 1)
- structure_loss: Kabsch-aligned MSE on coords + pair-distance auxiliary
- confidence_loss: encourages confidence ~ 1 - normalized_error

All losses respect a pad mask (B, L).
"""
from __future__ import annotations
import torch
import torch.nn as nn
import torch.nn.functional as F


def mlm_loss(logits: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    """labels with -100 are ignored."""
    return F.cross_entropy(
        logits.reshape(-1, logits.size(-1)),
        labels.reshape(-1),
        ignore_index=-100,
    )


def _kabsch_align(P: torch.Tensor, Q: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Align P to Q via Kabsch, per batch element.

    P, Q: [B, L, 3]; mask: [B, L] (True for valid residues)
    Returns P_aligned: [B, L, 3]
    """
    B = P.shape[0]
    orig_dtype = P.dtype
    out = torch.zeros_like(P)
    for b in range(B):
        m = mask[b]
        if m.sum() < 3:
            out[b] = P[b]
            continue
        # SVD on CUDA does not support bfloat16/float16 — promote to float32 around it.
        p = P[b, m].float()
        q = Q[b, m].float()
        pc, qc = p.mean(0), q.mean(0)
        pp, qq = p - pc, q - qc
        H = pp.T @ qq
        U, S, Vt = torch.linalg.svd(H)
        d = torch.sign(torch.det(Vt.T @ U.T))
        D = torch.diag(torch.tensor([1.0, 1.0, d], device=p.device, dtype=p.dtype))
        R = Vt.T @ D @ U.T
        p_full = P[b].float()
        p_full_centered = p_full - pc
        aligned = (R @ p_full_centered.T).T + qc
        out[b] = aligned.to(orig_dtype)
    return out


def structure_loss(
    pred_coords: torch.Tensor,
    true_coords: torch.Tensor,
    valid_mask: torch.Tensor,
    pair_weight: float = 0.5,
) -> dict:
    """Kabsch-aligned MSE + pair-distance loss.

    pred_coords, true_coords: [B, L, 3]
    valid_mask: [B, L] (True where residue has ground truth)
    Returns: dict with 'coord_mse', 'pair_mse', 'total'.
    """
    aligned = _kabsch_align(pred_coords, true_coords, valid_mask)
    diff_sq = ((aligned - true_coords) ** 2).sum(-1)  # [B, L]
    coord_mse = (diff_sq * valid_mask).sum() / valid_mask.sum().clamp(min=1)

    # Pair-distance auxiliary on a subset to save memory
    B, L, _ = pred_coords.shape
    pair_total = torch.tensor(0.0, device=pred_coords.device)
    pair_count = 0
    for b in range(B):
        m = valid_mask[b]
        n = int(m.sum().item())
        if n < 4:
            continue
        idx = torch.nonzero(m, as_tuple=True)[0]
        # Take up to 64 random pairs
        k = min(64, n * (n - 1) // 2)
        if k == 0:
            continue
        ii = idx[torch.randint(0, n, (k,), device=pred_coords.device)]
        jj = idx[torch.randint(0, n, (k,), device=pred_coords.device)]
        dp = torch.norm(pred_coords[b, ii] - pred_coords[b, jj], dim=-1)
        dt = torch.norm(true_coords[b, ii] - true_coords[b, jj], dim=-1)
        pair_total = pair_total + F.mse_loss(dp, dt, reduction="sum")
        pair_count = pair_count + k
    pair_mse = pair_total / max(pair_count, 1)

    return {
        "coord_mse": coord_mse,
        "pair_mse": pair_mse,
        "total": coord_mse + pair_weight * pair_mse,
    }


def confidence_loss(
    confidence: torch.Tensor,
    pred_coords: torch.Tensor,
    true_coords: torch.Tensor,
    valid_mask: torch.Tensor,
    error_scale: float = 5.0,
) -> torch.Tensor:
    """Encourage confidence ~= exp(-||pred - true|| / error_scale).

    confidence: [B, L] in [0, 1]
    """
    with torch.no_grad():
        # NOTE: this uses the model's pre-alignment predictions for the error
        # signal; we accept this as a simple proxy. A more principled
        # alternative would use Kabsch-aligned errors here too.
        err = torch.norm(pred_coords.detach() - true_coords, dim=-1)  # [B, L]
        target = torch.exp(-err / error_scale)
    bce = F.binary_cross_entropy(confidence, target, reduction="none")
    return (bce * valid_mask).sum() / valid_mask.sum().clamp(min=1)
