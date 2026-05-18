"""Defensive tests for structure_loss and Kabsch alignment.

These tests catch the kind of bug that would silently waste a 3h pod run:
- Kabsch alignment producing wrong rotation
- Loss being non-zero on a perfectly-rotated copy of the target
- Mask being applied incorrectly
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "jobs_gpu"))

import numpy as np
import pytest
import torch
from loss import _kabsch_align, structure_loss, confidence_loss, mlm_loss


def test_kabsch_aligns_rotated_copy_to_zero_error():
    """Apply a known 90° rotation to a tetrahedron. Kabsch should recover it
    and the aligned coords should match the original to machine precision."""
    torch.manual_seed(0)
    # 4 residues, batch=1
    P_orig = torch.tensor([[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]])
    # Rotation around z-axis 90° + translation
    R = torch.tensor([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    t = torch.tensor([5.0, -3.0, 2.0])
    P_rotated = (P_orig.squeeze(0) @ R.T + t).unsqueeze(0)

    mask = torch.ones(1, 4, dtype=torch.bool)
    aligned = _kabsch_align(P_rotated, P_orig, mask)
    err = (aligned - P_orig).abs().max().item()
    assert err < 1e-5, f"Kabsch failed to recover known rotation; max abs err = {err}"


def test_structure_loss_zero_on_identical_inputs():
    """If pred == true, coord_mse must be 0 (pair_mse may be small due to random sampling)."""
    torch.manual_seed(0)
    pred = torch.randn(2, 16, 3)
    true = pred.clone()
    mask = torch.ones(2, 16, dtype=torch.bool)
    out = structure_loss(pred, true, mask)
    assert out["coord_mse"].item() < 1e-6
    assert out["pair_mse"].item() < 1e-6


def test_structure_loss_zero_on_rotated_copy():
    """A rotated+translated copy should give ~0 coord_mse after Kabsch."""
    torch.manual_seed(1)
    true = torch.randn(1, 16, 3)
    # Apply a deterministic rotation
    R = torch.tensor([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    pred = (true.squeeze(0) @ R.T + torch.tensor([3., 2., 1.])).unsqueeze(0)
    mask = torch.ones(1, 16, dtype=torch.bool)
    out = structure_loss(pred, true, mask)
    assert out["coord_mse"].item() < 1e-4, f"Kabsch-aligned coord_mse not ~0: {out['coord_mse'].item()}"


def test_structure_loss_positive_on_perturbed_input():
    """A perturbed copy should give positive coord_mse roughly equal to the perturbation variance."""
    torch.manual_seed(2)
    true = torch.randn(1, 16, 3)
    pred = true + 0.5 * torch.randn(1, 16, 3)
    mask = torch.ones(1, 16, dtype=torch.bool)
    out = structure_loss(pred, true, mask)
    # The perturbation variance is 0.25 per axis, total squared ~0.75 per residue
    # After Kabsch alignment we expect the global rotation/translation degrees of
    # freedom to be removed, so the residual should be a bit less.
    assert 0.3 < out["coord_mse"].item() < 1.5


def test_structure_loss_respects_mask():
    """Residues with valid=False should not contribute to coord_mse."""
    torch.manual_seed(3)
    true = torch.randn(1, 16, 3)
    pred = true.clone()
    # Corrupt residues 10-15 — but mask them out
    pred[0, 10:] = 999.0
    mask = torch.zeros(1, 16, dtype=torch.bool)
    mask[0, :10] = True
    out = structure_loss(pred, true, mask)
    assert out["coord_mse"].item() < 1e-4, "Masked-out residues should not contribute"


def test_confidence_loss_zero_on_perfect_predictor():
    """If predictions are perfect, target = exp(0/5) = 1; confidence sigmoid≈1 → BCE ≈ 0."""
    torch.manual_seed(4)
    true = torch.randn(2, 8, 3)
    pred = true.clone()
    # Force confidence very close to 1
    confidence = torch.full((2, 8), 0.999)
    mask = torch.ones(2, 8, dtype=torch.bool)
    loss = confidence_loss(confidence, pred, true, mask)
    assert loss.item() < 0.05, f"Confidence loss too high on perfect predictor: {loss.item()}"


def test_confidence_loss_higher_when_overconfident_on_bad_predictions():
    """Bad predictions + high confidence → high loss; bad predictions + low confidence → lower loss."""
    torch.manual_seed(5)
    true = torch.randn(1, 8, 3)
    bad_pred = true + 20.0 * torch.randn(1, 8, 3)  # very far off
    mask = torch.ones(1, 8, dtype=torch.bool)
    conf_high = torch.full((1, 8), 0.99)
    conf_low = torch.full((1, 8), 0.01)
    loss_high = confidence_loss(conf_high, bad_pred, true, mask)
    loss_low = confidence_loss(conf_low, bad_pred, true, mask)
    assert loss_high.item() > loss_low.item(), \
        f"Overconfidence on bad predictions should hurt more: high={loss_high.item()}, low={loss_low.item()}"


def test_mlm_loss_ignores_minus_100():
    """Labels of -100 are ignored; remaining positions form the loss."""
    torch.manual_seed(6)
    logits = torch.randn(2, 4, 8)
    labels = torch.tensor([[3, -100, 4, -100], [-100, 5, -100, 6]])
    loss = mlm_loss(logits, labels)
    assert loss.isfinite().item()
    # Should be approximately log(8) ≈ 2.08 for random logits on 4 valid positions
    assert 1.0 < loss.item() < 3.5
