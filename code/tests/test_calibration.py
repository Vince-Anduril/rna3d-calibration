"""Unit tests for calibration metrics. Pure numpy, runs fast on Mac CPU."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "code"))

import numpy as np
import pytest
from calibration import (
    expected_calibration_error, maximum_calibration_error,
    brier_score, errors_to_correctness, bootstrap_ci, reliability_table,
)


def test_ece_perfect_calibration_is_zero():
    confs = np.concatenate([np.full(10, 0.3), np.full(10, 0.8)])
    correct = np.concatenate([
        np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
        np.array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0]),
    ])
    assert expected_calibration_error(confs, correct, n_bins=2) == pytest.approx(0.0, abs=1e-12)


def test_ece_systematic_overconfidence():
    confs = np.full(10, 0.9)
    correct = np.array([1.0, 1, 1, 1, 1, 0, 0, 0, 0, 0])  # 0.5 correct
    assert expected_calibration_error(confs, correct, n_bins=1) == pytest.approx(0.4)


def test_ece_weighted_average():
    confs = np.concatenate([np.full(8, 0.2), np.full(2, 0.9)])
    correct = np.concatenate([
        np.array([1, 1, 1, 1, 0, 0, 0, 0]),  # 0.5 in bin A
        np.array([1, 1]),                      # 1.0 in bin B
    ])
    # (8/10)*|0.2-0.5| + (2/10)*|0.9-1.0| = 0.24 + 0.02 = 0.26
    assert expected_calibration_error(confs, correct, n_bins=2) == pytest.approx(0.26)


def test_mce_picks_worst_bin():
    confs = np.concatenate([np.full(8, 0.2), np.full(2, 0.9)])
    correct = np.concatenate([
        np.array([1, 1, 1, 1, 0, 0, 0, 0]),
        np.array([1, 1]),
    ])
    # bin A gap = 0.3, bin B gap = 0.1 → MCE = 0.3
    assert maximum_calibration_error(confs, correct, n_bins=2) == pytest.approx(0.3)


def test_brier_score():
    confs = np.array([0.9, 0.1, 0.5])
    correct = np.array([1.0, 0, 1])
    # (0.01 + 0.01 + 0.25) / 3 = 0.09
    assert brier_score(confs, correct) == pytest.approx(0.09)


def test_errors_to_correctness():
    e = np.array([2.0, 5.0, 1.0, 10.0])
    c = errors_to_correctness(e, threshold_A=4.0)
    assert np.allclose(c, [1, 0, 1, 0])


def test_bootstrap_ci_brier():
    confs = np.full(100, 0.7)
    correct = np.array([1.0] * 70 + [0.0] * 30)  # 70% correct → perfectly calibrated
    point, lo, hi = bootstrap_ci(brier_score, confs, correct, n_boot=200, seed=1)
    # Brier for constant 0.7 vs 70% correct = 0.7*(1-0.7)^2 + 0.3*(0-0.7)^2 = 0.7*0.09+0.3*0.49=0.063+0.147=0.21
    assert 0.18 < point < 0.24
    assert lo < point <= hi


def test_reliability_table_shape():
    confs = np.random.RandomState(0).rand(200)
    correct = (confs > 0.5).astype(float)
    rt = reliability_table(confs, correct, n_bins=10)
    assert rt["n_bins"] == 10
    assert rt["n_total"] == 200
    assert len(rt["rows"]) == 10
    total = sum(r["count"] for r in rt["rows"])
    assert total == 200


def test_error_cases():
    with pytest.raises(ValueError, match="empty"):
        expected_calibration_error(np.array([]), np.array([]), n_bins=5)
    with pytest.raises(ValueError, match="shape mismatch"):
        expected_calibration_error(np.array([0.5]), np.array([1, 0]), n_bins=2)
    with pytest.raises(ValueError, match="n_bins"):
        expected_calibration_error(np.array([0.5]), np.array([1.0]), n_bins=0)
