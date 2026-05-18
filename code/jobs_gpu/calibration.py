"""Calibration metrics for continuous-error / scalar-confidence pairs.

Convention used here:
    confidence c in [0, 1]
    realized error e >= 0 (in Angstroms, per-residue)

We convert (c, e) into a binary-correctness frame by thresholding:
    correct := e <= threshold_A
That mirrors the classification calibration setting in Guo et al. 2017.

Provided:
    - expected_calibration_error
    - maximum_calibration_error
    - brier_score (squared error between confidence and correctness)
    - bootstrap_ci helper for any scalar metric

All functions accept numpy arrays and return python floats (for clean json /
markdown export).
"""
from __future__ import annotations
import numpy as np
from typing import Callable


def _validate(confidences: np.ndarray, correctness: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c = np.asarray(confidences, dtype=float)
    y = np.asarray(correctness, dtype=float)
    if c.size == 0 or y.size == 0:
        raise ValueError("empty input")
    if c.shape != y.shape:
        raise ValueError(f"shape mismatch: confidences {c.shape} vs correctness {y.shape}")
    return c, y


def expected_calibration_error(
    confidences: np.ndarray, correctness: np.ndarray, n_bins: int = 10
) -> float:
    """Binned ECE = sum_b (|B_b| / N) * |mean(conf in B_b) - mean(correct in B_b)|."""
    c, y = _validate(confidences, correctness)
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1, got {n_bins}")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = np.clip(np.digitize(c, edges) - 1, 0, n_bins - 1)
    n = c.size
    ece = 0.0
    for b in range(n_bins):
        mask = bins == b
        if mask.sum() == 0:
            continue
        ece += (mask.sum() / n) * abs(c[mask].mean() - y[mask].mean())
    return float(ece)


def maximum_calibration_error(
    confidences: np.ndarray, correctness: np.ndarray, n_bins: int = 10
) -> float:
    """MCE = max_b |mean(conf in B_b) - mean(correct in B_b)| (over non-empty bins)."""
    c, y = _validate(confidences, correctness)
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1, got {n_bins}")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = np.clip(np.digitize(c, edges) - 1, 0, n_bins - 1)
    gaps = []
    for b in range(n_bins):
        mask = bins == b
        if mask.sum() == 0:
            continue
        gaps.append(abs(c[mask].mean() - y[mask].mean()))
    return float(max(gaps)) if gaps else 0.0


def brier_score(confidences: np.ndarray, correctness: np.ndarray) -> float:
    """Brier score = mean((conf - correct)^2)."""
    c, y = _validate(confidences, correctness)
    return float(((c - y) ** 2).mean())


def errors_to_correctness(errors: np.ndarray, threshold_A: float) -> np.ndarray:
    """Binarize errors using a Å threshold. True if error <= threshold."""
    e = np.asarray(errors, dtype=float)
    return (e <= threshold_A).astype(float)


def bootstrap_ci(
    fn: Callable[[np.ndarray, np.ndarray], float],
    confidences: np.ndarray,
    correctness: np.ndarray,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Return (point_estimate, lower_ci, upper_ci) for a metric fn(c, y)."""
    c, y = _validate(confidences, correctness)
    point = fn(c, y)
    rng = np.random.default_rng(seed)
    boots = []
    n = c.size
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots.append(fn(c[idx], y[idx]))
    boots = np.array(boots)
    lo = float(np.percentile(boots, 100 * alpha / 2))
    hi = float(np.percentile(boots, 100 * (1 - alpha / 2)))
    return float(point), lo, hi


def reliability_table(
    confidences: np.ndarray, correctness: np.ndarray, n_bins: int = 10
) -> dict:
    """Per-bin diagnostics for plotting reliability diagrams."""
    c, y = _validate(confidences, correctness)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = np.clip(np.digitize(c, edges) - 1, 0, n_bins - 1)
    rows = []
    for b in range(n_bins):
        mask = bins == b
        rows.append({
            "bin": b,
            "lo": float(edges[b]),
            "hi": float(edges[b + 1]),
            "count": int(mask.sum()),
            "mean_confidence": float(c[mask].mean()) if mask.any() else float("nan"),
            "mean_correctness": float(y[mask].mean()) if mask.any() else float("nan"),
            "gap": float(abs(c[mask].mean() - y[mask].mean())) if mask.any() else float("nan"),
        })
    return {"n_bins": n_bins, "n_total": int(c.size), "rows": rows}
