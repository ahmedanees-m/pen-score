"""Bootstrap confidence intervals on scorecard rankings.

Same methodology as GENOME-ATLAS and MECH-CLASS:
    1000x bootstrap resampling on the editor universe
    Seed = 42 (fixed across all PEN-STACK papers)
    Report 2.5th and 97.5th percentiles as the 95% CI

Applied to:
    - Individual axis scores
    - PenScore rankings
    - Pre-registered prediction tests
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_SEED = 42
_N_BOOTSTRAP = 1000


def bootstrap_ranking_ci(
    scorecard: pd.DataFrame,
    score_col: str = "PenScore",
    n_bootstrap: int = _N_BOOTSTRAP,
    seed: int = _SEED,
    ci: float = 0.95,
) -> pd.DataFrame:
    """Compute bootstrap CI on rankings for each editor.

    Parameters
    ----------
    scorecard:
        DataFrame with at least ``editor_id`` and ``score_col``.
    score_col:
        Column to bootstrap ranks on.
    n_bootstrap:
        Number of bootstrap iterations.
    seed:
        Random seed (fixed to 42 per PEN-STACK convention).
    ci:
        Confidence interval width (0.95 = 95% CI).

    Returns
    -------
    DataFrame with columns: editor_id, rank_mean, rank_ci_lower, rank_ci_upper,
    score_mean, score_ci_lower, score_ci_upper.
    """
    rng = np.random.default_rng(seed)
    n_editors = len(scorecard)
    editors = scorecard["editor_id"].values
    scores: np.ndarray = np.asarray(scorecard[score_col].fillna(0.0), dtype=np.float64)

    rank_matrix = np.zeros((n_bootstrap, n_editors), dtype=np.float32)

    for b in range(n_bootstrap):
        # Resample editors with replacement to form a bootstrap universe
        idx = rng.integers(0, n_editors, size=n_editors)
        boot_scores = scores[idx]
        # For each original editor i, their rank = how many bootstrap scores
        # strictly exceed their own score + 1.  This keeps each editor's
        # identity fixed while letting the reference set vary, giving a
        # meaningful per-editor rank distribution.
        for i in range(n_editors):
            rank_matrix[b, i] = float(np.sum(boot_scores > scores[i]) + 1)

    alpha = (1.0 - ci) / 2.0
    rows = []
    for i, eid in enumerate(editors):
        rows.append(
            {
                "editor_id": eid,
                "rank_mean": float(rank_matrix[:, i].mean()),
                "rank_ci_lower": float(np.quantile(rank_matrix[:, i], alpha)),
                "rank_ci_upper": float(np.quantile(rank_matrix[:, i], 1 - alpha)),
                # Score statistics: point estimates (no within-editor variance
                # is available from a single-measurement scorecard)
                "score_mean": float(scores[i]),
                "score_ci_lower": float(scores[i]),
                "score_ci_upper": float(scores[i]),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_axis_ci(
    values: np.ndarray,
    n_bootstrap: int = _N_BOOTSTRAP,
    seed: int = _SEED,
    ci: float = 0.95,
) -> tuple[float, float, float]:
    """Bootstrap CI on a single axis value array.

    Returns
    -------
    (mean, ci_lower, ci_upper)
    """
    rng = np.random.default_rng(seed)
    alpha = (1.0 - ci) / 2.0
    boot_means = np.array(
        [rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_bootstrap)]
    )
    lo = float(np.quantile(boot_means, alpha))
    hi = float(np.quantile(boot_means, 1 - alpha))
    return float(boot_means.mean()), lo, hi
