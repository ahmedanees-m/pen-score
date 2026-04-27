"""PenScore composite function.

Weighted linear combination of 8 axis scores.  Weights are use-case specific
and loaded from use_case_profiles.yaml.  Missing axes (None) are excluded from
the weighted sum; the remaining weights are renormalised so they still sum to 1.

Default weights (human_therapeutic_aav_insertion use case, v0.1.1):
    w_DSB    = 0.24
    w_Spec   = 0.14
    w_Cargo  = 0.19
    w_Deliv  = 0.19
    w_Immuno = 0.09
    w_Prog   = 0.05
    w_Mature = 0.05
    w_Energy = 0.05

All weights sum to 1.0 over the 8 axes.  When an axis is None, the weight
is redistributed proportionally across the available axes.
S_Energy added in v0.1.1; weights were rescaled from v0.1.0 x 0.95.
"""

from __future__ import annotations

_AXES = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature", "S_Energy"]


def compute_pen_score(
    axis_scores: dict[str, float | None],
    weights: dict[str, float],
    exclude_axes: list[str] | None = None,
) -> tuple[float | None, list[str]]:
    """Compute weighted PenScore.

    Parameters
    ----------
    axis_scores:
        Dict mapping axis name -> score in [0, 1] or None.
    weights:
        Dict mapping axis name -> weight.  Should sum to ~1.0.
    exclude_axes:
        Axis names to exclude entirely from the composite.  Their weights are
        redistributed proportionally across remaining available axes.
        Excluded axes are NOT reported in the missing list.

    Returns
    -------
    (pen_score, missing_axes) where pen_score is in [0, 1] or None if all
    axes are missing or excluded.
    """
    excluded: frozenset[str] = frozenset(exclude_axes or [])
    available_w = 0.0
    weighted_sum = 0.0
    missing: list[str] = []

    for ax in _AXES:
        if ax in excluded:
            continue  # explicitly excluded; weight redistributed to remaining axes
        val = axis_scores.get(ax)
        w = weights.get(ax, 0.0)
        if val is None:
            missing.append(ax)
        else:
            available_w += w
            weighted_sum += w * val

    if available_w == 0.0:
        return None, missing

    # Renormalise: missing + excluded weights redistribute to available axes
    pen_score = weighted_sum / available_w
    return round(pen_score, 4), missing
