"""Editor ranking and selection decision support."""

from __future__ import annotations

import pandas as pd


def rank_editors(
    scorecard: pd.DataFrame,
    use_case: str,
    weights: dict[str, float],
    filters: dict[str, float] | None = None,
    top_k: int | None = None,
) -> pd.DataFrame:
    """Rank editors in the scorecard by PenScore for a given use case.

    Parameters
    ----------
    scorecard:
        DataFrame with columns: editor_id, canonical_accession, S_DSB, S_Spec,
        S_Cargo, S_Deliv, S_Immuno, S_Prog, S_Mature, PenScore.
    use_case:
        String label for the use case (for display).
    weights:
        Axis weights dict.
    filters:
        Optional minimum value per axis (e.g. {'S_Deliv': 0.5} to restrict
        to AAV-deliverable editors only).
    top_k:
        Return only the top-k editors.

    Returns
    -------
    Sorted DataFrame with added columns: use_case, rank, reasoning.
    """
    df = scorecard.copy()

    # Apply filters
    if filters:
        for col, min_val in filters.items():
            if col in df.columns:
                df = df[df[col].fillna(0.0) >= min_val]

    # Sort by PenScore descending
    df = df.sort_values("PenScore", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    df["use_case"] = use_case

    # Build human-readable reasoning column
    df["reasoning"] = df.apply(lambda row: _build_reasoning(row, weights), axis=1)

    if top_k is not None:
        df = df.head(top_k)

    return df


def _build_reasoning(row: pd.Series, weights: dict[str, float]) -> list[str]:
    """Build axis-by-axis reasoning for an editor ranking."""
    reasons: list[str] = []
    axis_labels = {
        "S_DSB": "DSB avoidance",
        "S_Spec": "specificity",
        "S_Cargo": "cargo capacity",
        "S_Deliv": "AAV deliverability",
        "S_Immuno": "low immunogenicity",
        "S_Prog": "programmability",
        "S_Mature": "therapeutic maturity",
    }
    for ax, label in axis_labels.items():
        val = row.get(ax)
        w = weights.get(ax, 0.0)
        if val is None:
            reasons.append(f"{label}: N/A (weight {w:.2f})")
        elif val >= 0.8:
            reasons.append(f"Excellent {label}: {val:.2f} (weight {w:.2f})")
        elif val >= 0.5:
            reasons.append(f"Moderate {label}: {val:.2f} (weight {w:.2f})")
        else:
            reasons.append(f"Poor {label}: {val:.2f} (weight {w:.2f})")
    return reasons
