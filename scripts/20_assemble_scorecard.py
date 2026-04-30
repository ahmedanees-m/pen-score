"""Assemble PEN-SCORE Composite Scorecard.

Reads all 8 axis parquet files (v0.1.1), merges on editor_id, and computes the
composite PenScore per axis_definitions.yaml v1.0.1.

Composite formula (axis_definitions composite.PenScore):
    PenScore = sum(w_i * S_i for i in available_axes)
             / sum(w_i          for i in available_axes)
    Missing axes (NaN) are excluded from both numerator and denominator.
    Weights are renormalised automatically over available axes per editor.

Default weights (human_therapeutic_aav_insertion use case, v0.1.1 8-axis):
    S_DSB    0.24
    S_Spec   0.14
    S_Cargo  0.19
    S_Deliv  0.19
    S_Immuno 0.09
    S_Prog   0.05
    S_Mature 0.05
    S_Energy 0.05

Inputs  (auto-detected from ~/pen-stack/data/pen-score/axes/*/):
    cargo_scores.parquet
    deliv_scores.parquet
    dsb_scores.parquet
    energy_scores.parquet
    immuno_scores.parquet
    mature_scores.parquet
    prog_scores.parquet
    spec_scores.parquet

Output:
    ~/pen-stack/data/pen-score/scorecard.parquet
    ~/pen-stack/data/pen-score/scorecard.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/20_assemble_scorecard.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO     = Path(__file__).resolve().parent.parent
DATA     = REPO / "pen_score" / "data"
AXES_DIR = Path.home() / "pen-stack" / "data" / "pen-score" / "axes"
OUT_DIR  = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecards"
OUT_DIR.mkdir(parents=True, exist_ok=True)
# Also write to the flat pen-score dir for backwards compat
FLAT_DIR = Path.home() / "pen-stack" / "data" / "pen-score"

# Axis parquet filenames and their score columns  (v0.1.1 - 8 axes)
AXIS_FILES: dict[str, tuple[str, float]] = {
    "dsb_scores.parquet":    ("S_DSB",    0.24),
    "spec_scores.parquet":   ("S_Spec",   0.14),
    "cargo_scores.parquet":  ("S_Cargo",  0.19),
    "deliv_scores.parquet":  ("S_Deliv",  0.19),
    "immuno_scores.parquet": ("S_Immuno", 0.09),
    "prog_scores.parquet":   ("S_Prog",   0.05),
    "mature_scores.parquet": ("S_Mature", 0.05),
    "energy_scores.parquet": ("S_Energy", 0.05),
}

# Axis subdirectory names
AXIS_SUBDIR: dict[str, str] = {
    "dsb_scores.parquet":    "dsb",
    "spec_scores.parquet":   "spec",
    "cargo_scores.parquet":  "cargo",
    "deliv_scores.parquet":  "deliv",
    "immuno_scores.parquet": "immuno",
    "prog_scores.parquet":   "prog",
    "mature_scores.parquet": "mature",
    "energy_scores.parquet": "energy",
}


def compute_penscore(row: pd.Series,
                     score_cols: list[str],
                     weights: dict[str, float]) -> float | None:
    """Weighted average of available (non-NaN) axis scores."""
    num = 0.0
    den = 0.0
    for col in score_cols:
        v = row[col]
        if pd.notna(v):
            w = weights[col]
            num += w * float(v)
            den += w
    if den == 0.0:
        return None
    return round(num / den, 4)


def _load_editor_metadata() -> pd.DataFrame:
    """Pull organism, mechanism_bucket, year_discovered, primary_reference from YAML."""
    import yaml
    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    rows = []
    for ed in universe["editors"]:
        rows.append({
            "editor_id":        ed["id"],
            "organism":         ed.get("organism", ""),
            "mechanism_bucket": ed.get("mechanism_bucket", ""),
            "year_discovered":  ed.get("year_discovered"),
            "primary_reference": ed.get("primary_reference", ""),
        })
    return pd.DataFrame(rows)


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE Composite Scorecard Assembly (v0.1.1, 8 axes)")
    print(f"  Axes dir: {AXES_DIR}")
    print(f"  Output:   {OUT_DIR}")
    print("=" * 65)

    # Load and merge all axis parquets
    merged: pd.DataFrame | None = None
    score_cols: list[str] = []
    weights: dict[str, float] = {}

    for fname, (score_col, weight) in AXIS_FILES.items():
        subdir = AXIS_SUBDIR[fname]
        path = AXES_DIR / subdir / fname
        if not path.exists():
            print(f"  WARNING: {path} not found - axis {score_col} excluded")
            continue

        df_ax = pd.read_parquet(path)[["editor_id", score_col]]
        n_ok = df_ax[score_col].notna().sum()
        print(f"  Loaded {score_col:<10}  w={weight:.2f}  "
              f"({n_ok}/{len(df_ax)} computed)  <- {path.name}")

        if merged is None:
            # Start with a full editor list from the first axis
            merged = pd.read_parquet(path)[["editor_id", "canonical_accession"]]
            merged = merged.merge(df_ax, on="editor_id", how="left")
        else:
            merged = merged.merge(df_ax, on="editor_id", how="left")

        score_cols.append(score_col)
        weights[score_col] = weight

    if merged is None or merged.empty:
        raise RuntimeError("No axis parquets loaded - nothing to assemble.")

    print(f"\n  {len(merged)} editors  |  {len(score_cols)} axes loaded\n")

    # Compute composite PenScore
    merged["PenScore"] = merged.apply(
        lambda row: compute_penscore(row, score_cols, weights), axis=1
    )

    # Count available axes per editor
    merged["n_axes"] = merged[score_cols].notna().sum(axis=1)

    # Report
    sorted_df = merged.sort_values("PenScore", ascending=False, na_position="last")

    print(f"{'Rank':<5} {'Editor':<22} {'PenScore':<10} "
          + "  ".join(f"{c[2:]:<7}" for c in score_cols)
          + "  n_ax")
    print("-" * (5 + 22 + 10 + 9 * len(score_cols) + 5))

    for rank, (_, row) in enumerate(sorted_df.iterrows(), 1):
        ps = f"{row['PenScore']:.4f}" if pd.notna(row["PenScore"]) else " None"
        axis_vals = "  ".join(
            f"{row[c]:.4f}" if pd.notna(row[c]) else "  None"
            for c in score_cols
        )
        print(f"{rank:<5} {row['editor_id']:<22} {ps:<10}  {axis_vals}  {int(row['n_axes'])}")

    # Merge editor metadata from YAML
    meta = _load_editor_metadata()
    merged = meta.merge(merged, on="editor_id", how="right")

    # Write outputs
    col_order = (["editor_id", "canonical_accession",
                  "organism", "mechanism_bucket", "year_discovered", "primary_reference",
                  "PenScore", "n_axes"]
                 + score_cols)
    out_df = merged[col_order]

    # Primary output: scorecards/public_scorecard.parquet
    out_df.to_parquet(OUT_DIR / "public_scorecard.parquet", index=False)
    out_df.to_csv(OUT_DIR / "public_scorecard.csv", index=False)
    # Backwards-compat flat copy
    out_df.to_parquet(FLAT_DIR / "scorecard.parquet", index=False)
    out_df.to_csv(FLAT_DIR / "scorecard.csv", index=False)

    n_full  = (out_df["n_axes"] == len(score_cols)).sum()
    n_part  = ((out_df["n_axes"] > 0) & (out_df["n_axes"] < len(score_cols))).sum()
    n_none  = (out_df["n_axes"] == 0).sum()
    top     = sorted_df.iloc[0]

    print(f"\nWritten -> {OUT_DIR}/public_scorecard.parquet (.csv)  [primary]")
    print(f"         -> {FLAT_DIR}/scorecard.parquet (.csv)  [compat]")
    print(f"  Full ({len(score_cols)}/{len(score_cols)} axes): {n_full}  |  "
          f"Partial: {n_part}  |  No score: {n_none}")
    print(f"  Top editor: {top['editor_id']}  PenScore={top['PenScore']:.4f}")
    print(f"  Weights used: "
          + ", ".join(f"{c}={weights[c]}" for c in score_cols))


if __name__ == "__main__":
    main()
