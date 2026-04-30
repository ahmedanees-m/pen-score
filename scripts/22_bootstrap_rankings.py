"""Bootstrap CIs on PenScore Rankings.

Compute 1000x bootstrap confidence intervals on each editor's PenScore and
rank position (noise-model bootstrap: Gaussian perturbation σ=0.02 per axis,
seed=42, per axis_definitions composite.PenScore).

Method: for each bootstrap iteration, add N(0, 0.02) noise to each available
axis score, clip to [0,1], recompute weighted PenScore, rank all editors.
Report 2.5/97.5 percentile CI on both score and rank.

Outputs:
    ~/pen-stack/data/pen-score/scorecards/bootstrap_rankings.parquet
    ~/pen-stack/data/pen-score/scorecards/bootstrap_rankings.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/22_bootstrap_rankings.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SCORECARD = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
OUT_DIR   = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecards"
OUT_DIR.mkdir(parents=True, exist_ok=True)

AXES = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature"]

# Default weights (human_therapeutic_aav_insertion use case, axis_definitions v1.0.0)
WEIGHTS = {
    "S_DSB":    0.25,
    "S_Spec":   0.15,
    "S_Cargo":  0.20,
    "S_Deliv":  0.20,
    "S_Immuno": 0.10,
    "S_Prog":   0.05,
    "S_Mature": 0.05,
}

N_BOOT  = 1000
SEED    = 42
SIGMA   = 0.02    # axis measurement uncertainty (fixed)


def penscore_from_matrix(scores: np.ndarray,
                         weights: np.ndarray,
                         available: np.ndarray) -> np.ndarray:
    """
    scores:    (n_editors, n_axes)  - may contain NaN for missing axes
    weights:   (n_axes,)
    available: (n_editors, n_axes) bool - True where axis was computed

    Returns:   (n_editors,) PenScore, with renormalised weights per editor.
    """
    n_editors = scores.shape[0]
    out = np.full(n_editors, np.nan)
    for i in range(n_editors):
        mask = available[i]
        if mask.sum() == 0:
            continue
        w = weights[mask]
        s = scores[i, mask]
        out[i] = (w * s).sum() / w.sum()
    return out


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE Bootstrap CIs on Rankings")
    print(f"  N_BOOT={N_BOOT}  SEED={SEED}  σ={SIGMA}")
    print("=" * 60)

    df = pd.read_parquet(SCORECARD)
    editors = df["editor_id"].values
    n_editors = len(editors)
    n_axes = len(AXES)

    # Build arrays
    score_mat  = df[AXES].values.astype(float)          # (n_editors, n_axes)
    avail_mat  = ~np.isnan(score_mat)                   # (n_editors, n_axes)
    w_arr      = np.array([WEIGHTS[a] for a in AXES])   # (n_axes,)

    # Fill NaN with 0 for noise addition (we'll mask back with available)
    score_mat_filled = np.where(avail_mat, score_mat, 0.0)

    rng = np.random.default_rng(SEED)
    boot_scores = np.zeros((N_BOOT, n_editors))
    boot_ranks  = np.zeros((N_BOOT, n_editors))

    for b in range(N_BOOT):
        noise  = rng.normal(0.0, SIGMA, size=(n_editors, n_axes))
        noisy  = np.clip(score_mat_filled + noise * avail_mat, 0.0, 1.0)
        # Restore NaN for unavailable axes
        noisy  = np.where(avail_mat, noisy, np.nan)
        ps     = penscore_from_matrix(noisy, w_arr, avail_mat)
        boot_scores[b] = ps
        # Rank: 1=best; handle NaN by giving them last rank
        order  = np.argsort(np.where(np.isnan(ps), -np.inf, ps))[::-1]
        ranks  = np.empty(n_editors)
        ranks[order] = np.arange(1, n_editors + 1)
        boot_ranks[b] = ranks

    # Summarise
    out_df = pd.DataFrame({
        "editor_id":       editors,
        "penscore_observed": df["PenScore"].values,
        "penscore_mean":   np.nanmean(boot_scores, axis=0).round(4),
        "penscore_lo95":   np.nanquantile(boot_scores, 0.025, axis=0).round(4),
        "penscore_hi95":   np.nanquantile(boot_scores, 0.975, axis=0).round(4),
        "rank_observed":   df["PenScore"].rank(ascending=False, method="min").astype(int).values,
        "rank_mean":       np.nanmean(boot_ranks, axis=0).round(1),
        "rank_lo95":       np.nanquantile(boot_ranks, 0.025, axis=0).round(1),
        "rank_hi95":       np.nanquantile(boot_ranks, 0.975, axis=0).round(1),
        "n_axes":          df["n_axes"].values,
    }).sort_values("penscore_observed", ascending=False)

    # Report
    print(f"\n{'Editor':<22} {'Score':>7} {'Score 95%CI':>18}  "
          f"{'Rank':>5} {'Rank 95%CI':>14}  axes")
    print("-" * 75)
    for _, row in out_df.iterrows():
        ps = f"{row['penscore_observed']:.4f}" if pd.notna(row['penscore_observed']) else "  None"
        ci_s = f"[{row['penscore_lo95']:.4f},{row['penscore_hi95']:.4f}]"
        ci_r = f"[{row['rank_lo95']:.0f},{row['rank_hi95']:.0f}]"
        print(f"  {row['editor_id']:<20} {ps}  {ci_s}  {row['rank_observed']:>4}  {ci_r}  {int(row['n_axes'])}")

    out_df.to_parquet(OUT_DIR / "bootstrap_rankings.parquet", index=False)
    out_df.to_csv(OUT_DIR / "bootstrap_rankings.csv", index=False)
    print(f"\nWritten -> {OUT_DIR}/bootstrap_rankings.parquet (.csv)")
    print(f"  Bootstrap params: N={N_BOOT}, seed={SEED}, σ={SIGMA} per axis")


if __name__ == "__main__":
    main()
