"""Test Pre-registered Prediction P2: IS621 top-3 programmable DSB-free.

Pre-registration (pre_registration.yaml v1.0.2):
    Prediction: IS621 ranks in top-3 of programmable DSB-free systems
    Threshold:  Top 3 of ~5-10 editors in that subset
    Use case:   human_therapeutic_aav_insertion

Subset definition:
    "Programmable DSB-free systems" = editors with S_DSB >= 0.9 (DSB-free class)
    AND S_Prog = 1.0 (fully programmable - RNA-guided or CAST-type).
    Excludes site-specific recombinases (Cre, Bxb1, phiC31, Lambda_Int; S_Prog=0.0).

Verdict logic:
    IS621 rank within subset <= 3 -> PASS
    IS621 rank within subset  > 3 -> FAIL

Outputs:
    ~/pen-stack/data/pen-score/predictions/P2_is621_result.json
    ~/pen-stack/data/pen-score/predictions/P2_is621_result.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/31_test_pred_2_is621.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SCORECARD = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
BOOT      = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecards" / "bootstrap_rankings.parquet"
OUT_DIR   = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DSB_FREE_THRESHOLD = 0.85   # S_DSB >= 0.9 (bucket value)
PROG_THRESHOLD     = 0.99   # S_Prog = 1.0 (programmable)
TOP_N = 3

PREDICTION = {
    "id": 2,
    "editor": "IS621",
    "prediction": "Ranks in top-3 of programmable DSB-free systems",
    "threshold": f"Top {TOP_N} of ~5-10 editors in that subset",
    "use_case": "human_therapeutic_aav_insertion",
}


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE P2 - IS621 top-3 programmable DSB-free")
    print("=" * 65)

    df   = pd.read_parquet(SCORECARD)
    boot = pd.read_parquet(BOOT)

    # Define programmable DSB-free subset
    subset = df[(df["S_DSB"] >= DSB_FREE_THRESHOLD) &
                (df["S_Prog"] >= PROG_THRESHOLD)].copy()
    subset = subset.sort_values("PenScore", ascending=False, na_position="last")
    subset["subset_rank"] = range(1, len(subset) + 1)

    print(f"\nProgrammable DSB-free subset (S_DSB >= {DSB_FREE_THRESHOLD} AND S_Prog >= {PROG_THRESHOLD}):")
    print(f"  n = {len(subset)} editors\n")
    for _, row in subset.iterrows():
        marker = "  <-- IS621" if row["editor_id"] == "IS621" else ""
        print(f"  #{int(row['subset_rank']):<3} {row['editor_id']:<22}  "
              f"PenScore={row['PenScore']:.4f}  S_DSB={row['S_DSB']:.1f}  "
              f"S_Prog={row['S_Prog']:.1f}{marker}")

    # Evaluate P2
    is621_rows = subset[subset["editor_id"] == "IS621"]

    if is621_rows.empty:
        verdict  = "NOT_EVALUABLE"
        rank_obs = None
        reason   = "IS621 excluded from programmable DSB-free subset"
    else:
        rank_obs = int(is621_rows.iloc[0]["subset_rank"])
        verdict  = "PASS" if rank_obs <= TOP_N else "FAIL"
        reason   = f"IS621 rank #{rank_obs} in subset (threshold: top {TOP_N})"

    # Bootstrap CI for IS621
    boot_is = boot[boot["editor_id"] == "IS621"]
    rank_ci_lo = int(boot_is["rank_lo95"].iloc[0]) if not boot_is.empty else None
    rank_ci_hi = int(boot_is["rank_hi95"].iloc[0]) if not boot_is.empty else None

    print(f"\n{'='*50}")
    print(f"P2 VERDICT: {verdict}")
    print(f"  {reason}")
    if rank_ci_lo and rank_ci_hi:
        print(f"  Overall bootstrap rank CI (1000x): [{rank_ci_lo}, {rank_ci_hi}]")
    print(f"{'='*50}\n")

    # Write outputs
    result = {
        **PREDICTION,
        "verdict":      verdict,
        "rank_observed": rank_obs,
        "rank_ci_lo95":  rank_ci_lo,
        "rank_ci_hi95":  rank_ci_hi,
        "subset_size":   len(subset),
        "subset_def":    (f"S_DSB >= {DSB_FREE_THRESHOLD} (DSB-free) AND "
                          f"S_Prog >= {PROG_THRESHOLD} (programmable)"),
        "reason": reason,
        "subset_rankings": [
            {
                "subset_rank": int(r["subset_rank"]),
                "editor_id":   r["editor_id"],
                "PenScore":    float(r["PenScore"]) if pd.notna(r["PenScore"]) else None,
                "S_Prog":      float(r["S_Prog"]),
            }
            for _, r in subset.iterrows()
        ],
    }

    out_json = OUT_DIR / "P2_is621_result.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")

    subset[["subset_rank", "editor_id", "PenScore", "S_DSB", "S_Prog"]].to_csv(
        OUT_DIR / "P2_is621_result.csv", index=False
    )

    print(f"Written -> {out_json}")
    print(f"         -> {OUT_DIR / 'P2_is621_result.csv'}")


if __name__ == "__main__":
    main()
