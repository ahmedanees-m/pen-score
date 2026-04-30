"""Test Pre-registered Prediction P1: evoCAST top-5 DSB-free integrases.

Pre-registration (pre_registration.yaml v1.0.2):
    Prediction: evoCAST ranks in top-5 of AAV-deliverable DSB-free integrases
    Threshold:  Top 5 of ~10-15 editors in that subset
    Use case:   human_therapeutic_aav_insertion

Subset definition:
    "DSB-free integrases" = editors with S_DSB >= 0.9 (mechanism_bucket contains
    DSB_FREE_TRANSEST_RECOMBINASE). All such editors are AAV-compatible by mechanism
    class (Tn7-like and serine/tyrosine recombinases all fit within AAV cargo limits).
    Subset is NOT filtered by S_Deliv to avoid excluding evoCAST (sentinel).

Verdict logic:
    evoCAST rank within subset <= 5 -> PASS
    evoCAST rank within subset  > 5 -> FAIL
    evoCAST excluded from subset     -> NOT EVALUABLE

Note on evoCAST S_Deliv sentinel:
    evoCAST sequence not in public databases (Witte 2025 Science SI Table S1 required).
    S_Deliv and S_Immuno are missing; PenScore computed from 5/7 axes (renormalized).
    P1 result is marked "provisional" pending sentinel resolution. If S_Deliv < 0.5
    when resolved, evoCAST may be excluded from the AAV-deliverable subset.

Outputs:
    ~/pen-stack/data/pen-score/predictions/P1_evocast_result.json
    ~/pen-stack/data/pen-score/predictions/P1_evocast_result.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/30_test_pred_1_evocast.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SCORECARD = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
BOOT      = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecards" / "bootstrap_rankings.parquet"
OUT_DIR   = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DSB_FREE_THRESHOLD = 0.85   # S_DSB >= this -> DSB-free (bucket value is 0.9)
TOP_N = 5                   # pre-registered top-5 threshold

PREDICTION = {
    "id": 1,
    "editor": "evoCAST",
    "prediction": "Ranks in top-5 of AAV-deliverable DSB-free integrases",
    "threshold": f"Top {TOP_N} of ~10-15 editors in that subset",
    "use_case": "human_therapeutic_aav_insertion",
}


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE P1 - evoCAST top-5 DSB-free integrases")
    print("=" * 65)

    df   = pd.read_parquet(SCORECARD)
    boot = pd.read_parquet(BOOT)

    # Define DSB-free integrase subset
    subset = df[df["S_DSB"] >= DSB_FREE_THRESHOLD].copy()
    subset = subset.sort_values("PenScore", ascending=False, na_position="last")
    subset["subset_rank"] = range(1, len(subset) + 1)

    print(f"\nDSB-free integrase subset (S_DSB >= {DSB_FREE_THRESHOLD}):")
    print(f"  n = {len(subset)} editors\n")
    for _, row in subset.iterrows():
        marker = "  <-- evoCAST" if row["editor_id"] == "evoCAST" else ""
        deliv_str = f"{row['S_Deliv']:.4f}" if pd.notna(row.get("S_Deliv")) else "  None"
        print(f"  #{int(row['subset_rank']):<3} {row['editor_id']:<22}  "
              f"PenScore={row['PenScore']:.4f}  n_axes={int(row['n_axes'])}  "
              f"S_Deliv={deliv_str}{marker}")

    # Evaluate P1
    evo_rows = subset[subset["editor_id"] == "evoCAST"]

    if evo_rows.empty:
        verdict   = "NOT_EVALUABLE"
        rank_obs  = None
        reason    = "evoCAST excluded from DSB-free subset (S_DSB < threshold)"
    else:
        rank_obs = int(evo_rows.iloc[0]["subset_rank"])
        if rank_obs <= TOP_N:
            verdict = "PASS"
            reason  = f"evoCAST rank #{rank_obs} in subset (threshold: top {TOP_N})"
        else:
            verdict = "FAIL"
            reason  = f"evoCAST rank #{rank_obs} in subset (threshold: top {TOP_N})"

    # Bootstrap CI for evoCAST rank overall (for context)
    boot_evo = boot[boot["editor_id"] == "evoCAST"]
    rank_ci_lo = int(boot_evo["rank_lo95"].iloc[0]) if not boot_evo.empty else None
    rank_ci_hi = int(boot_evo["rank_hi95"].iloc[0]) if not boot_evo.empty else None

    # Provisional flag (sentinel)
    provisional = pd.isna(evo_rows.iloc[0]["S_Deliv"]) if not evo_rows.empty else True

    print(f"\n{'='*50}")
    print(f"P1 VERDICT: {verdict}")
    print(f"  {reason}")
    if provisional:
        print(f"  *** PROVISIONAL: S_Deliv=None (sentinel). Result may change")
        print(f"      if resolved S_Deliv < 0.5 excludes evoCAST from subset. ***")
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
        "subset_def":    f"S_DSB >= {DSB_FREE_THRESHOLD} (DSB_FREE_TRANSEST_RECOMBINASE)",
        "provisional":   provisional,
        "provisional_reason": (
            "S_Deliv=None (evoCAST sentinel; Witte 2025 Science SI Table S1 required)"
            if provisional else None
        ),
        "reason": reason,
        "subset_rankings": [
            {
                "subset_rank": int(r["subset_rank"]),
                "editor_id":   r["editor_id"],
                "PenScore":    float(r["PenScore"]) if pd.notna(r["PenScore"]) else None,
                "n_axes":      int(r["n_axes"]),
            }
            for _, r in subset.iterrows()
        ],
    }

    out_json = OUT_DIR / "P1_evocast_result.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")

    out_csv = OUT_DIR / "P1_evocast_result.csv"
    subset[["subset_rank", "editor_id", "PenScore", "n_axes", "S_Deliv"]].to_csv(
        out_csv, index=False
    )

    print(f"Written -> {out_json}")
    print(f"         -> {out_csv}")


if __name__ == "__main__":
    main()
