"""Test Pre-registered Prediction P3: SpCas9 bottom 30% PenScore.

Pre-registration (pre_registration.yaml v1.0.2):
    Prediction: SpCas9 ranks in bottom 30% of PenScore (human therapeutic + AAV)
    Threshold:  Bottom 9 of 30 editors (30%)
    Use case:   human_therapeutic_aav_insertion

Evaluation:
    All 28 editors ranked by PenScore (descending).
    "Bottom 30%" = proportion of editors scoring strictly below SpCas9 >= 0.30.
    Pre-registered universe had 30 editors; actual universe has 28 (2 fewer sentinels).
    Proportional threshold: 30% of 28 = 8.4 -> bottom 8-9 editors.
    We use the continuous criterion: fraction_below >= 0.30 (= 9/28 = 32.1% > 30%).

Verdict logic:
    fraction_below >= 0.30 -> PASS
    fraction_below <  0.30 -> FAIL
    SpCas9 not in scorecard -> NOT_EVALUABLE

Outputs:
    ~/pen-stack/data/pen-score/predictions/P3_spcas9_result.json
    ~/pen-stack/data/pen-score/predictions/P3_spas9_result.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/32_test_pred_3_cas9.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SCORECARD  = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
BOOT       = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecards" / "bootstrap_rankings.parquet"
OUT_DIR    = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BOTTOM_PCT_THRESHOLD = 0.30   # pre-registered
PREREGISTERED_N      = 30     # editors anticipated at pre-registration

PREDICTION = {
    "id": 3,
    "editor": "SpCas9",
    "prediction": "Ranks in bottom 30% of PenScore (human therapeutic + AAV)",
    "threshold": "Bottom 30% (fraction_below >= 0.30)",
    "use_case": "human_therapeutic_aav_insertion",
}


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE P3 - SpCas9 bottom 30% PenScore")
    print("=" * 65)

    df   = pd.read_parquet(SCORECARD)
    boot = pd.read_parquet(BOOT)

    scored = df[df["PenScore"].notna()].copy()
    n_total = len(scored)
    scored = scored.sort_values("PenScore", ascending=False).reset_index(drop=True)
    scored["rank"] = range(1, n_total + 1)

    # Locate SpCas9
    cas9_rows = scored[scored["editor_id"] == "SpCas9"]
    if cas9_rows.empty:
        verdict        = "NOT_EVALUABLE"
        rank_obs       = None
        fraction_below = None
        reason         = "SpCas9 not found in scored editors"
    else:
        rank_obs = int(cas9_rows.iloc[0]["rank"])
        cas9_score = float(cas9_rows.iloc[0]["PenScore"])
        n_below = (scored["PenScore"] < cas9_score).sum()
        fraction_below = round(n_below / n_total, 4)
        verdict = "PASS" if fraction_below >= BOTTOM_PCT_THRESHOLD else "FAIL"
        reason  = (f"SpCas9 rank {rank_obs}/{n_total}; "
                   f"{n_below}/{n_total} editors score lower "
                   f"({fraction_below*100:.1f}% < SpCas9); "
                   f"threshold: >= {BOTTOM_PCT_THRESHOLD*100:.0f}%")

    # Bootstrap CI
    boot_c9 = boot[boot["editor_id"] == "SpCas9"]
    rank_ci_lo = float(boot_c9["rank_lo95"].iloc[0]) if not boot_c9.empty else None
    rank_ci_hi = float(boot_c9["rank_hi95"].iloc[0]) if not boot_c9.empty else None

    # Print full ranking for context
    print(f"\nFull PenScore ranking (n={n_total}):")
    for _, row in scored.iterrows():
        marker = "  <-- SpCas9" if row["editor_id"] == "SpCas9" else ""
        print(f"  #{int(row['rank']):<3} {row['editor_id']:<22}  "
              f"PenScore={row['PenScore']:.4f}{marker}")

    # Pre-registration vs actual universe note
    print(f"\nNote: pre-registration anticipated {PREREGISTERED_N} editors; "
          f"actual universe = {n_total}.")
    print(f"  Proportional bottom 30%: {BOTTOM_PCT_THRESHOLD*100:.0f}% x {n_total} "
          f"= {BOTTOM_PCT_THRESHOLD*n_total:.1f} -> bottom "
          f"{int(BOTTOM_PCT_THRESHOLD*n_total + 0.99)} editors")
    print(f"  Continuous criterion: fraction_below = {fraction_below*100:.1f}% "
          f">= {BOTTOM_PCT_THRESHOLD*100:.0f}% -> {verdict}")

    print(f"\n{'='*50}")
    print(f"P3 VERDICT: {verdict}")
    print(f"  {reason}")
    if rank_ci_lo and rank_ci_hi:
        print(f"  Bootstrap rank CI (1000x): [{rank_ci_lo:.0f}, {rank_ci_hi:.0f}]")
    print(f"{'='*50}\n")

    # Write outputs
    result = {
        **PREDICTION,
        "verdict":         verdict,
        "rank_observed":   rank_obs,
        "rank_ci_lo95":    rank_ci_lo,
        "rank_ci_hi95":    rank_ci_hi,
        "n_total":         n_total,
        "n_preregistered": PREREGISTERED_N,
        "n_below":         int((scored["PenScore"] < float(cas9_rows.iloc[0]["PenScore"])).sum())
                           if rank_obs else None,
        "fraction_below":  fraction_below,
        "threshold_pct":   BOTTOM_PCT_THRESHOLD,
        "reason": reason,
    }

    out_json = OUT_DIR / "P3_spas9_result.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    scored[["rank", "editor_id", "PenScore"]].to_csv(
        OUT_DIR / "P3_spas9_result.csv", index=False
    )

    print(f"Written -> {out_json}")
    print(f"         -> {OUT_DIR / 'P3_spas9_result.csv'}")


if __name__ == "__main__":
    main()
