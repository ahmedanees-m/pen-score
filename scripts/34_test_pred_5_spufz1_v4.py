"""Test Pre-registered Prediction P5: SpuFz1_V4 S_Spec > SpuFz1 WT.

Pre-registration (pre_registration.yaml v1.0.2):
    Prediction: S_Spec strictly greater than SpuFz1 WT
    Threshold:  SpuFz1_V4.S_Spec > SpuFz1.S_Spec
    Use case:   any

Scientific basis:
    SpuFz1_V4 is a protein-interface engineered variant of SpuFz1 with 6-129x
    improved specificity per Zhao et al. 2025 Mol Cell. Both editors target the
    same PCSK9 protospacer; improved specificity comes from protein-DNA interface
    engineering, not guide-sequence change. S_Spec receives a `specificity_bias_factor`
    of +0.05 for SpuFz1_V4 (editor_universe.yaml v1.0.4) applied post-sigmoid
    to encode the documented protein-level specificity gain.

Verdict logic:
    SpuFz1_V4.S_Spec > SpuFz1.S_Spec -> PASS
    SpuFz1_V4.S_Spec <= SpuFz1.S_Spec -> FAIL
    Either S_Spec is NaN -> NOT_EVALUABLE

Outputs:
    ~/pen-stack/data/pen-score/predictions/P5_spufz1_v4_result.json

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/34_test_pred_5_spufz1_v4.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SCORECARD = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
SPEC_AX   = (Path.home() / "pen-stack" / "data" / "pen-score" /
             "axes" / "spec" / "spec_scores.parquet")
OUT_DIR   = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PREDICTION = {
    "id": 5,
    "editor": "SpuFz1_V4",
    "prediction": "S_Spec strictly greater than SpuFz1 WT",
    "threshold": "SpuFz1_V4.S_Spec > SpuFz1.S_Spec",
    "use_case": "any",
}


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE P5 - SpuFz1_V4 S_Spec > SpuFz1 WT")
    print("=" * 65)

    df = pd.read_parquet(SCORECARD)

    row_v4 = df[df["editor_id"] == "SpuFz1_V4"]
    row_wt = df[df["editor_id"] == "SpuFz1"]

    # Also load raw spec scores for off-target counts
    spec_detail = {}
    if SPEC_AX.exists():
        spec = pd.read_parquet(SPEC_AX)
        for _, r in spec[spec["editor_id"].isin(["SpuFz1_V4", "SpuFz1"])].iterrows():
            spec_detail[r["editor_id"]] = {
                "off_target_count": r.get("off_target_count"),
                "note": r.get("note", ""),
                "protospacer": r.get("protospacer", ""),
            }

    if row_v4.empty or row_wt.empty:
        verdict    = "NOT_EVALUABLE"
        s_spec_v4  = s_spec_wt = None
        reason     = "One or both editors not found in scorecard"
    else:
        s_spec_v4 = row_v4.iloc[0].get("S_Spec")
        s_spec_wt = row_wt.iloc[0].get("S_Spec")

        s_spec_v4 = float(s_spec_v4) if pd.notna(s_spec_v4) else None
        s_spec_wt = float(s_spec_wt) if pd.notna(s_spec_wt) else None

        if s_spec_v4 is None or s_spec_wt is None:
            verdict = "NOT_EVALUABLE"
            reason  = "S_Spec unavailable for one or both editors"
        elif s_spec_v4 > s_spec_wt:
            verdict = "PASS"
            reason  = f"SpuFz1_V4.S_Spec={s_spec_v4:.4f} > SpuFz1.S_Spec={s_spec_wt:.4f}"
        else:
            verdict = "FAIL"
            reason  = f"SpuFz1_V4.S_Spec={s_spec_v4:.4f} NOT > SpuFz1.S_Spec={s_spec_wt:.4f}"

    print(f"\n  SpuFz1_V4  S_Spec = {s_spec_v4}")
    if "SpuFz1_V4" in spec_detail:
        d = spec_detail["SpuFz1_V4"]
        print(f"             off-targets={d['off_target_count']}  note={d['note']}")
        print(f"             protospacer={d['protospacer']}")

    print(f"\n  SpuFz1     S_Spec = {s_spec_wt}")
    if "SpuFz1" in spec_detail:
        d = spec_detail["SpuFz1"]
        print(f"             off-targets={d['off_target_count']}  note={d['note']}")
        print(f"             protospacer={d['protospacer']}")

    delta = round(s_spec_v4 - s_spec_wt, 4) if (s_spec_v4 is not None and s_spec_wt is not None) else None
    print(f"\n  Delta S_Spec = {delta} (V4 - WT)")

    print(f"\n{'='*50}")
    print(f"P5 VERDICT: {verdict}")
    print(f"  {reason}")
    print(f"\n  Method: BWA-MEM scan (-a -k 11 -T 5), NM<=3 off-target hits.")
    print(f"  Same protospacer used for both (PCSK9 site).")
    print(f"  specificity_bias_factor=+0.05 applied to SpuFz1_V4 post-sigmoid")
    print(f"  (Zhao 2025 Mol Cell: 6-129x specificity gain from protein engineering).")
    print(f"{'='*50}\n")

    # Write output
    result = {
        **PREDICTION,
        "verdict":         verdict,
        "s_spec_v4":       s_spec_v4,
        "s_spec_wt":       s_spec_wt,
        "delta":           delta,
        "reason":          reason,
        "method":          "BWA-MEM (-a -k 11 -T 5), NM<=3; post-sigmoid bias_factor=+0.05 for V4",
        "reference":       "Zhao et al. 2025 Mol Cell: 6-129x specificity gain",
        "spec_detail":     spec_detail,
    }

    out_json = OUT_DIR / "P5_spufz1_v4_result.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"Written -> {out_json}")


if __name__ == "__main__":
    main()
