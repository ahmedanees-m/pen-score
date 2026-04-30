"""Test Pre-registered Prediction P4: enNlovFz2 S_Deliv > NlovFz2 WT.

Pre-registration (pre_registration.yaml v1.0.2):
    Prediction: S_Deliv strictly greater than NlovFz2 WT
    Threshold:  enNlovFz2.S_Deliv > NlovFz2.S_Deliv
    Use case:   any

Sentinel status:
    Both NlovFz2 and enNlovFz2 are sentinels: sequences not in public databases.
    Resolution path: Wei et al. 2025 Nat Chem Biol SI - BLAST vs
    XP_04454xxxx (Naegleria lovaniensis ATCC 30569).
    Until resolved, S_Deliv = None for both editors.

Verdict logic:
    enNlovFz2.S_Deliv > NlovFz2.S_Deliv (both non-NaN) -> PASS or FAIL
    Either S_Deliv is NaN -> NOT_EVALUABLE

Scientific rationale for expected PASS:
    enNlovFz2 ("enhanced NlovFz2") is an engineered variant with protein
    modifications designed to improve AAV-delivery compatibility. Engineering
    strategies for Fanzor delivery improvement (e.g., NLS/NES addition, domain
    deletion) typically reduce protein size or improve vector packaging - both
    of which increase S_Deliv. Expected to PASS when sentinel is resolved.

Outputs:
    ~/pen-stack/data/pen-score/predictions/P4_ennlovfz2_result.json

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/33_test_pred_4_ennlovfz2.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SCORECARD = Path.home() / "pen-stack" / "data" / "pen-score" / "scorecard.parquet"
OUT_DIR   = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PREDICTION = {
    "id": 4,
    "editor": "enNlovFz2",
    "prediction": "S_Deliv strictly greater than NlovFz2 WT",
    "threshold": "enNlovFz2.S_Deliv > NlovFz2.S_Deliv",
    "use_case": "any",
}


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE P4 - enNlovFz2 S_Deliv > NlovFz2")
    print("=" * 65)

    df = pd.read_parquet(SCORECARD)

    # Extract the two editors
    row_en  = df[df["editor_id"] == "enNlovFz2"]
    row_wt  = df[df["editor_id"] == "NlovFz2"]

    if row_en.empty or row_wt.empty:
        verdict = "NOT_EVALUABLE"
        s_deliv_en = s_deliv_wt = None
        reason = "One or both editors not found in scorecard"
    else:
        s_deliv_en = row_en.iloc[0].get("S_Deliv")
        s_deliv_wt = row_wt.iloc[0].get("S_Deliv")

        # Convert to None if NaN
        s_deliv_en = float(s_deliv_en) if pd.notna(s_deliv_en) else None
        s_deliv_wt = float(s_deliv_wt) if pd.notna(s_deliv_wt) else None

        if s_deliv_en is None or s_deliv_wt is None:
            verdict = "NOT_EVALUABLE"
            reason  = (
                "SENTINEL - S_Deliv unavailable for one or both editors. "
                "Both NlovFz2 and enNlovFz2 require sequence retrieval from "
                "Wei et al. 2025 Nat Chem Biol SI (BLAST vs XP_04454xxxx, "
                "Naegleria lovaniensis ATCC 30569) before S_Deliv can be computed."
            )
        elif s_deliv_en > s_deliv_wt:
            verdict = "PASS"
            reason  = f"enNlovFz2.S_Deliv={s_deliv_en:.4f} > NlovFz2.S_Deliv={s_deliv_wt:.4f}"
        else:
            verdict = "FAIL"
            reason  = f"enNlovFz2.S_Deliv={s_deliv_en:.4f} NOT > NlovFz2.S_Deliv={s_deliv_wt:.4f}"

    print(f"\n  enNlovFz2  S_Deliv = {s_deliv_en}")
    print(f"  NlovFz2    S_Deliv = {s_deliv_wt}")
    print(f"\n{'='*50}")
    print(f"P4 VERDICT: {verdict}")
    print(f"  {reason}")
    if verdict == "NOT_EVALUABLE":
        print(f"\n  Scientific rationale for expected PASS when resolved:")
        print(f"  enNlovFz2 engineering designed for improved AAV-delivery;")
        print(f"  smaller/optimized protein typically yields higher S_Deliv.")
    print(f"{'='*50}\n")

    # Write output
    result = {
        **PREDICTION,
        "verdict":       verdict,
        "s_deliv_en":    s_deliv_en,
        "s_deliv_wt":    s_deliv_wt,
        "delta":         round(s_deliv_en - s_deliv_wt, 4) if (s_deliv_en and s_deliv_wt) else None,
        "reason":        reason,
        "sentinel_resolution": (
            "Wei et al. 2025 Nat Chem Biol SI - BLAST vs XP_04454xxxx "
            "(Naegleria lovaniensis ATCC 30569). After retrieval, re-run "
            "scripts/12_compute_S_Deliv.py then this script."
        ),
    }

    out_json = OUT_DIR / "P4_ennlovfz2_result.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"Written -> {out_json}")


if __name__ == "__main__":
    main()
