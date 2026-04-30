"""Summarise all pre-registered prediction results.

Reads P1-P5 JSON outputs from ~/pen-stack/data/pen-score/predictions/
and produces a single summary table.

Outputs:
    ~/pen-stack/data/pen-score/predictions/prediction_summary.json
    ~/pen-stack/data/pen-score/predictions/prediction_summary.csv

Run:
    python3 ~/pen-stack/code/repos/pen-score/scripts/35_summarise_predictions.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

OUT_DIR = Path.home() / "pen-stack" / "data" / "pen-score" / "predictions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RESULT_FILES = [
    ("P1", "P1_evocast_result.json"),
    ("P2", "P2_is621_result.json"),
    ("P3", "P3_spas9_result.json"),
    ("P4", "P4_ennlovfz2_result.json"),
    ("P5", "P5_spufz1_v4_result.json"),
]

PASS_EMOJI   = "PASS"
FAIL_EMOJI   = "FAIL"
PEND_EMOJI   = "NOT_EVALUABLE"


def main() -> None:
    print("=" * 65)
    print("PEN-SCORE Pre-registered Prediction Summary")
    print("=" * 65)

    rows = []
    for pid, fname in RESULT_FILES:
        fpath = OUT_DIR / fname
        if not fpath.exists():
            print(f"  WARNING: {fname} not found - run script {int(pid[1:])+29}_*.py first")
            continue
        data = json.loads(fpath.read_text(encoding="utf-8"))
        rows.append({
            "ID":        pid,
            "Editor":    data.get("editor", ""),
            "Prediction": data.get("prediction", ""),
            "Verdict":   data.get("verdict", "MISSING"),
            "Reason":    data.get("reason", ""),
        })

    if not rows:
        print("No results found. Run P1-P5 scripts first.")
        return

    df = pd.DataFrame(rows)

    # Counts (cast to native int to avoid numpy int64 JSON serialisation errors)
    n_pass  = int((df["Verdict"] == "PASS").sum())
    n_fail  = int((df["Verdict"] == "FAIL").sum())
    n_pend  = int((df["Verdict"] == "NOT_EVALUABLE").sum())
    n_eval  = n_pass + n_fail

    print(f"\n{'ID':<5} {'Editor':<14} {'Verdict':<16} Reason")
    print("-" * 80)
    for _, r in df.iterrows():
        print(f"  {r['ID']:<4} {r['Editor']:<14} {r['Verdict']:<16} {r['Reason'][:55]}")

    print(f"\n{'='*50}")
    print(f"Summary: {n_pass}/{n_eval} evaluated predictions PASS")
    print(f"  PASS:           {n_pass}")
    print(f"  FAIL:           {n_fail}")
    print(f"  NOT_EVALUABLE:  {n_pend}  (pending sentinel resolution)")

    # Outcome lookup
    score_key = f"{n_pass}/5" if n_pend == 0 else f"{n_pass}/{n_eval} (+ {n_pend} pending)"
    policy = {
        "5/5": "Strong claim supported",
        "4/5": "Report which prediction failed",
        "3/5": "reframe as framework + lessons",
    }

    if n_fail == 0 and n_pend == 1:
        pub_policy = (f"4/4 evaluated pass + 1 pending (P4 sentinel). "
                      f"Strong claim supported; P4 pending resolution noted.")
    elif n_pass == 5:
        pub_policy = policy["5/5"]
    elif n_pass == 4:
        pub_policy = policy["4/5"]
    elif n_pass == 3:
        pub_policy = policy["3/5"]
    else:
        pub_policy = "Consult the pre-registration."

    print(f"\nOutcome: {pub_policy}")
    print(f"{'='*50}\n")

    # Write outputs
    summary = {
        "n_total":         len(RESULT_FILES),
        "n_evaluated":     n_eval,
        "n_pass":          n_pass,
        "n_fail":          n_fail,
        "n_not_evaluable": n_pend,
        "outcome": pub_policy,
        "predictions": rows,
    }

    # Convert any numpy int64 in rows to native Python int/str
    def _make_serializable(obj):
        if hasattr(obj, "item"):   # numpy scalar
            return obj.item()
        return obj

    def _clean(d):
        return {k: _make_serializable(v) for k, v in d.items()}

    summary["predictions"] = [_clean(r) for r in summary["predictions"]]

    out_json = OUT_DIR / "prediction_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    df.to_csv(OUT_DIR / "prediction_summary.csv", index=False)

    print(f"Written -> {out_json}")
    print(f"         -> {OUT_DIR / 'prediction_summary.csv'}")


if __name__ == "__main__":
    main()
