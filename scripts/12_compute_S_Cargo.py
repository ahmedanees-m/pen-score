"""Compute S_Cargo (Cargo Capacity axis).

Formula: log10(cargo_capacity_bp) / log10(1e6), clipped to [0, 1].
Source:  editor_universe.yaml `cargo_capacity_bp` field (literature-curated).
No external compute - pure YAML lookup.

Run:
    docker run --rm \
        -v ~/pen-stack/data:/data \
        -v ~/pen-stack/code/repos/pen-score:/pkg \
        -w /pkg pen-stack/biophysics:0.1.0 \
        python scripts/12_compute_S_Cargo.py \
        2>&1 | tee ~/pen-stack/logs/pen-score/S_Cargo_$(date +%Y%m%d).log
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path("/data/pen-score/axes/cargo")
OUT.mkdir(parents=True, exist_ok=True)

LOG_MAX = 6.0  # log10(1_000_000 bp)


def cargo_score(cargo_bp: int | None) -> float | None:
    if not cargo_bp or cargo_bp <= 0:
        return None
    return round(max(0.0, min(1.0, math.log10(cargo_bp) / LOG_MAX)), 4)


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Cargo computation")
    print("=" * 60)

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]
        cargo_bp = ed.get("cargo_capacity_bp")
        sc = cargo_score(cargo_bp)
        rows.append({
            "editor_id": eid,
            "canonical_accession": acc,
            "cargo_capacity_bp": cargo_bp,
            "S_Cargo": sc,
        })
        print(f"  {eid:<22} {str(cargo_bp or 'N/A'):>10} bp  S_Cargo={sc}")

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "cargo_scores.parquet", index=False)
    df.to_csv(OUT / "cargo_scores.csv", index=False)
    print(f"\nWritten -> {OUT}/cargo_scores.parquet (.csv)")
    print(f"  Non-null scores: {df['S_Cargo'].notna().sum()}/{len(df)}")


if __name__ == "__main__":
    main()
