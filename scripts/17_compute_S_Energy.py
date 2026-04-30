"""Compute S_Energy (Energy Independence) for all curated editors.

S_Energy is a binary axis (0.0 or 1.0) that scores whether an editor requires
cellular ATP hydrolysis for its core editing reaction.  The Walker A/B motif
scan is the primary method.  Multi-subunit ATPase systems use the
``walker_motif_override`` field in editor_universe.yaml.

Output
------
``data/pen-score/axes/energy/S_Energy.parquet``  (editor_id | S_Energy)

Usage
-----
    python scripts/17_compute_S_Energy.py
    python scripts/17_compute_S_Energy.py --output /path/to/S_Energy.parquet
    python scripts/17_compute_S_Energy.py --dry-run  # print scores, no write
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running as a script from the repo root.
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd

from pen_score.axes import energy as energy_axis
from pen_score.data.loader import load_editor_universe


def main(output: str | None = None, dry_run: bool = False) -> pd.DataFrame:
    editors = load_editor_universe()

    rows: list[dict] = []
    for ed in editors:
        override = getattr(ed, "walker_motif_override", None)
        s_energy = energy_axis.score(
            accession=ed.canonical_accession,
            walker_motif_override=override,
        )
        status = (
            f"{s_energy:.1f}"
            if s_energy is not None
            else "None (sentinel)"
        )
        print(f"  {ed.id:<20}  accession={ed.canonical_accession:<20}  S_Energy={status}")
        rows.append({"editor_id": ed.id, "S_Energy": s_energy})

    df = pd.DataFrame(rows)

    # Summary
    n_total = len(df)
    n_none = df["S_Energy"].isna().sum()
    n_indep = (df["S_Energy"] == 1.0).sum()
    n_dep = (df["S_Energy"] == 0.0).sum()
    print(f"\n--- Summary ---")
    print(f"  Total editors  : {n_total}")
    print(f"  S_Energy = 1.0 (energy-independent): {n_indep}")
    print(f"  S_Energy = 0.0 (ATP-dependent)     : {n_dep}")
    print(f"  S_Energy = None (sentinel)          : {n_none}")

    if not dry_run:
        out_path = Path(output) if output else Path("data/pen-score/axes/energy/S_Energy.parquet")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out_path, index=False)
        print(f"\n[OK] Written to {out_path}")

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", "-o", default=None, help="Output parquet path.")
    parser.add_argument("--dry-run", action="store_true", help="Print only; do not write.")
    args = parser.parse_args()
    main(output=args.output, dry_run=args.dry_run)
