"""Compute S_Prog (Programmability axis).

Formula: Binary: 1.0 if RNA-guided, 0.0 if site-specific att-site recombinase.
Source:  editor_universe.yaml `rna_guided` field (Boolean); cross-checked with
         MECH-CLASS Tier B programmability sub-class.

No external compute - pure YAML lookup.

Run:
    docker run --rm \
        -v ~/pen-stack/data:/data \
        -v ~/pen-stack/code/repos/pen-score:/pkg \
        -w /pkg pen-stack/biophysics:0.1.0 \
        python scripts/15_compute_S_Prog.py \
        2>&1 | tee ~/pen-stack/logs/pen-score/S_Prog_$(date +%Y%m%d).log
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path("/data/pen-score/axes/prog")
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Prog computation")
    print("=" * 60)

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]
        rna = ed.get("rna_guided", False)
        sp = 1.0 if rna else 0.0
        rows.append({
            "editor_id": eid,
            "canonical_accession": acc,
            "rna_guided": rna,
            "S_Prog": sp,
        })
        print(f"  {eid:<22} rna_guided={str(rna):<5}  S_Prog={sp}")

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "prog_scores.parquet", index=False)
    df.to_csv(OUT / "prog_scores.csv", index=False)
    print(f"\nWritten -> {OUT}/prog_scores.parquet (.csv)")
    rna_count = sum(1 for r in rows if r["rna_guided"])
    print(f"  RNA-guided (S_Prog=1.0): {rna_count}  |  Site-specific (S_Prog=0.0): {len(rows)-rna_count}")


if __name__ == "__main__":
    main()
