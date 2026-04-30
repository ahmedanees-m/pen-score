"""Compute S_Deliv (AAV Deliverability axis).

Formula: sigmoid(0.005 * (900 - total_aa))
Source:  Protein length from UniProt REST API per accession.
         NO_UNIPROT entries use protein_length_aa from editor_universe.yaml.
         REQUIRES_STEP7 sentinels are skipped (S_Deliv = NaN).

Run:
    docker run --rm \
        -v ~/pen-stack/data:/data \
        -v ~/pen-stack/code/repos/pen-score:/pkg \
        -w /pkg pen-stack/biophysics:0.1.0 \
        python scripts/13_compute_S_Deliv.py \
        2>&1 | tee ~/pen-stack/logs/pen-score/S_Deliv_$(date +%Y%m%d).log
"""

from __future__ import annotations

import math
import time
from pathlib import Path

import pandas as pd
import requests
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path("/data/pen-score/axes/deliv")
OUT.mkdir(parents=True, exist_ok=True)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
REQ_DELAY_S = 0.25


def sigmoid(aa: int) -> float:
    return round(1.0 / (1.0 + math.exp(0.005 * (aa - 900))), 4)


def fetch_length(accession: str) -> int:
    # Strip engineered-variant suffix (e.g. A0A0L0H5U9_V4 -> A0A0L0H5U9)
    base = accession.split("_")[0] if "_" in accession and not accession.startswith("A0A") else accession
    r = requests.get(f"{UNIPROT_BASE}/{base}.json", timeout=15)
    if r.status_code != 200:
        raise ValueError(f"UniProt HTTP {r.status_code} for {base}")
    return int(r.json()["sequence"]["length"])


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Deliv computation")
    print("=" * 60)

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]

        # Sentinel - sequence not yet resolved
        if acc.startswith("REQUIRES"):
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "protein_length_aa": None, "S_Deliv": None, "note": "REQUIRES_STEP7",
            })
            print(f"  {eid:<22} SKIP (sentinel)")
            continue

        # NO_UNIPROT - use hardcoded length
        if acc == "NO_UNIPROT":
            length = ed.get("protein_length_aa")
            sc = sigmoid(length) if length else None
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "protein_length_aa": length, "S_Deliv": sc, "note": "NO_UNIPROT",
            })
            print(f"  {eid:<22} {length} aa  S_Deliv={sc}  [NO_UNIPROT]")
            continue

        # Normal UniProt accession
        try:
            length = fetch_length(acc)
            sc = sigmoid(length)
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "protein_length_aa": length, "S_Deliv": sc, "note": "",
            })
            print(f"  {eid:<22} {length} aa  S_Deliv={sc}")
        except Exception as exc:
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "protein_length_aa": None, "S_Deliv": None, "note": str(exc),
            })
            print(f"  {eid:<22} ERROR: {exc}")
        time.sleep(REQ_DELAY_S)

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "deliv_scores.parquet", index=False)
    df.to_csv(OUT / "deliv_scores.csv", index=False)
    n_ok = df["S_Deliv"].notna().sum()
    n_skip = df["note"].str.startswith("REQUIRES").sum()
    print(f"\nWritten -> {OUT}/deliv_scores.parquet (.csv)")
    print(f"  Computed: {n_ok}  |  Skipped (sentinel): {n_skip}  |  Errors: {df['S_Deliv'].isna().sum()-n_skip}")


if __name__ == "__main__":
    main()
