"""Compute S_Mature (Therapeutic Maturity axis).

Formula:
    raw_count = PubMed hit count for (editor_terms) AND (clinical OR therapeutic ...)
    score = log10(raw_count + 1) / log10(max_count + 1)

Normalised over the full editor universe so SpCas9 ~ 1.0.
Two-pass: first fetch all counts, then normalize.

Run:
    docker run --rm \
        -v ~/pen-stack/data:/data \
        -v ~/pen-stack/code/repos/pen-score:/pkg \
        -w /pkg pen-stack/biophysics:0.1.0 \
        python scripts/16_compute_S_Mature.py \
        2>&1 | tee ~/pen-stack/logs/pen-score/S_Mature_$(date +%Y%m%d).log
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import pandas as pd
import requests
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path("/data/pen-score/axes/mature")
OUT.mkdir(parents=True, exist_ok=True)

NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINICAL_QUERY = "clinical[tw] OR preclinical[tw] OR therapeutic[tw] OR gene_therapy[tw]"
REQ_DELAY_S = 0.34


def pubmed_count(terms: list[str]) -> int:
    q = "(" + " OR ".join(terms) + ") AND (" + CLINICAL_QUERY + ")"
    params = {"db": "pubmed", "term": q, "retmax": "0", "retmode": "json"}
    r = requests.get(NCBI_ESEARCH, params=params, timeout=15)
    r.raise_for_status()
    return int(r.json()["esearchresult"]["count"])


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Mature computation")
    print("=" * 60)

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    # Pass 1 - fetch PubMed counts
    print("Pass 1: fetching PubMed citation counts ...")
    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]
        terms = ed.get("references_used_for_pubmed") or [eid]
        try:
            count = pubmed_count(terms)
            rows.append({"editor_id": eid, "canonical_accession": acc,
                         "pubmed_count": count, "S_Mature": None})
            print(f"  {eid:<22} {count} hits")
        except Exception as exc:
            rows.append({"editor_id": eid, "canonical_accession": acc,
                         "pubmed_count": 0, "S_Mature": None})
            print(f"  {eid:<22} ERROR: {exc}")
        sys.stdout.flush()
        time.sleep(REQ_DELAY_S)

    # Pass 2 - normalize
    max_count = max(r["pubmed_count"] for r in rows) if rows else 1
    print(f"\nmax_count = {max_count}  (normalisation denominator)")
    print("\nPass 2: normalizing ...")
    for row in rows:
        c = row["pubmed_count"]
        row["S_Mature"] = round(math.log10(c + 1) / math.log10(max_count + 1), 4) if max_count > 0 else 0.0
        print(f"  {row['editor_id']:<22} count={c}  S_Mature={row['S_Mature']}")

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "mature_scores.parquet", index=False)
    df.to_csv(OUT / "mature_scores.csv", index=False)
    print(f"\nWritten -> {OUT}/mature_scores.parquet (.csv)")
    print(f"  Computed: {df['S_Mature'].notna().sum()}/{len(df)}")


if __name__ == "__main__":
    main()
