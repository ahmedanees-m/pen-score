"""Compute S_DSB (DSB Avoidance axis).

Formula:
    base = 1.0 - P(DSB_NUCLEASE)   # from MECH-CLASS Tier A
    if composite_flag: score = min(1.0, base + 0.1)
    else:              score = base

    The composite_flag is True only for IS110-class editors with BOTH
    PF01548 (IS110 DEDD RuvC-fold N-terminal) AND PF02371 (IS110 Tnp C-terminal)
    Pfam domains co-occurring. This gate was introduced in mech-class v0.5.1 to
    prevent false-positive composite flags on dual-domain non-IS110 proteins
    (e.g., SpCas9 RuvC+HNH). See SCORE_PROVENANCE.md the Upstream Dependencies section.

Primary method: MECH-CLASS >=v0.5.1 Tier A classifier.
  Requires mech-class installed and ESM-2 embeddings (slow, CPU-only).

Fallback (if mech_class unavailable): bucket heuristic from mechanism_bucket
  field in editor_universe.yaml - encodes the MECH-CLASS v0.5.1-aligned output.
  DSB_NUCLEASE -> 0.0 | DSB_FREE_TRANSEST_RECOMBINASE -> 0.9 | TRANSPOSASE -> 0.5

  Note: DSB_NUCLEASE bucket maps to 0.0 (not 0.1). The composite bonus (+0.1)
  requires the IS110-specific PF01548 and PF02371 Pfam gate which no DSB_NUCLEASE
  editor in our universe satisfies. Aligns bucket heuristic with v0.5.1 model.

Run with mech-class (recommended):
    docker run --rm \
        -v ~/pen-stack/data:/data \
        -v ~/pen-stack/code/repos/pen-score:/pkg \
        -v ~/pen-stack/code/repos/mech-class:/mech_class_repo \
        -w /pkg pen-stack/biophysics:0.1.0 \
        bash -c "pip install lightgbm --quiet && PYTHONPATH=/mech_class_repo \
                 python scripts/10_compute_S_DSB.py" \
        2>&1 | tee ~/pen-stack/logs/pen-score/S_DSB_$(date +%Y%m%d).log
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path("/data/pen-score/axes/dsb")
OUT.mkdir(parents=True, exist_ok=True)

# Try to import mech_class from local repo mount
sys.path.insert(0, "/mech_class_repo")


def _bucket_heuristic(mechanism_bucket: str) -> float:
    """Fallback S_DSB from mechanism_bucket field.

    Aligned with mech-class v0.5.1: DSB_NUCLEASE -> 0.0.
    The composite bonus (+0.1) requires PF01548 and PF02371 (IS110 Pfam gate),
    which no DSB_NUCLEASE editor in the current universe satisfies.
    Previously 0.1 in v0.5.0 (coincided with FP composite bonus on SpCas9).
    """
    if "DSB_NUCLEASE" in mechanism_bucket:
        return 0.0   # v0.5.1: no composite bonus for non-IS110 nucleases
    if "DSB_FREE" in mechanism_bucket or "RECOMBINASE" in mechanism_bucket:
        return 0.9
    if "TRANSPOSASE" in mechanism_bucket:
        return 0.5
    return 0.5


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_DSB computation")
    print("=" * 60)

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    # Attempt to load mech-class model
    predictor = None
    use_model = False
    try:
        from mech_class.api import Predictor  # type: ignore[import]
        import requests as _req  # noqa: F401 (sanity check)

        predictor = Predictor.load(model_dir="/data/models")
        print("mech_class Predictor loaded - using Tier A model (>=v0.5.1 required)\n")
        use_model = True
    except ImportError as exc:
        print(f"mech_class not available ({exc}) - using bucket heuristic (v0.5.1-aligned)\n")
    except Exception as exc:
        print(f"mech_class load failed ({exc}) - using bucket heuristic (v0.5.1-aligned)\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]
        bucket = ed.get("mechanism_bucket", "")
        sentinel = acc.startswith("REQUIRES") or acc == "NO_UNIPROT"

        sc: float | None = None
        note = ""

        if use_model and not sentinel:
            # Strip variant suffix for UniProt lookup
            lookup = acc.split("_")[0] if "_" in acc and not acc.startswith("A0A") else acc
            try:
                # fetch_sequence from UniProt for ESM-2 embedding
                import requests as req_mod
                r = req_mod.get(f"https://rest.uniprot.org/uniprotkb/{lookup}.json", timeout=15)
                if r.status_code != 200:
                    raise ValueError(f"UniProt HTTP {r.status_code}")
                seq = r.json()["sequence"]["value"]
                pfam_ids = [x["id"] for x in r.json().get("uniProtKBCrossReferences", [])
                            if x.get("database") == "Pfam"]

                pred = predictor.predict_from_sequence(
                    accession=lookup, sequence=seq, pfam_hits=pfam_ids
                )
                p_dsb = pred.tier_a_probabilities.get("DSB_NUCLEASE", 0.0)  # type: ignore[attr-defined]
                base = 1.0 - p_dsb
                sc = round(min(1.0, base + 0.1) if pred.composite else float(base), 4)
                note = f"mech_class p_dsb={p_dsb:.3f}"
            except Exception as exc:
                note = f"model_err: {str(exc)[:60]}"

        if sc is None:
            sc = _bucket_heuristic(bucket)
            note = note or "bucket_heuristic"

        rows.append({
            "editor_id": eid, "canonical_accession": acc,
            "mechanism_bucket": bucket, "S_DSB": sc, "note": note,
        })
        print(f"  {eid:<22} S_DSB={sc:.4f}  [{note[:55]}]")

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "dsb_scores.parquet", index=False)
    df.to_csv(OUT / "dsb_scores.csv", index=False)
    model_n = sum(1 for r in rows if "mech_class" in r["note"])
    heur_n = len(rows) - model_n
    print(f"\nWritten -> {OUT}/dsb_scores.parquet (.csv)")
    print(f"  Tier A model: {model_n}  |  Bucket heuristic: {heur_n}")


if __name__ == "__main__":
    main()
