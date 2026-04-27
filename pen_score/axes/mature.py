"""S_Mature - Therapeutic Maturity axis.

Formula (from axis_definitions.yaml):
    raw_count = PubMed citation count for clinical/preclinical terms
    score = log10(raw_count + 1) / log10(max_count_over_universe + 1)

where max_count_over_universe is computed across all editors (normalises
well-known editors like SpCas9 to 1.0 and novel editors to ~0.0).

PubMed queries via NCBI E-utilities:
    esearch.fcgi?db=pubmed&term=<editor_name>+AND+(clinical+OR+preclinical+
        OR+therapeutic+OR+gene+therapy)&retmax=0&usehistory=y

Run script 16_compute_S_Mature.py to batch-fetch counts and cache them.
Live computation is available but slower (1 API call per editor, ~1 s each).

Requires: requests (core dep, always available).

Unit tests:
    SpCas9  -> ~1.0  (thousands of citations)
    IS621   -> ~0.1-0.3  (recent, few clinical citations)
    PE2     -> ~0.6-0.8  (growing clinical literature)
"""

from __future__ import annotations

import math
import warnings
from pathlib import Path

_NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_MAX_COUNT_FALLBACK = 10_000  # approximate upper bound for normalisation before full run


def _query_pubmed_count(search_terms: list[str]) -> int:
    """Return PubMed hit count for (editor_terms) AND (clinical terms)."""
    import requests

    term_str = " OR ".join(f'"{t}"' for t in search_terms)
    clinical_terms = 'clinical[tw] OR preclinical[tw] OR therapeutic[tw] OR "gene therapy"[tw]'
    query = f"({term_str}) AND ({clinical_terms})"
    url = f"{_NCBI_EUTILS_BASE}/esearch.fcgi"
    params: dict[str, str | int] = {"db": "pubmed", "term": query, "retmax": 0, "retmode": "json"}
    resp = requests.get(url, params=params, timeout=15)
    resp.raise_for_status()
    return int(resp.json()["esearchresult"]["count"])


def score(
    accession: str,
    search_terms: list[str] | None = None,
    max_count: int = _MAX_COUNT_FALLBACK,
    cached_parquet: Path | None = None,
) -> float | None:
    """Compute S_Mature for an editor.

    Parameters
    ----------
    accession:
        UniProt accession or editor id.
    search_terms:
        PubMed search terms for this editor (e.g. ['SpCas9', 'Cas9', 'CRISPR-Cas9']).
        Taken from editor_universe.yaml `references_used_for_pubmed` field if None.
    max_count:
        Universe-wide maximum citation count for normalisation.
        Set correctly by script 16 after running all editors.
    cached_parquet:
        Path to pre-computed ``mature_scores.parquet`` from script 16.

    Returns
    -------
    float in [0, 1] or None if not computable.
    """
    if cached_parquet and Path(cached_parquet).exists():
        import pandas as pd

        df = pd.read_parquet(cached_parquet)
        row = df[df["canonical_accession"] == accession]
        if not row.empty:
            return float(row["S_Mature"].iloc[0])

    if search_terms is None:
        # Try to look up from universe yaml
        from pen_score.data.loader import load_editor_universe

        for ed in load_editor_universe():
            if ed.canonical_accession == accession or ed.id == accession:
                search_terms = ed.references_used_for_pubmed or [ed.id]
                break

    if not search_terms:
        warnings.warn(
            f"S_Mature: no search terms for {accession}. Add 'references_used_for_pubmed' "
            "in editor_universe.yaml.",
            stacklevel=2,
        )
        return None

    try:
        count = _query_pubmed_count(search_terms)
        return round(math.log10(count + 1) / math.log10(max_count + 1), 4)
    except Exception as exc:
        warnings.warn(f"S_Mature PubMed query failed for {accession}: {exc}", stacklevel=2)
        return None
