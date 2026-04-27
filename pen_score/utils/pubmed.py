"""NCBI E-utilities client for S_Mature PubMed citation counts."""

from __future__ import annotations

import time

import requests

_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_TOOL = "pen-score"
_EMAIL = "ahmedaneesm@gmail.com"


def fetch_pubmed_count(search_terms: list[str], extra_terms: str = "gene therapy") -> int:
    """Return PubMed hit count for an editor + therapeutic context query.

    Parameters
    ----------
    search_terms:
        Editor name synonyms (e.g. ['SpCas9', 'Cas9', 'CRISPR-Cas9']).
    extra_terms:
        Additional context (default: 'gene therapy').

    Returns
    -------
    Integer hit count.
    """
    term_str = " OR ".join(f'"{t}"' for t in search_terms)
    query = f"({term_str}) AND ({extra_terms})"
    url = f"{_BASE}/esearch.fcgi"
    params: dict[str, str | int] = {
        "db": "pubmed",
        "term": query,
        "retmax": 0,
        "retmode": "json",
        "tool": _TOOL,
        "email": _EMAIL,
    }
    for attempt in range(4):
        resp = requests.get(url, params=params, timeout=20)
        if resp.status_code == 429:
            time.sleep(2**attempt)
            continue
        resp.raise_for_status()
        return int(resp.json()["esearchresult"]["count"])
    resp.raise_for_status()
    return 0


def batch_fetch_counts(
    editor_term_map: dict[str, list[str]],
    extra_terms: str = "gene therapy",
    sleep_between: float = 0.34,  # NCBI limit: 3 req/s without API key
) -> dict[str, int]:
    """Fetch PubMed counts for multiple editors.

    Parameters
    ----------
    editor_term_map:
        Dict mapping editor_id -> list of search terms.
    extra_terms:
        Extra context terms for all queries.
    sleep_between:
        Seconds to sleep between requests.

    Returns
    -------
    Dict mapping editor_id -> count.
    """
    counts: dict[str, int] = {}
    for editor_id, terms in editor_term_map.items():
        counts[editor_id] = fetch_pubmed_count(terms, extra_terms=extra_terms)
        time.sleep(sleep_between)
    return counts
