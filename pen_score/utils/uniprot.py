"""UniProt REST API client.

Reused pattern from mech-class MECH-CLASS.  Key fixes inherited:
    - /stream with TSV + sequence returns HTTP 400 (payload too large).
      Use /search JSON + cursor pagination instead.
    - Retry on 429 with exponential back-off (UniProt rate limit: 50 req/s).
"""

from __future__ import annotations

import time

import requests

_BASE = "https://rest.uniprot.org/uniprotkb"
_HEADERS = {"Accept": "application/json"}
_RETRY_DELAYS = [1, 2, 4, 8]  # seconds


def _get_with_retry(url: str, params: dict | None = None) -> dict:
    for delay in _RETRY_DELAYS:
        resp = requests.get(url, params=params, headers=_HEADERS, timeout=30)
        if resp.status_code == 429:
            time.sleep(delay)
            continue
        resp.raise_for_status()
        return resp.json()
    resp.raise_for_status()
    return {}


def fetch_entry(accession: str) -> dict:
    """Fetch full UniProt entry for an accession."""
    return _get_with_retry(f"{_BASE}/{accession}")


def fetch_sequence(accession: str) -> str:
    """Return the canonical FASTA sequence for a UniProt accession."""
    entry = fetch_entry(accession)
    return entry["sequence"]["value"]


def fetch_sequence_length(accession: str) -> int:
    """Return the sequence length for a UniProt accession."""
    entry = fetch_entry(accession)
    return int(entry["sequence"]["length"])


def fetch_pfam_hits(accession: str) -> list[str]:
    """Return a list of Pfam accession strings for a UniProt entry."""
    entry = fetch_entry(accession)
    hits: list[str] = []
    for db_ref in entry.get("uniProtKBCrossReferences", []):
        if db_ref.get("database") == "Pfam":
            hits.append(db_ref["id"])
    return hits


def search_by_pfam(pfam_id: str, reviewed: bool = True, max_results: int = 500) -> list[dict]:
    """Return UniProt entries containing a given Pfam domain."""
    results: list[dict] = []
    review_filter = "reviewed:true" if reviewed else ""
    query_parts = [f"database:(type:pfam id:{pfam_id})"]
    if review_filter:
        query_parts.append(review_filter)
    query = " AND ".join(query_parts)

    url = f"{_BASE}/search"
    params: dict = {
        "query": query,
        "fields": "accession,protein_name,length,organism_name",
        "format": "json",
        "size": min(max_results, 500),
    }
    while True:
        data = _get_with_retry(url, params)
        results.extend(data.get("results", []))
        link = data.get("link")
        if not link or len(results) >= max_results:
            break
        url = link
        params = {}

    return results[:max_results]
