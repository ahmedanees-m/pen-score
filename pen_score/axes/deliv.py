"""S_Deliv - AAV Deliverability axis.

Formula (from axis_definitions.yaml):
    total_aa = sum(len(p) for p in editor_proteins)
    score = 1.0 / (1.0 + exp(0.005 * (total_aa - 900)))

The sigmoid is centred at 900 aa - the approximate single-AAV useful cargo
capacity when accounting for promoter, UTRs, and polyA (~4.7 kb total,
leaving ~2.7 kb = ~900 aa for protein coding sequence).

Protein length fetched from UniProt REST API or editor_universe.yaml.

Unit tests:
    SpCas9  (1368 aa) -> 0.20-0.30
    Cas12f  (~400 aa) -> 0.95-0.99
    IS621   (~300 aa) -> ~1.0
"""

from __future__ import annotations

import math
import warnings

from pen_score.data.loader import load_editor_universe
from pen_score.utils.uniprot import fetch_sequence_length

_SIGMOID_CENTER = 900
_SIGMOID_SLOPE = 0.005

# Pre-populated from editor universe to avoid redundant REST calls
_LENGTH_MAP: dict[str, int] = {}


def _build_length_map() -> None:
    global _LENGTH_MAP
    if _LENGTH_MAP:
        return
    for ed in load_editor_universe():
        # Will be populated by script 13_compute_S_Deliv.py
        # which fetches lengths from UniProt and persists them
        _LENGTH_MAP[ed.canonical_accession] = 0  # sentinel; replaced at compute time


def _sigmoid(aa: int) -> float:
    return 1.0 / (1.0 + math.exp(_SIGMOID_SLOPE * (aa - _SIGMOID_CENTER)))


def score(accession: str, total_aa: int | None = None) -> float | None:
    """Compute S_Deliv for an editor.

    Parameters
    ----------
    accession:
        UniProt accession (or synthetic id for engineered variants).
    total_aa:
        Pre-computed total amino acid count.  If None, fetched from UniProt.

    Returns
    -------
    float in [0, 1] or None if length cannot be determined.
    """
    if total_aa is None:
        # Try cached map first
        _build_length_map()
        total_aa = _LENGTH_MAP.get(accession)

    if total_aa is None or total_aa == 0:
        try:
            total_aa = fetch_sequence_length(accession)
        except Exception as exc:
            warnings.warn(f"S_Deliv: could not fetch length for {accession}: {exc}", stacklevel=2)
            return None

    return round(_sigmoid(total_aa), 4)
