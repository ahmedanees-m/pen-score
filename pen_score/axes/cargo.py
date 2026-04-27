"""S_Cargo - Cargo Capacity axis.

Formula (from axis_definitions.yaml):
    log_cargo = log10(cargo_capacity_bp)
    score = (log_cargo - 0) / (log10(1e6) - 0)   # normalise over 1 bp .. 1 Mb
    score = clip(score, 0.0, 1.0)

Source of cargo_capacity_bp: editor_universe.yaml (literature-curated).
No external compute required; pure lookup.

Unit tests:
    PE2  (200 bp)  -> 0.38-0.42
    SpCas9+HDR (1000 bp) -> ~0.50
    IS621 (1 Mb)   -> 1.0
"""

from __future__ import annotations

import math

from pen_score.data.loader import load_editor_universe

_LOG_MIN = 0.0  # log10(1 bp)
_LOG_MAX = 6.0  # log10(1 000 000 bp)

# Build lookup once at module import
_CARGO_MAP: dict[str, int] = {}


def _build_cargo_map() -> None:
    global _CARGO_MAP
    if _CARGO_MAP:
        return
    for ed in load_editor_universe():
        _CARGO_MAP[ed.canonical_accession] = ed.cargo_capacity_bp
        # Also index by editor id for convenience
        _CARGO_MAP[ed.id] = ed.cargo_capacity_bp


def score(accession: str) -> float | None:
    """Compute S_Cargo for a given UniProt accession or editor id.

    Returns
    -------
    float in [0, 1] or None if the editor is not found in the universe.
    """
    _build_cargo_map()
    cargo_bp = _CARGO_MAP.get(accession)
    if cargo_bp is None:
        return None
    if cargo_bp <= 0:
        return 0.0
    log_cargo = math.log10(cargo_bp)
    raw = (log_cargo - _LOG_MIN) / (_LOG_MAX - _LOG_MIN)
    return round(max(0.0, min(1.0, raw)), 4)
