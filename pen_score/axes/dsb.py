"""S_DSB - DSB Avoidance axis.

Formula (from axis_definitions.yaml):
    base = 1.0 - P(DSB_NUCLEASE)          # from MECH-CLASS Tier A
    if composite_flag (IS110 detected):
        score = min(1.0, base + 0.1)
    else:
        score = base

IS110 Tier-A hard gate (mech-class v0.5.2+):
    if tier_a_gate_override is True:
        score = 1.0                         # domain-evidence confidence 0.99
    This gate fires when PF01548 and PF02371 are both present, ensuring
    novel IS110 proteins (e.g. IS621) are not mis-scored by OOD ML features.

Requires the [mech-class] optional extra.  Returns None if mech-class is
not installed or prediction fails.
"""

from __future__ import annotations

import warnings


def score(accession: str, pfam_hits: list[str] | None = None) -> float | None:
    """Compute S_DSB for a given UniProt accession.

    Parameters
    ----------
    accession:
        UniProt accession (e.g. 'Q99ZW2' for SpCas9).
    pfam_hits:
        Optional pre-fetched Pfam accession list.  If None, fetched via
        UniProt REST API inside mech_class.

    Returns
    -------
    float in [0, 1] or None if mech-class is unavailable.
    """
    try:
        from mech_class.api import Predictor
    except ImportError:
        warnings.warn(
            "mech-class is not installed. S_DSB returns None. "
            "Install with: pip install pen-score[mech-class]",
            stacklevel=2,
        )
        return None

    try:
        predictor = Predictor.load()
        pred = predictor.predict_from_sequence(accession=accession, pfam_hits=pfam_hits or [])

        # IS110 Tier-A hard gate (mech-class v0.5.2+).
        # When PF01548 and PF02371 are both present the gate overrides the ML
        # probability to prevent OOD mis-scoring of novel IS110 proteins.
        # IS110-family bridge recombinases have domain-evidence DSB_FREE
        # confidence >= 0.90; with the IS110 composite bonus (+0.1) this
        # always yields S_DSB = 1.0.
        if getattr(pred, "tier_a_gate_override", False):
            return 1.0

        p_dsb_nuclease = pred.tier_a_probabilities.get("DSB_NUCLEASE", 0.0)
        base = 1.0 - p_dsb_nuclease
        if pred.composite:
            return min(1.0, base + 0.1)
        return float(base)
    except Exception as exc:
        warnings.warn(f"S_DSB computation failed for {accession}: {exc}", stacklevel=2)
        return None
