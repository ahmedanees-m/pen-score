"""S_Prog - Programmability axis.

Formula (from axis_definitions.yaml):
    Binary flag from MECH-CLASS Tier B:
        1.0  if editor is RNA-guided (rna_guided=True in editor_universe.yaml)
        0.5  if editor has retargeting potential but not RNA-guided (e.g., TALE)
        0.0  if editor is fixed to a specific att-site or target sequence

Source: editor_universe.yaml `rna_guided` field (primary) and MECH-CLASS Tier B
`programmable` sub-class flag (secondary, for proteins not in the universe).

This axis is intentionally NOT orthogonal to mechanism class - it consumes
MECH-CLASS output directly (same as S_DSB).  Both axes are derived from the
same upstream model; orthogonality between them is not required.

Unit tests:
    SpCas9  -> 1.0  (RNA-guided)
    IS621   -> 1.0  (RNA-guided bridge recombinase)
    Cre     -> 0.0  (fixed att-site, not programmable)
    Bxb1    -> 0.0  (fixed attP/attB, not programmable)
"""

from __future__ import annotations

from pen_score.data.loader import load_editor_universe

# Pre-built lookup: canonical_accession -> S_Prog
_PROG_MAP: dict[str, float] = {}


def _build_prog_map() -> None:
    global _PROG_MAP
    if _PROG_MAP:
        return
    for ed in load_editor_universe():
        _PROG_MAP[ed.canonical_accession] = 1.0 if ed.rna_guided else 0.0
        _PROG_MAP[ed.id] = 1.0 if ed.rna_guided else 0.0


def score(accession: str) -> float | None:
    """Compute S_Prog for an editor.

    Returns
    -------
    1.0 if RNA-guided, 0.0 if site-specific, None if not in the universe.
    """
    _build_prog_map()
    return _PROG_MAP.get(accession)
