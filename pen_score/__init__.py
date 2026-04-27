"""pen-score: Multi-axis scoring framework for programmable genome editors.

Part of PEN-STACK - a unified computational infrastructure for non-destructive
genome engineering.

Public API (v0.1.3+)
---------------------
* :class:`~pen_score.api.Scorer` - multi-axis scoring engine.
* :func:`~pen_score.api.get_editor_metadata` - return PEN-COMPARE v3.2
  metadata (``intrinsic_cargo_mechanism``, ``cell_based_evidence``, alias
  resolution) for any editor by ID or deprecated alias.
* :class:`~pen_score.api.EditorMetadata` - frozen result dataclass.
"""

from __future__ import annotations

try:
    from pen_score._version import __version__
except ImportError:
    __version__ = "unknown"

from pen_score.api import EditorMetadata, Scorer, get_editor_metadata

__all__ = [
    "__version__",
    "EditorMetadata",
    "Scorer",
    "get_editor_metadata",
]
