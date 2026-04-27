"""Public API: Scorer class for multi-axis editor scoring.

v0.1.3 additions
----------------
* :func:`get_editor_metadata` - return PEN-COMPARE v3.2-compatible metadata for a
  single editor by ID or alias.  Resolves deprecated aliases (e.g. "IS622" ->
  "ISCro4") with a :class:`DeprecationWarning`.
* :class:`EditorMetadata` - frozen dataclass for the returned record.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from pydantic import BaseModel, Field

from pen_score.data.loader import load_editor_universe, load_use_case_profiles

# Valid axis names for exclude_axes validation
_VALID_AXES: frozenset[str] = frozenset(
    ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature", "S_Energy"]
)


# ---------------------------------------------------------------------------
# EditorMetadata - v0.1.3 / PEN-COMPARE v3.2 metadata API
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EditorMetadata:
    """Metadata record for a single curated editor (pen-score v0.1.3+).

    All fields are sourced from ``editor_universe.yaml`` schema v1.0.7.

    Attributes
    ----------
    editor_id:
        Canonical editor identifier (e.g. "ISCro4").
    canonical_name:
        Same as *editor_id*.  Exists for symmetry with alias resolution.
    aliases:
        Deprecated alternative names that resolve to this editor
        (e.g. ``["IS622"]`` for ISCro4).
    uniprot:
        UniProt canonical accession, or a sentinel string such as
        "REQUIRES_STEP7" / "NO_UNIPROT" if unresolved.
    organism:
        Source organism name string.
    intrinsic_cargo_mechanism:
        **True** if the editor carries/inserts cargo as part of its catalytic
        mechanism (IS110 bridge recombinases, CAST transposases, site-specific
        recombinases).  **False** if cargo requires an external HDR donor
        template (SpCas9+HDR, prime editors).
        Required by PEN-COMPARE v3.2 Gate 3.
    cell_based_evidence:
        **True** if peer-reviewed mammalian cell activity data exists at >1 %
        editing efficiency.  **False** for in-vitro / E. coli only.
        Required by PEN-COMPARE v3.2 TRUE_WRITER tier.
    cell_based_sources:
        List of supporting citations (DOI strings or author-year notes).
    """

    editor_id: str
    canonical_name: str
    aliases: list[str]
    uniprot: str | None
    organism: str | None
    intrinsic_cargo_mechanism: bool
    cell_based_evidence: bool
    cell_based_sources: list[str]


def get_editor_metadata(editor_id: str) -> EditorMetadata:
    """Return v3.2-compatible metadata for a single editor.

    Looks up *editor_id* in the pen-score editor universe.  If the exact ID
    is not found, the function searches the ``aliases`` field of every editor
    and resolves the match with a :class:`DeprecationWarning` pointing to the
    canonical name.

    Parameters
    ----------
    editor_id:
        Canonical editor ID (e.g. ``"ISCro4"``) or a deprecated alias
        (e.g. ``"IS622"``).

    Returns
    -------
    EditorMetadata
        Frozen metadata record.

    Raises
    ------
    KeyError
        If *editor_id* is not found as an ID or alias in
        ``editor_universe.yaml`` v1.0.7.

    Examples
    --------
    >>> from pen_score import get_editor_metadata
    >>> md = get_editor_metadata("ISCro4")
    >>> md.intrinsic_cargo_mechanism
    True
    >>> md.cell_based_evidence
    True
    >>> md = get_editor_metadata("IS621")
    >>> md.cell_based_evidence   # E. coli + cryo-EM only; no robust human-cell data
    False
    """
    editors = load_editor_universe()

    # Direct lookup by canonical ID
    for ed in editors:
        if ed.id == editor_id:
            return _build_metadata(ed)

    # Alias lookup with deprecation warning
    for ed in editors:
        if editor_id in (ed.aliases or []):
            warnings.warn(
                f"'{editor_id}' is a deprecated alias for '{ed.id}'. "
                f"Use '{ed.id}' (canonical name) in new code.",
                DeprecationWarning,
                stacklevel=2,
            )
            return _build_metadata(ed)

    raise KeyError(
        f"Editor '{editor_id}' not found in pen-score editor_universe "
        f"v1.0.7 (neither as a canonical ID nor as an alias)."
    )


def _build_metadata(ed: Any) -> EditorMetadata:
    """Build an :class:`EditorMetadata` from an :class:`~pen_score.data.loader.EditorEntry`."""
    return EditorMetadata(
        editor_id=ed.id,
        canonical_name=ed.id,
        aliases=list(ed.aliases or []),
        uniprot=ed.canonical_accession,
        organism=ed.organism,
        intrinsic_cargo_mechanism=bool(ed.intrinsic_cargo_mechanism),
        cell_based_evidence=bool(ed.cell_based_evidence),
        cell_based_sources=list(ed.cell_based_sources or []),
    )


class AxisScores(BaseModel):
    S_DSB: float | None = Field(None, ge=0.0, le=1.0)
    S_Spec: float | None = Field(None, ge=0.0, le=1.0)
    S_Cargo: float | None = Field(None, ge=0.0, le=1.0)
    S_Deliv: float | None = Field(None, ge=0.0, le=1.0)
    S_Immuno: float | None = Field(None, ge=0.0, le=1.0)
    S_Prog: float | None = Field(None, ge=0.0, le=1.0)
    S_Mature: float | None = Field(None, ge=0.0, le=1.0)
    S_Energy: float | None = Field(None, ge=0.0, le=1.0)


class ScoringResult(BaseModel):
    accession: str
    editor_id: str | None = None
    use_case: str
    axes: AxisScores
    pen_score: float | None = Field(None, ge=0.0, le=1.0)
    reasoning: list[str] = Field(default_factory=list)
    axes_missing: list[str] = Field(default_factory=list)


class Scorer:
    """Multi-axis scorer for programmable genome editors.

    Usage::

        scorer = Scorer.load()
        result = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
        )
        print(result.pen_score)   # 0.9290
        print(result.reasoning)  # axis-level strength/weakness bullets

        top5 = scorer.select_editor(
            use_case="human_therapeutic_aav_insertion",
            top_k=5,
            require_dsb_free=True,
        )
    """

    def __init__(self, scorecard_path: Path | None = None) -> None:
        self._scorecard_path = scorecard_path
        self._scorecard: pd.DataFrame | None = None
        self._use_case_profiles = load_use_case_profiles()

    @classmethod
    def load(cls, scorecard_path: Path | None = None) -> Scorer:
        """Load a Scorer. If scorecard_path is given, load pre-computed scores."""
        return cls(scorecard_path=scorecard_path)

    # ------------------------------------------------------------------
    # Core scoring
    # ------------------------------------------------------------------

    def score_editor(
        self,
        accession: str,
        use_case: str = "human_therapeutic_aav_insertion",
        exclude_axes: list[str] | None = None,
    ) -> ScoringResult:
        """Score a single editor by UniProt accession.

        Parameters
        ----------
        accession:
            UniProt canonical accession.
        use_case:
            One of the pre-registered use-case profile keys.
        exclude_axes:
            Axis names to exclude from the composite PenScore.  Weights of
            remaining axes are renormalised to sum to 1.0.  Useful for
            biophysical-only comparisons (e.g. ``exclude_axes=['S_Mature']``
            gives a fair score for newly-characterised editors that lack
            clinical literature regardless of their biophysical performance).
            Valid values: S_DSB, S_Spec, S_Cargo, S_Deliv, S_Immuno, S_Prog,
            S_Mature, S_Energy.

        Raises
        ------
        ValueError
            If any name in *exclude_axes* is not a recognised axis.
        """
        if exclude_axes:
            invalid = [ax for ax in exclude_axes if ax not in _VALID_AXES]
            if invalid:
                raise ValueError(
                    f"Unknown axis name(s) in exclude_axes: {invalid}. "
                    f"Valid axes: {sorted(_VALID_AXES)}"
                )

        if (
            self._scorecard is not None
            and accession in self._scorecard["canonical_accession"].values
        ):
            row = self._scorecard[self._scorecard["canonical_accession"] == accession].iloc[0]
            axes = AxisScores(
                S_DSB=row.get("S_DSB"),
                S_Spec=row.get("S_Spec"),
                S_Cargo=row.get("S_Cargo"),
                S_Deliv=row.get("S_Deliv"),
                S_Immuno=row.get("S_Immuno"),
                S_Prog=row.get("S_Prog"),
                S_Mature=row.get("S_Mature"),
                S_Energy=row.get("S_Energy"),
            )
        else:
            axes = self._compute_axes_live(accession)

        weights = self._use_case_profiles.get(use_case, self._default_weights())
        pen_score, missing = self._composite(axes, weights, exclude_axes=exclude_axes)
        reasoning = self._generate_reasoning(axes, weights)

        return ScoringResult(
            accession=accession,
            use_case=use_case,
            axes=axes,
            pen_score=pen_score,
            reasoning=reasoning,
            axes_missing=missing,
        )

    def get_scorecard(self, use_case: str = "human_therapeutic_aav_insertion") -> pd.DataFrame:
        """Return scores for the full editor universe as a DataFrame."""
        if self._scorecard_path and self._scorecard_path.exists():
            return pd.read_parquet(self._scorecard_path)
        editors = load_editor_universe()
        rows: list[dict[str, Any]] = []
        for ed in editors:
            result = self.score_editor(ed.canonical_accession, use_case=use_case)
            row: dict[str, Any] = {
                "editor_id": ed.id,
                "canonical_accession": ed.canonical_accession,
            }
            row.update(result.axes.model_dump())
            row["PenScore"] = result.pen_score
            rows.append(row)
        return pd.DataFrame(rows)

    def select_editor(
        self,
        use_case: str,
        top_k: int = 5,
        require_dsb_free: bool = False,
        filters: dict[str, float] | None = None,
    ) -> pd.DataFrame:
        """Return top-k ranked editors for a use case.

        Parameters
        ----------
        use_case:
            One of the five pre-registered use-case profile keys.
        top_k:
            Number of editors to return (default 5).
        require_dsb_free:
            If True, restrict candidates to editors with S_DSB >= 0.85
            (IS110 integrases, site-specific recombinases, compact transposases).
        filters:
            Additional column-level minimum-value filters applied before ranking.
        """
        sc = self.get_scorecard(use_case=use_case)
        if require_dsb_free:
            sc = sc[sc["S_DSB"].notna() & (sc["S_DSB"] >= 0.85)]
        if filters:
            for col, val in filters.items():
                if col in sc.columns:
                    sc = sc[sc[col] >= val]
        return sc.nlargest(top_k, "PenScore").reset_index(drop=True)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_axes_live(
        self, accession: str, walker_motif_override: bool | None = None
    ) -> AxisScores:
        """Compute all axes live for a single accession (lazy imports)."""
        from pen_score.axes import cargo, deliv, dsb, energy, mature, prog, spec

        return AxisScores(
            S_DSB=dsb.score(accession),
            S_Spec=spec.score(accession),
            S_Cargo=cargo.score(accession),
            S_Deliv=deliv.score(accession),
            S_Immuno=None,  # requires netMHCpan -- run offline via script 14
            S_Prog=prog.score(accession),
            S_Mature=mature.score(accession),
            S_Energy=energy.score(accession, walker_motif_override=walker_motif_override),
        )

    @staticmethod
    def _generate_reasoning(axes: AxisScores, weights: dict[str, float]) -> list[str]:
        """Generate human-readable reasoning bullets for a scoring result."""
        _AXES = [
            "S_DSB",
            "S_Spec",
            "S_Cargo",
            "S_Deliv",
            "S_Immuno",
            "S_Prog",
            "S_Mature",
            "S_Energy",
        ]
        _LABELS: dict[str, str] = {
            "S_DSB": "DSB safety",
            "S_Spec": "Specificity",
            "S_Cargo": "Cargo capacity",
            "S_Deliv": "Deliverability",
            "S_Immuno": "Immunogenicity",
            "S_Prog": "Programmability",
            "S_Mature": "Tech maturity",
            "S_Energy": "Energy independence",
        }
        bullets: list[str] = []
        for ax in _AXES:
            v: float | None = getattr(axes, ax, None)
            label: str = _LABELS[ax]
            w: float = weights.get(ax, 0.0)
            if v is None:
                bullets.append(f"[missing] {label} (axis not computed; excluded from composite)")
            elif v >= 0.75:
                bullets.append(f"[strength] {label}: {v:.4f} (weight={w:.2f})")
            elif v <= 0.35:
                bullets.append(f"[weakness] {label}: {v:.4f} (weight={w:.2f})")
        return bullets

    @staticmethod
    def _default_weights() -> dict[str, float]:
        return {
            "S_DSB": 0.24,
            "S_Spec": 0.14,
            "S_Cargo": 0.19,
            "S_Deliv": 0.19,
            "S_Immuno": 0.09,
            "S_Prog": 0.05,
            "S_Mature": 0.05,
            "S_Energy": 0.05,
        }

    @staticmethod
    def _composite(
        axes: AxisScores,
        weights: dict[str, float],
        exclude_axes: list[str] | None = None,
    ) -> tuple[float | None, list[str]]:
        """Weighted composite PenScore; returns (score, list_of_missing_axes).

        Parameters
        ----------
        axes:
            Axis scores (None = not computed).
        weights:
            Per-axis weights from the use-case profile.
        exclude_axes:
            Axis names to exclude entirely.  Their weights are redistributed
            proportionally across the remaining axes.  Excluded axes are NOT
            reported in the missing list.
        """
        excluded: frozenset[str] = frozenset(exclude_axes or [])
        total_w = 0.0
        total_s = 0.0
        missing: list[str] = []
        for ax, w in weights.items():
            if ax in excluded:
                continue  # explicitly excluded - skip without recording as missing
            val: float | None = getattr(axes, ax, None)
            if val is None:
                missing.append(ax)
            else:
                total_w += w
                total_s += w * val
        if total_w == 0.0:
            return None, missing
        # Renormalise: missing + excluded axes redistribute weight to available axes
        return round(total_s / total_w, 4), missing
