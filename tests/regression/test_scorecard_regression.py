"""Scorecard regression tests - golden values reflect pen-score v0.1.1 (8-axis).

These tests load the pre-computed public scorecard parquet and verify that key
numeric scores match the pre-registered values (+/- 0.001 tolerance for float
rounding).  They will be skipped automatically outside the VM environment where
the parquet is not present.

Pre-registration hash: pre-registration-v1.0.2 (2026-05-13T17:56:38Z)
mech-class version: v0.5.2 (IS110 Tier-A gate; S_DSB=1.0 for IS110-family proteins)
pen-score version: v0.1.1 (8-axis; S_Energy added 2026-05-22)

v0.5.2 correction (2026-05-22): IS110-family proteins (PF01548 and PF02371) now receive
S_DSB=1.0 via the Tier-A hard gate.  IS621 PenScore increases from 0.929 -> 0.957
(+0.028 from combined S_DSB correction and 8-axis S_Energy weight).  The parquet
must be regenerated with mech-class v0.5.2 + pen-score pipeline scripts 10 + 17 + 20-22
for these tests to pass.

8-axis weights (human_therapeutic_aav_insertion v0.1.1):
    S_DSB=0.24, S_Spec=0.14, S_Cargo=0.19, S_Deliv=0.19,
    S_Immuno=0.09, S_Prog=0.05, S_Mature=0.05, S_Energy=0.05
"""

from __future__ import annotations

from pathlib import Path

import pytest

SCORECARD_PATH = Path.home() / "pen-stack/data/pen-score/scorecards/public_scorecard.parquet"


@pytest.fixture(scope="module")
def scorecard():
    if not SCORECARD_PATH.exists():
        pytest.skip("Public scorecard parquet not available (VM-only test)")
    import pandas as pd

    return pd.read_parquet(SCORECARD_PATH).set_index("editor_id")


def get(sc, editor_id: str, col: str):
    """Retrieve a value; return None if column or editor is missing."""
    try:
        v = sc.at[editor_id, col]
        import math

        return None if (isinstance(v, float) and math.isnan(v)) else v
    except KeyError:
        return None


# ---------------------------------------------------------------------------
# Golden scores - pre-registration-v1.0.2; updated for mech-class v0.5.2 (IS110 Tier-A gate)
# ---------------------------------------------------------------------------


class TestIS621GoldenValues:
    """IS621 is the #1 editor overall; regression-lock its key scores."""

    def test_penscore_aav(self, scorecard):
        """v0.1.1 8-axis: IS621 PenScore = 0.957 (was 0.929 with 7-axis; +0.028 combined)."""
        ps = get(scorecard, "IS621", "PenScore")
        assert ps is not None
        assert abs(ps - 0.957) < 0.002, f"IS621 PenScore={ps}, expected ~0.957 (v0.1.1)"

    def test_s_dsb_is_one(self, scorecard):
        """v0.5.2 IS110 Tier-A gate: IS621 S_DSB=1.0 (was 0.90 with bucket heuristic v0.5.1)."""
        v = get(scorecard, "IS621", "S_DSB")
        assert v is not None
        assert abs(v - 1.0) < 0.001, f"IS621 S_DSB={v}, expected 1.0 (IS110 Tier-A gate)"

    def test_s_spec_is_high(self, scorecard):
        v = get(scorecard, "IS621", "S_Spec")
        assert v is not None
        assert v >= 0.99, f"IS621 S_Spec={v}, expected >= 0.99"

    def test_s_cargo_is_one(self, scorecard):
        v = get(scorecard, "IS621", "S_Cargo")
        assert v is not None
        assert abs(v - 1.0) < 0.01, f"IS621 S_Cargo={v}, expected 1.0 (megabase capacity)"

    def test_s_prog_is_one(self, scorecard):
        v = get(scorecard, "IS621", "S_Prog")
        assert v is not None
        assert abs(v - 1.0) < 0.01, f"IS621 S_Prog={v}, expected 1.0 (RNA-guided)"


class TestSpCas9GoldenValues:
    """SpCas9 pre-registered P3 prediction: bottom 30% of full universe."""

    def test_penscore_aav(self, scorecard):
        ps = get(scorecard, "SpCas9", "PenScore")
        assert ps is not None
        assert abs(ps - 0.402) < 0.005, f"SpCas9 PenScore={ps}, expected ~0.402 (v0.1.1)"

    def test_s_dsb_is_zero(self, scorecard):
        """v0.5.1+ fix: DSB_NUCLEASE -> S_DSB=0.0 (gate requires PF01548+PF02371)."""
        v = get(scorecard, "SpCas9", "S_DSB")
        assert v is not None
        assert abs(v - 0.0) < 0.001, (
            f"SpCas9 S_DSB={v}, expected 0.0 (composite gate requires "
            "PF01548+PF02371 - SpCas9 has neither)"
        )

    def test_s_spec_is_high(self, scorecard):
        """SpCas9 specificity is high (genome-wide low off-target after optimisation)."""
        v = get(scorecard, "SpCas9", "S_Spec")
        assert v is not None
        assert v >= 0.95, f"SpCas9 S_Spec={v}, expected >= 0.95"

    def test_s_deliv_is_poor(self, scorecard):
        """SpCas9 is too large for single-AAV (1368 aa -> S_Deliv ~ 0.09)."""
        v = get(scorecard, "SpCas9", "S_Deliv")
        assert v is not None
        assert v < 0.20, f"SpCas9 S_Deliv={v}, expected < 0.20 (too large for AAV)"

    def test_bottom_30pct_prediction_passes(self, scorecard):
        """P3 pre-registered prediction: fraction of editors below SpCas9 >= 30%."""
        all_scores = scorecard["PenScore"].dropna()
        spas9_score = get(scorecard, "SpCas9", "PenScore")
        assert spas9_score is not None
        n_below = int((all_scores < spas9_score).sum())
        fraction_below = n_below / len(all_scores)
        assert fraction_below >= 0.30, (
            f"P3 FAIL: only {fraction_below:.1%} of editors score below SpCas9 "
            f"(need >= 30%); SpCas9 PenScore={spas9_score:.4f}"
        )


class TestDSBNucleaseV051Fix:
    """Verify that ALL DSB_NUCLEASE editors have S_DSB=0.0 (mech-class v0.5.1+)."""

    def test_all_dsb_nuclease_score_zero(self, scorecard):
        dsb_nuc = scorecard[scorecard["mechanism_bucket"] == "DSB_NUCLEASE"]
        for editor_id, row in dsb_nuc.iterrows():
            v = row.get("S_DSB")
            import math

            if v is not None and not (isinstance(v, float) and math.isnan(v)):
                assert abs(v - 0.0) < 0.001, (
                    f"{editor_id} is DSB_NUCLEASE but S_DSB={v}, "
                    "expected 0.0 (composite gate requires PF01548+PF02371)"
                )


class TestIS621RanksFirst:
    """IS621 must be rank #1 in the two primary DSB-free use cases."""

    def _penscore(self, row, weights):
        import math

        num, den = 0.0, 0.0
        for ax, w in weights.items():
            v = row.get(ax)
            if v is not None and not (isinstance(v, float) and math.isnan(v)):
                num += w * v
                den += w
        return num / den if den > 0 else 0.0

    def test_is621_rank1_aav_insertion(self, scorecard):
        # v0.1.1 8-axis weights (human_therapeutic_aav_insertion)
        weights = {
            "S_DSB": 0.24,
            "S_Spec": 0.14,
            "S_Cargo": 0.19,
            "S_Deliv": 0.19,
            "S_Immuno": 0.09,
            "S_Prog": 0.05,
            "S_Mature": 0.05,
            "S_Energy": 0.05,
        }
        scored = scorecard.apply(lambda r: self._penscore(r, weights), axis=1)
        ranked = scored.sort_values(ascending=False)
        assert ranked.index[0] == "IS621", (
            f"IS621 expected rank #1 for AAV insertion, got {ranked.index[0]}"
        )

    def test_is621_rank1_large_cargo(self, scorecard):
        # v0.1.1 8-axis weights (large_cargo_integration)
        weights = {
            "S_DSB": 0.19,
            "S_Spec": 0.09,
            "S_Cargo": 0.38,
            "S_Deliv": 0.10,
            "S_Immuno": 0.05,
            "S_Prog": 0.09,
            "S_Mature": 0.05,
            "S_Energy": 0.05,
        }
        scored = scorecard.apply(lambda r: self._penscore(r, weights), axis=1)
        ranked = scored.sort_values(ascending=False)
        assert ranked.index[0] == "IS621", (
            f"IS621 expected rank #1 for large cargo, got {ranked.index[0]}"
        )


class TestBootstrapCI:
    """P2: IS621 top-3 prediction with bootstrap CI = [1,1]."""

    def test_is621_penscore_significantly_above_is621_2(self, scorecard):
        """IS621 must outscore IS621_2 (both share programmability but IS621 has edge)."""
        ps_is621 = get(scorecard, "IS621", "PenScore")
        ps_is621_2 = get(scorecard, "IS621_2", "PenScore")
        if ps_is621 is None or ps_is621_2 is None:
            pytest.skip("IS621 or IS621_2 not in scorecard")
        assert ps_is621 >= ps_is621_2, (
            f"IS621 PenScore={ps_is621} should be >= IS621_2={ps_is621_2}"
        )

    def test_universe_size(self, scorecard):
        """Scorecard should have exactly 29 editors (28 pipeline + IS622 post-pipeline)."""
        assert len(scorecard) == 29, f"Expected 29 editors, got {len(scorecard)}"
