"""Unit tests for pen_score.api - Scorer class methods that require no external deps."""

from __future__ import annotations

from pathlib import Path

import pytest


class TestScorerComposite:
    """Tests for Scorer._composite (mirrors composite.compute_pen_score)."""

    def test_full_axis_set(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=1.0,
            S_Spec=0.7,
            S_Cargo=1.0,
            S_Deliv=0.95,
            S_Immuno=0.8,
            S_Prog=1.0,
            S_Mature=0.3,
            S_Energy=1.0,
        )
        score, missing = Scorer._composite(axes, default_weights)
        assert missing == []
        assert score is not None
        assert 0.85 <= score <= 0.95

    def test_missing_axis_excluded(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=1.0,
            S_Spec=1.0,
            S_Cargo=1.0,
            S_Deliv=1.0,
            S_Immuno=None,
            S_Prog=1.0,
            S_Mature=1.0,
        )
        score, missing = Scorer._composite(axes, default_weights)
        assert "S_Immuno" in missing
        assert score is not None
        assert abs(score - 1.0) < 0.001

    def test_all_zeros_returns_zero(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=0.0,
            S_Spec=0.0,
            S_Cargo=0.0,
            S_Deliv=0.0,
            S_Immuno=0.0,
            S_Prog=0.0,
            S_Mature=0.0,
            S_Energy=0.0,
        )
        score, missing = Scorer._composite(axes, default_weights)
        assert missing == []
        assert score == 0.0


class TestScorerReasoning:
    """Tests for Scorer._generate_reasoning - pure logic, no external deps."""

    def test_strength_bullet_generated(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(S_DSB=1.0, S_Spec=0.9, S_Cargo=1.0)
        bullets = Scorer._generate_reasoning(axes, default_weights)
        strength_bullets = [b for b in bullets if "[strength]" in b]
        assert len(strength_bullets) >= 2

    def test_weakness_bullet_generated(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(S_DSB=0.0, S_Spec=0.3, S_Cargo=0.2)
        bullets = Scorer._generate_reasoning(axes, default_weights)
        weakness_bullets = [b for b in bullets if "[weakness]" in b]
        assert len(weakness_bullets) >= 2

    def test_missing_axis_bullet_generated(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(S_DSB=None, S_Spec=None)
        bullets = Scorer._generate_reasoning(axes, default_weights)
        missing_bullets = [b for b in bullets if "[missing]" in b]
        assert len(missing_bullets) >= 2

    def test_mid_range_produces_no_bullet(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=0.5,
            S_Spec=0.5,
            S_Cargo=0.5,
            S_Deliv=0.5,
            S_Immuno=0.5,
            S_Prog=0.5,
            S_Mature=0.5,
            S_Energy=0.5,
        )
        bullets = Scorer._generate_reasoning(axes, default_weights)
        assert len(bullets) == 0

    def test_returns_list(self, default_weights):
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores()
        result = Scorer._generate_reasoning(axes, default_weights)
        assert isinstance(result, list)


class TestExcludeAxes:
    """Tests for the exclude_axes parameter of Scorer.score_editor."""

    def test_exclude_s_mature_raises_penscore(self):
        """Excluding S_Mature should raise PenScore when S_Mature was 0.

        IS622 has S_Mature=0.0; excluding it should increase PenScore.
        We simulate this with a direct _composite call.
        """
        from pen_score.api import AxisScores, Scorer

        # IS622-like axes: strong biophysical performance, S_Mature=0
        axes = AxisScores(
            S_DSB=1.0,
            S_Spec=0.99,
            S_Cargo=1.0,
            S_Deliv=0.95,
            S_Immuno=0.76,
            S_Prog=1.0,
            S_Mature=0.0,
            S_Energy=1.0,
        )
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
        full_score, _ = Scorer._composite(axes, weights)
        bio_score, _ = Scorer._composite(axes, weights, exclude_axes=["S_Mature"])
        assert bio_score is not None
        assert full_score is not None
        assert bio_score > full_score, (
            f"Biophysical score ({bio_score}) should exceed full score ({full_score}) "
            "when S_Mature=0.0 is excluded"
        )
        assert bio_score > 0.90, f"IS622-like biophysical PenScore should be >0.90, got {bio_score}"

    def test_exclude_axes_weight_renormalises_to_one(self):
        """After excluding axes, the effective weight fraction used should be <1 but score is
        renormalised to [0,1]."""
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=1.0,
            S_Spec=1.0,
            S_Cargo=1.0,
            S_Deliv=1.0,
            S_Immuno=1.0,
            S_Prog=1.0,
            S_Mature=1.0,
            S_Energy=1.0,
        )
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
        # All axes = 1.0, so regardless of which axes are excluded the score is still 1.0
        score, missing = Scorer._composite(axes, weights, exclude_axes=["S_Mature", "S_Immuno"])
        assert score is not None
        assert abs(score - 1.0) < 1e-9, (
            f"All-ones axes should give 1.0 regardless of exclusion, got {score}"
        )
        assert missing == [], f"No missing axes expected, got {missing}"

    def test_exclude_axes_invalid_raises(self):
        """Passing an unrecognised axis name should raise ValueError."""
        from pen_score.api import Scorer

        scorer = Scorer.load()
        with pytest.raises(ValueError, match="Unknown axis"):
            scorer.score_editor(
                accession="A0A2X3M8B0",
                exclude_axes=["S_Nonexistent"],
            )

    def test_exclude_all_returns_none(self):
        """Excluding all 8 axes should return None PenScore."""
        from pen_score.api import AxisScores, Scorer

        axes = AxisScores(
            S_DSB=1.0,
            S_Spec=1.0,
            S_Cargo=1.0,
            S_Deliv=1.0,
            S_Immuno=1.0,
            S_Prog=1.0,
            S_Mature=1.0,
            S_Energy=1.0,
        )
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
        all_axes = [
            "S_DSB",
            "S_Spec",
            "S_Cargo",
            "S_Deliv",
            "S_Immuno",
            "S_Prog",
            "S_Mature",
            "S_Energy",
        ]
        score, _ = Scorer._composite(axes, weights, exclude_axes=all_axes)
        assert score is None, f"Excluding all axes should return None, got {score}"


class TestSelectEditorDSBFree:
    """Tests for the require_dsb_free flag - uses scorecard parquet."""

    _SC_PATH = Path.home() / "pen-stack/data/pen-score/scorecards/public_scorecard.parquet"

    @pytest.fixture(scope="class")
    def scorer(self):
        if not self._SC_PATH.exists():
            pytest.skip("Public scorecard parquet not available in this environment")
        from pen_score.api import Scorer

        return Scorer(scorecard_path=self._SC_PATH)

    def test_dsb_free_flag_filters_correctly(self, scorer):
        result = scorer.select_editor(
            use_case="human_therapeutic_aav_insertion",
            top_k=10,
            require_dsb_free=True,
        )
        for _, row in result.iterrows():
            assert row["S_DSB"] >= 0.85, (
                f"{row['editor_id']} has S_DSB={row['S_DSB']} but require_dsb_free=True"
            )

    def test_dsb_free_reduces_candidate_pool(self, scorer):
        all_editors = scorer.select_editor(
            use_case="human_therapeutic_aav_insertion",
            top_k=28,
            require_dsb_free=False,
        )
        dsb_free = scorer.select_editor(
            use_case="human_therapeutic_aav_insertion",
            top_k=28,
            require_dsb_free=True,
        )
        assert len(dsb_free) < len(all_editors)
        assert 10 <= len(dsb_free) <= 15

    def test_top_k_respects_limit(self, scorer):
        result = scorer.select_editor(
            use_case="large_cargo_integration",
            top_k=3,
            require_dsb_free=True,
        )
        assert len(result) <= 3
