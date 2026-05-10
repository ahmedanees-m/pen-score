"""Extended Scorer tests - covers scorecard-path branch, get_scorecard, select_editor."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest


@pytest.fixture
def sample_scorecard() -> pd.DataFrame:
    """5-editor scorecard simulating the full public scorecard schema."""
    return pd.DataFrame(
        [
            {
                "editor_id": "IS621",
                "canonical_accession": "A0A7C9VKZ0",
                "S_DSB": 1.0,
                "S_Spec": 0.99,
                "S_Cargo": 1.0,
                "S_Deliv": 0.80,
                "S_Immuno": 0.76,
                "S_Prog": 1.0,
                "S_Mature": 0.20,
                "S_Energy": 1.0,
                "PenScore": 0.957,
            },
            {
                "editor_id": "ISCro4",
                "canonical_accession": "D2TGM5",
                "S_DSB": 1.0,
                "S_Spec": 0.98,
                "S_Cargo": 1.0,
                "S_Deliv": 0.82,
                "S_Immuno": 0.77,
                "S_Prog": 1.0,
                "S_Mature": 0.05,
                "S_Energy": 1.0,
                "PenScore": 0.930,
            },
            {
                "editor_id": "CAST_VK",
                "canonical_accession": "A0A2X3M8B1",
                "S_DSB": 1.0,
                "S_Spec": 0.60,
                "S_Cargo": 0.90,
                "S_Deliv": 0.70,
                "S_Immuno": 0.55,
                "S_Prog": 0.0,
                "S_Mature": 0.08,
                "S_Energy": 0.0,
                "PenScore": 0.720,
            },
            {
                "editor_id": "SpCas9",
                "canonical_accession": "Q99ZW2",
                "S_DSB": 0.0,
                "S_Spec": 0.75,
                "S_Cargo": 0.5,
                "S_Deliv": 0.30,
                "S_Immuno": 0.0,
                "S_Prog": 1.0,
                "S_Mature": 1.0,
                "S_Energy": 0.0,
                "PenScore": 0.402,
            },
            {
                "editor_id": "PE2",
                "canonical_accession": "A0A2X3M8B0",
                "S_DSB": 0.0,
                "S_Spec": 0.85,
                "S_Cargo": 0.40,
                "S_Deliv": 0.25,
                "S_Immuno": 0.42,
                "S_Prog": 1.0,
                "S_Mature": 0.60,
                "S_Energy": 0.0,
                "PenScore": 0.450,
            },
        ]
    )


class TestScorerFromScorecard:
    """Tests for Scorer.score_editor when a pre-computed scorecard is loaded."""

    def test_score_from_scorecard_returns_result(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()
        scorer._scorecard = sample_scorecard

        result = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
        )

        assert result.accession == "A0A7C9VKZ0"
        assert result.axes.S_DSB == 1.0
        assert result.axes.S_Spec == 0.99
        assert result.pen_score is not None

    def test_score_from_scorecard_uses_weights(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()
        scorer._scorecard = sample_scorecard

        result = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
        )

        assert result.use_case == "human_therapeutic_aav_insertion"
        assert result.pen_score is not None
        assert 0.0 <= result.pen_score <= 1.0

    def test_score_from_scorecard_penscore_close_to_reference(self, sample_scorecard):
        """IS621 PenScore from the pre-computed row should be roughly correct when
        re-derived from axis values + default weights."""
        from pen_score.api import Scorer

        scorer = Scorer.load()
        scorer._scorecard = sample_scorecard

        result = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
        )

        # The scorer re-derives from axes using weights, so it won't exactly match
        # the pre-computed PenScore column (which may have slightly different weights),
        # but it should be within a reasonable range for IS621.
        assert result.pen_score is not None
        assert result.pen_score > 0.70

    def test_score_editor_exclude_axes_invalid_raises(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()
        scorer._scorecard = sample_scorecard

        with pytest.raises(ValueError, match="Unknown axis"):
            scorer.score_editor(
                accession="A0A7C9VKZ0",
                exclude_axes=["S_Invalid"],
            )

    def test_score_editor_exclude_s_mature(self, sample_scorecard):
        """Excluding S_Mature for IS621 (S_Mature=0.20) should raise the PenScore."""
        from pen_score.api import Scorer

        scorer = Scorer.load()
        scorer._scorecard = sample_scorecard

        full = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
        )
        bio = scorer.score_editor(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
            exclude_axes=["S_Mature"],
        )

        # IS621 has S_Mature=0.20 which slightly drags down the full score.
        # Excluding it should give equal or higher score.
        assert bio.pen_score is not None
        assert full.pen_score is not None
        assert bio.pen_score >= full.pen_score - 0.001  # allow float rounding


class TestGetScorecardFromParquet:
    """Tests for Scorer.get_scorecard when scorecard_path is given."""

    def test_loads_from_parquet_when_available(self, sample_scorecard, tmp_path):
        from pen_score.api import Scorer

        parquet_path = tmp_path / "test_scorecard.parquet"
        sample_scorecard.to_parquet(parquet_path, index=False)

        scorer = Scorer.load(scorecard_path=parquet_path)
        sc = scorer.get_scorecard(use_case="human_therapeutic_aav_insertion")

        assert len(sc) == len(sample_scorecard)
        assert set(sc.columns) == set(sample_scorecard.columns)

    def test_scorecard_has_correct_editors(self, sample_scorecard, tmp_path):
        from pen_score.api import Scorer

        parquet_path = tmp_path / "sc.parquet"
        sample_scorecard.to_parquet(parquet_path, index=False)

        scorer = Scorer.load(scorecard_path=parquet_path)
        sc = scorer.get_scorecard()

        assert "IS621" in sc["editor_id"].values
        assert "SpCas9" in sc["editor_id"].values


class TestSelectEditorMocked:
    """Tests for Scorer.select_editor using a mocked get_scorecard."""

    def test_top_k_returns_correct_count(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=3,
            )

        assert len(result) == 3

    def test_sorted_by_penscore(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=5,
            )

        scores = result["PenScore"].tolist()
        assert scores == sorted(scores, reverse=True)

    def test_top1_is_highest_penscore(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=5,
            )

        assert result.iloc[0]["editor_id"] == "IS621"

    def test_require_dsb_free_filters_nucleases(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=10,
                require_dsb_free=True,
            )

        editor_ids = result["editor_id"].tolist()
        assert "SpCas9" not in editor_ids
        assert "PE2" not in editor_ids

    def test_require_dsb_free_keeps_is621(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=10,
                require_dsb_free=True,
            )

        assert "IS621" in result["editor_id"].values
        assert "ISCro4" in result["editor_id"].values

    def test_additional_filters_applied(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=10,
                filters={"S_Prog": 1.0},
            )

        for _, row in result.iterrows():
            assert row["S_Prog"] >= 1.0, f"{row['editor_id']} should have S_Prog >= 1.0"

    def test_index_reset(self, sample_scorecard):
        from pen_score.api import Scorer

        scorer = Scorer.load()

        with patch.object(scorer, "get_scorecard", return_value=sample_scorecard):
            result = scorer.select_editor(
                use_case="human_therapeutic_aav_insertion",
                top_k=3,
            )

        assert list(result.index) == list(range(len(result)))


class TestScorerDefaultWeights:
    """Tests for Scorer._default_weights."""

    def test_weights_sum_to_one(self):
        from pen_score.api import Scorer

        weights = Scorer._default_weights()
        assert abs(sum(weights.values()) - 1.0) < 1e-9

    def test_all_axes_present(self):
        from pen_score.api import Scorer

        weights = Scorer._default_weights()
        expected = {"S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature", "S_Energy"}
        assert set(weights.keys()) == expected

    def test_s_dsb_highest_weight(self):
        from pen_score.api import Scorer

        weights = Scorer._default_weights()
        assert weights["S_DSB"] == max(weights.values())
