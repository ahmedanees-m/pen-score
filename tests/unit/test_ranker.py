"""Unit tests for pen_score.scorer.ranker - pure pandas, no external deps."""

from __future__ import annotations

import pandas as pd
import pytest


@pytest.fixture
def sample_scorecard() -> pd.DataFrame:
    """Minimal scorecard DataFrame with 5 editors for ranking tests."""
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
                "PenScore": 0.930,
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
                "PenScore": 0.402,
            },
            {
                "editor_id": "evoCAST",
                "canonical_accession": "REQUIRES_STEP7_evoCAST",
                "S_DSB": 1.0,
                "S_Spec": 0.60,
                "S_Cargo": 0.90,
                "S_Deliv": None,  # sentinel
                "S_Immuno": None,
                "S_Prog": 0.0,
                "S_Mature": 0.10,
                "PenScore": 0.720,
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
                "PenScore": 0.450,
            },
        ]
    )


@pytest.fixture
def default_weights() -> dict[str, float]:
    return {
        "S_DSB": 0.24,
        "S_Spec": 0.14,
        "S_Cargo": 0.19,
        "S_Deliv": 0.19,
        "S_Immuno": 0.09,
        "S_Prog": 0.05,
        "S_Mature": 0.05,
    }


class TestRankEditors:
    """Tests for rank_editors()."""

    def test_sorted_by_penscore_descending(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(sample_scorecard, use_case="test_case", weights=default_weights)

        scores = result["PenScore"].tolist()
        assert scores == sorted(scores, reverse=True), "Should be sorted descending by PenScore"

    def test_rank_column_starts_at_one(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(sample_scorecard, use_case="test_case", weights=default_weights)

        assert result["rank"].iloc[0] == 1
        assert result["rank"].iloc[-1] == len(result)

    def test_use_case_column_set(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(
            sample_scorecard, use_case="human_therapeutic_aav_insertion", weights=default_weights
        )

        assert all(result["use_case"] == "human_therapeutic_aav_insertion")

    def test_reasoning_column_added(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(sample_scorecard, use_case="test_case", weights=default_weights)

        assert "reasoning" in result.columns
        assert all(isinstance(r, list) for r in result["reasoning"])

    def test_top_k_limits_rows(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(
            sample_scorecard, use_case="test_case", weights=default_weights, top_k=3
        )

        assert len(result) == 3
        assert result["PenScore"].iloc[0] >= result["PenScore"].iloc[-1]

    def test_filter_by_min_s_dsb(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(
            sample_scorecard,
            use_case="dsb_free_case",
            weights=default_weights,
            filters={"S_DSB": 0.85},
        )

        for _, row in result.iterrows():
            assert row["S_DSB"] >= 0.85, f"{row['editor_id']} should have S_DSB >= 0.85"

    def test_filter_removes_nucleases(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(
            sample_scorecard,
            use_case="dsb_free_case",
            weights=default_weights,
            filters={"S_DSB": 0.85},
        )

        editor_ids = result["editor_id"].tolist()
        assert "SpCas9" not in editor_ids
        assert "PE2" not in editor_ids

    def test_filter_on_nonexistent_column_ignored(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        # Should not raise even if column doesn't exist
        result = rank_editors(
            sample_scorecard,
            use_case="test",
            weights=default_weights,
            filters={"S_NonExistentAxis": 0.5},
        )

        assert len(result) == len(sample_scorecard)

    def test_original_dataframe_not_mutated(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        original_cols = set(sample_scorecard.columns)
        rank_editors(sample_scorecard, use_case="test", weights=default_weights)

        assert set(sample_scorecard.columns) == original_cols, "Original DF should not be mutated"

    def test_top1_is_is621(self, sample_scorecard, default_weights):
        from pen_score.scorer.ranker import rank_editors

        result = rank_editors(sample_scorecard, use_case="test", weights=default_weights)

        assert result.iloc[0]["editor_id"] == "IS621"


class TestBuildReasoning:
    """Tests for _build_reasoning()."""

    @pytest.fixture
    def weights(self) -> dict[str, float]:
        return {
            "S_DSB": 0.24,
            "S_Spec": 0.14,
            "S_Cargo": 0.19,
            "S_Deliv": 0.19,
            "S_Immuno": 0.09,
            "S_Prog": 0.05,
            "S_Mature": 0.05,
        }

    def test_excellent_for_high_scores(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series(
            {
                "S_DSB": 1.0,
                "S_Spec": 0.99,
                "S_Cargo": 1.0,
                "S_Deliv": 0.95,
                "S_Immuno": 0.90,
                "S_Prog": 1.0,
                "S_Mature": 0.80,
            }
        )
        reasons = _build_reasoning(row, weights)

        excellent = [r for r in reasons if r.startswith("Excellent")]
        assert len(excellent) >= 5

    def test_moderate_for_mid_scores(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series({"S_DSB": 0.6, "S_Spec": 0.55, "S_Cargo": 0.65})
        reasons = _build_reasoning(row, weights)

        moderate = [r for r in reasons if r.startswith("Moderate")]
        assert len(moderate) >= 2

    def test_poor_for_low_scores(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series({"S_DSB": 0.0, "S_Spec": 0.3})
        reasons = _build_reasoning(row, weights)

        poor = [r for r in reasons if r.startswith("Poor")]
        assert len(poor) >= 2

    def test_na_for_none_values(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series({"S_DSB": None, "S_Deliv": None})
        reasons = _build_reasoning(row, weights)

        na_reasons = [r for r in reasons if "N/A" in r]
        assert len(na_reasons) >= 2

    def test_returns_list(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series({"S_DSB": 0.5})
        result = _build_reasoning(row, weights)

        assert isinstance(result, list)

    def test_weight_shown_in_reason(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series({"S_DSB": 1.0})
        reasons = _build_reasoning(row, weights)

        dsb_reason = next((r for r in reasons if "DSB avoidance" in r), None)
        assert dsb_reason is not None
        assert "0.24" in dsb_reason  # weight shown

    def test_all_axes_covered(self, weights):
        from pen_score.scorer.ranker import _build_reasoning

        row = pd.Series(
            {
                "S_DSB": 0.5,
                "S_Spec": 0.5,
                "S_Cargo": 0.5,
                "S_Deliv": 0.5,
                "S_Immuno": 0.5,
                "S_Prog": 0.5,
                "S_Mature": 0.5,
            }
        )
        reasons = _build_reasoning(row, weights)

        # All 7 axes in weights should be represented
        assert len(reasons) == len(weights)
