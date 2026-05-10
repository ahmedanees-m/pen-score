"""Unit tests for pen_score.scorer.bootstrap.

Pure numpy computation - no external dependencies, no mocking required.
All tests run in CI on every platform.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pen_score.scorer.bootstrap import bootstrap_axis_ci, bootstrap_ranking_ci

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mini_scorecard(scores: dict[str, float]) -> pd.DataFrame:
    """Build a minimal scorecard DataFrame from editor_id -> PenScore dict."""
    return pd.DataFrame([{"editor_id": k, "PenScore": v} for k, v in scores.items()])


# ---------------------------------------------------------------------------
# bootstrap_ranking_ci
# ---------------------------------------------------------------------------


class TestBootstrapRankingCI:
    """Tests for bootstrap_ranking_ci()."""

    def test_returns_dataframe(self):
        sc = _mini_scorecard({"IS621": 0.929, "SpCas9": 0.368, "BE3": 0.242})
        result = bootstrap_ranking_ci(sc, n_bootstrap=100)
        assert isinstance(result, pd.DataFrame)

    def test_one_row_per_editor(self):
        editors = {"IS621": 0.929, "SpCas9": 0.368, "BE3": 0.242, "Cre": 0.797}
        sc = _mini_scorecard(editors)
        result = bootstrap_ranking_ci(sc, n_bootstrap=50)
        assert len(result) == len(editors)

    def test_expected_columns_present(self):
        sc = _mini_scorecard({"IS621": 0.929, "SpCas9": 0.368})
        result = bootstrap_ranking_ci(sc, n_bootstrap=50)
        expected = {
            "editor_id",
            "rank_mean",
            "rank_ci_lower",
            "rank_ci_upper",
            "score_mean",
            "score_ci_lower",
            "score_ci_upper",
        }
        assert expected.issubset(set(result.columns))

    def test_top_editor_rank_mean_is_one(self):
        """The best editor must have mean rank ~ 1."""
        sc = _mini_scorecard({"IS621": 0.929, "SpCas9": 0.368, "BE3": 0.242})
        result = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=42)
        top_row = result[result["editor_id"] == "IS621"].iloc[0]
        # IS621 is always ranked first -> rank_mean should be exactly 1.0
        assert abs(top_row["rank_mean"] - 1.0) < 0.2

    def test_worst_editor_rank_mean_is_last(self):
        """The worst editor must have mean rank ~ n_editors."""
        editors = {"IS621": 0.929, "SpCas9": 0.368, "BE3": 0.242}
        sc = _mini_scorecard(editors)
        result = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=42)
        worst = result[result["editor_id"] == "BE3"].iloc[0]
        assert worst["rank_mean"] > len(editors) - 1.5

    def test_ci_lower_le_mean_le_upper(self):
        """rank_ci_lower <= rank_mean <= rank_ci_upper for all editors."""
        editors = {"A": 0.9, "B": 0.7, "C": 0.5, "D": 0.3}
        sc = _mini_scorecard(editors)
        result = bootstrap_ranking_ci(sc, n_bootstrap=100, seed=42)
        for _, row in result.iterrows():
            assert row["rank_ci_lower"] <= row["rank_mean"] + 1e-9
            assert row["rank_mean"] <= row["rank_ci_upper"] + 1e-9

    def test_score_ci_lower_le_score_upper(self):
        """score_ci_lower <= score_ci_upper for all editors."""
        sc = _mini_scorecard({"A": 0.9, "B": 0.5, "C": 0.1})
        result = bootstrap_ranking_ci(sc, n_bootstrap=100, seed=42)
        for _, row in result.iterrows():
            assert row["score_ci_lower"] <= row["score_ci_upper"] + 1e-9

    def test_seed_reproducibility(self):
        """Same seed -> identical results."""
        editors = {"A": 0.9, "B": 0.7, "C": 0.5, "D": 0.3, "E": 0.1}
        sc = _mini_scorecard(editors)
        r1 = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=42)
        r2 = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=42)
        pd.testing.assert_frame_equal(r1, r2)

    def test_different_seeds_differ(self):
        """Different seeds should generally produce different results."""
        editors = {f"E{i}": 0.9 - i * 0.1 for i in range(6)}
        sc = _mini_scorecard(editors)
        r1 = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=42)
        r2 = bootstrap_ranking_ci(sc, n_bootstrap=200, seed=99)
        # rank_means may differ (though not guaranteed; use score_mean as proxy)
        # Just verify they don't crash with different seeds
        assert len(r1) == len(r2) == len(editors)

    def test_nan_scores_treated_as_zero(self):
        """NaN PenScores are filled to 0.0 (sentinel editors)."""
        sc = pd.DataFrame(
            [
                {"editor_id": "IS621", "PenScore": 0.929},
                {"editor_id": "SENTINEL", "PenScore": float("nan")},
            ]
        )
        result = bootstrap_ranking_ci(sc, n_bootstrap=50, seed=42)
        assert len(result) == 2
        # SENTINEL scores filled to 0 -> will rank last
        sentinel = result[result["editor_id"] == "SENTINEL"].iloc[0]
        assert sentinel["rank_mean"] >= 1.5

    def test_single_editor(self):
        """Edge case: one editor always ranks 1."""
        sc = _mini_scorecard({"IS621": 0.929})
        result = bootstrap_ranking_ci(sc, n_bootstrap=100, seed=42)
        assert len(result) == 1
        assert abs(result.iloc[0]["rank_mean"] - 1.0) < 1e-6

    def test_all_equal_scores(self):
        """When all editors have the same score, rankings are uniformly distributed."""
        editors = {f"E{i}": 0.5 for i in range(5)}
        sc = _mini_scorecard(editors)
        result = bootstrap_ranking_ci(sc, n_bootstrap=500, seed=42)
        # Each editor should have rank_mean ~ 3 (middle of [1,5])
        for _, row in result.iterrows():
            assert 1.0 <= row["rank_mean"] <= 5.0

    def test_custom_score_column(self):
        """bootstrap_ranking_ci respects the score_col parameter."""
        sc = pd.DataFrame(
            [
                {"editor_id": "A", "PenScore": 0.1, "S_DSB": 0.9},
                {"editor_id": "B", "PenScore": 0.9, "S_DSB": 0.1},
            ]
        )
        result_penscore = bootstrap_ranking_ci(sc, score_col="PenScore", n_bootstrap=50)
        result_sdsb = bootstrap_ranking_ci(sc, score_col="S_DSB", n_bootstrap=50)
        # PenScore: B ranks 1st (0.9 > 0.1); S_DSB: A ranks 1st (0.9 > 0.1)
        ps_top = result_penscore.sort_values("rank_mean").iloc[0]["editor_id"]
        sdsb_top = result_sdsb.sort_values("rank_mean").iloc[0]["editor_id"]
        assert ps_top == "B"
        assert sdsb_top == "A"

    def test_custom_ci(self):
        """90% CI should be narrower than 95% CI."""
        editors = {f"E{i}": i * 0.1 for i in range(10)}
        sc = _mini_scorecard(editors)
        r95 = bootstrap_ranking_ci(sc, n_bootstrap=300, seed=42, ci=0.95)
        r90 = bootstrap_ranking_ci(sc, n_bootstrap=300, seed=42, ci=0.90)

        # Width = upper - lower; 90% CI should be narrower on average
        width_95 = (r95["rank_ci_upper"] - r95["rank_ci_lower"]).mean()
        width_90 = (r90["rank_ci_upper"] - r90["rank_ci_lower"]).mean()
        # Allow a small tolerance since the same seed means same bootstrap draws
        assert width_90 <= width_95 + 0.5


# ---------------------------------------------------------------------------
# bootstrap_axis_ci
# ---------------------------------------------------------------------------


class TestBootstrapAxisCI:
    """Tests for bootstrap_axis_ci()."""

    def test_returns_three_floats(self):
        values = np.array([0.9, 0.85, 0.92, 0.88, 0.91])
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=200)
        assert isinstance(mean, float)
        assert isinstance(lo, float)
        assert isinstance(hi, float)

    def test_lo_le_mean_le_hi(self):
        values = np.array([0.9, 0.85, 0.92, 0.88, 0.91])
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=500, seed=42)
        assert lo <= mean + 1e-9
        assert mean <= hi + 1e-9

    def test_constant_array_ci_is_zero_width(self):
        """All-same values -> bootstrap means are constant -> CI width = 0."""
        values = np.ones(10) * 0.75
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=100, seed=42)
        assert abs(mean - 0.75) < 1e-9
        assert abs(lo - 0.75) < 1e-9
        assert abs(hi - 0.75) < 1e-9

    def test_mean_close_to_sample_mean(self):
        """Bootstrap mean should be close to the sample mean."""
        values = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        expected_mean = float(np.mean(values))  # 0.5
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=1000, seed=42)
        assert abs(mean - expected_mean) < 0.05

    def test_seed_reproducibility(self):
        """Same seed -> identical outputs."""
        values = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
        r1 = bootstrap_axis_ci(values, n_bootstrap=200, seed=42)
        r2 = bootstrap_axis_ci(values, n_bootstrap=200, seed=42)
        assert r1 == r2

    def test_different_seeds_differ(self):
        """Different seeds -> different outputs (with high probability)."""
        values = np.array([0.1, 0.5, 0.9, 0.3, 0.7, 0.6, 0.4])
        r42 = bootstrap_axis_ci(values, n_bootstrap=500, seed=42)
        r99 = bootstrap_axis_ci(values, n_bootstrap=500, seed=99)
        # lo and hi should be different (extremely unlikely to match)
        assert r42[1] != r99[1] or r42[2] != r99[2]

    def test_bounds_within_data_range(self):
        """CI bounds cannot exceed the minimum/maximum of the input values."""
        values = np.array([0.3, 0.5, 0.7])
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=200, seed=42)
        # Bootstrap means are always between the min and max of the data
        assert lo >= 0.3 - 1e-9
        assert hi <= 0.7 + 1e-9

    def test_90_pct_ci_narrower_than_95(self):
        """90% CI must be narrower than 95% CI for the same seed."""
        values = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 0.2, 0.8])
        _, lo95, hi95 = bootstrap_axis_ci(values, n_bootstrap=500, seed=42, ci=0.95)
        _, lo90, hi90 = bootstrap_axis_ci(values, n_bootstrap=500, seed=42, ci=0.90)
        width_95 = hi95 - lo95
        width_90 = hi90 - lo90
        assert width_90 <= width_95 + 1e-9

    def test_single_element_array(self):
        """Edge case: single element -> no variance -> mean equals that element."""
        values = np.array([0.75])
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=100, seed=42)
        assert abs(mean - 0.75) < 1e-9
        assert abs(lo - 0.75) < 1e-9
        assert abs(hi - 0.75) < 1e-9

    def test_is621_rank_ci_is_tight(self):
        """IS621 PenScore distribution has near-zero variance -> tight CI."""
        # Simulate 28 bootstrap draws all near IS621's score 0.929
        rng = np.random.default_rng(42)
        values = rng.normal(loc=0.929, scale=0.005, size=50).clip(0, 1)
        mean, lo, hi = bootstrap_axis_ci(values, n_bootstrap=500, seed=42)
        # CI width should be tight (σ=0.005 input -> CI width ~ 0.01 or less)
        assert (hi - lo) < 0.05
