"""Unit tests for score axis computations (no external dependencies required)."""

from __future__ import annotations

from pen_score.axes import cargo, deliv, prog


class TestSCargo:
    """S_Cargo is pure computation from YAML - no external deps."""

    def test_megabase_scores_one(self):
        # IS621: 1,000,000 bp -> log10(1e6)/log10(1e6) = 1.0
        s = cargo.score("IS621")
        assert s is not None
        assert abs(s - 1.0) < 0.01

    def test_pe2_scores_low(self):
        # PE2: 200 bp -> log10(200)/6 ~ 0.38
        s = cargo.score("PE2")
        assert s is not None
        assert 0.35 <= s <= 0.42

    def test_cas9_scores_mid(self):
        # SpCas9+HDR: 1000 bp -> log10(1000)/6 = 3/6 = 0.5
        s = cargo.score("SpCas9")
        assert s is not None
        assert 0.48 <= s <= 0.52

    def test_unknown_accession_returns_none(self):
        s = cargo.score("NONEXISTENT_XXXX")
        assert s is None

    def test_all_editors_score_between_0_and_1(self):
        from pen_score.data.loader import load_editor_universe

        for ed in load_editor_universe():
            s = cargo.score(ed.id)
            if s is not None:
                assert 0.0 <= s <= 1.0, f"{ed.id}: S_Cargo={s} out of range"


class TestSDeliv:
    """S_Deliv requires only the sigmoid formula - no external deps."""

    def test_compact_editor_scores_high(self):
        # IS621 (~300 aa) -> sigmoid(0.005 * (300 - 900)) ~ 0.998
        s = deliv.score("IS621", total_aa=300)
        assert s is not None
        assert s >= 0.95

    def test_large_editor_scores_low(self):
        # SpCas9 (1368 aa) -> sigmoid(0.005 * (1368 - 900)) ~ 0.10
        s = deliv.score("SpCas9", total_aa=1368)
        assert s is not None
        assert s <= 0.30

    def test_boundary_900aa_scores_half(self):
        # 900 aa -> sigmoid(0) = 0.5 exactly
        s = deliv.score("boundary", total_aa=900)
        assert s is not None
        assert abs(s - 0.5) < 0.01

    def test_scores_in_range(self):
        for aa in [100, 300, 600, 900, 1200, 1500, 2000]:
            s = deliv.score("test", total_aa=aa)
            assert s is not None
            assert 0.0 <= s <= 1.0

    def test_monotone_decreasing_with_size(self):
        """Larger editors must score lower on deliverability."""
        aas = [200, 500, 900, 1200, 1800]
        scores = [deliv.score("x", total_aa=aa) for aa in aas]
        for i in range(len(scores) - 1):
            assert scores[i] > scores[i + 1], (
                f"Deliverability not monotone: aa={aas[i]} score={scores[i]} "
                f"aa={aas[i + 1]} score={scores[i + 1]}"
            )


class TestSProg:
    """S_Prog is a binary lookup from editor_universe.yaml - no external deps."""

    def test_spcas9_is_rna_guided(self):
        s = prog.score("SpCas9")
        assert s == 1.0

    def test_is621_is_rna_guided(self):
        s = prog.score("IS621")
        assert s == 1.0

    def test_cre_is_not_programmable(self):
        s = prog.score("Cre")
        assert s == 0.0

    def test_bxb1_is_not_programmable(self):
        s = prog.score("Bxb1")
        assert s == 0.0

    def test_unknown_returns_none(self):
        s = prog.score("NONEXISTENT_XXXX")
        assert s is None

    def test_all_scores_are_valid_values(self):
        """All S_Prog values must be 0.0, 0.5, or 1.0 (discrete tiers)."""
        from pen_score.data.loader import load_editor_universe

        valid = {0.0, 0.5, 1.0}
        for ed in load_editor_universe():
            s = prog.score(ed.id)
            if s is not None:
                assert s in valid, f"{ed.id}: S_Prog={s} not in {valid}"


class TestComposite:
    """Tests for the composite PenScore computation function."""

    def test_all_axes_present_returns_weighted_sum(self, mock_axis_scores, default_weights):
        from pen_score.scorer.composite import compute_pen_score

        score, missing = compute_pen_score(mock_axis_scores, default_weights)
        assert missing == []
        assert score is not None
        # Manual (v0.1.1 weights): 1.0*0.24 + 0.7*0.14 + 1.0*0.19 + 0.95*0.19 + 0.8*0.09
        #   + 1.0*0.05 + 0.3*0.05 + 1.0*0.05
        # = 0.24 + 0.098 + 0.19 + 0.1805 + 0.072 + 0.05 + 0.015 + 0.05 = 0.8955
        assert abs(score - 0.8955) < 0.005

    def test_missing_axis_excluded_from_denominator(self, default_weights):
        from pen_score.scorer.composite import compute_pen_score

        # S_Immuno missing - denominator reduces from 1.0 to 0.90
        scores = {
            "S_DSB": 1.0,
            "S_Spec": 1.0,
            "S_Cargo": 1.0,
            "S_Deliv": 1.0,
            "S_Immuno": None,
            "S_Prog": 1.0,
            "S_Mature": 1.0,
        }
        score, missing = compute_pen_score(scores, default_weights)
        assert "S_Immuno" in missing
        # All available axes = 1.0 -> renormalised score should still be 1.0
        assert score is not None
        assert abs(score - 1.0) < 0.001

    def test_all_missing_returns_none(self, default_weights):
        from pen_score.scorer.composite import compute_pen_score

        _AX = [
            "S_DSB",
            "S_Spec",
            "S_Cargo",
            "S_Deliv",
            "S_Immuno",
            "S_Prog",
            "S_Mature",
            "S_Energy",
        ]
        scores = dict.fromkeys(_AX)
        score, missing = compute_pen_score(scores, default_weights)
        assert score is None
        assert len(missing) == 8

    def test_score_bounded_between_0_and_1(self, default_weights):
        from pen_score.scorer.composite import compute_pen_score

        for v in [0.0, 0.25, 0.5, 0.75, 1.0]:
            _AX2 = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature"]
            scores = dict.fromkeys(_AX2, v)
            score, _ = compute_pen_score(scores, default_weights)
            assert score is not None
            assert abs(score - v) < 0.001, f"Uniform score={v} -> PenScore={score}"

    def test_output_rounded_to_4dp(self, mock_axis_scores, default_weights):
        from pen_score.scorer.composite import compute_pen_score

        score, _ = compute_pen_score(mock_axis_scores, default_weights)
        assert score is not None
        # Should be rounded to 4 decimal places
        assert score == round(score, 4)
