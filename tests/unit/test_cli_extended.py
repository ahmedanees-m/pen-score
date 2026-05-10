"""Extended CLI tests using Click's CliRunner - no external scoring deps required."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
from click.testing import CliRunner


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def mock_scoring_result():
    """A ScoringResult-like mock for Scorer.score_editor."""
    from pen_score.api import AxisScores, ScoringResult

    return ScoringResult(
        accession="A0A7C9VKZ0",
        editor_id="IS621",
        use_case="human_therapeutic_aav_insertion",
        axes=AxisScores(
            S_DSB=1.0,
            S_Spec=0.99,
            S_Cargo=1.0,
            S_Deliv=0.80,
            S_Immuno=0.76,
            S_Prog=1.0,
            S_Mature=0.20,
            S_Energy=1.0,
        ),
        pen_score=0.957,
        reasoning=["[strength] DSB safety: 1.0000 (weight=0.24)"],
        axes_missing=[],
    )


@pytest.fixture
def mock_scorecard() -> pd.DataFrame:
    """Minimal 3-row scorecard DataFrame."""
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
        ]
    )


class TestScoreEditorCommand:
    """Tests for the `pen-score score-editor` CLI command."""

    def test_score_editor_basic(self, runner, mock_scoring_result):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.score_editor.return_value = mock_scoring_result

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(main, ["score-editor", "A0A7C9VKZ0"])

        assert result.exit_code == 0, f"CLI failed: {result.output}"
        assert "A0A7C9VKZ0" in result.output

    def test_score_editor_with_use_case(self, runner, mock_scoring_result):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.score_editor.return_value = mock_scoring_result

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-editor", "A0A7C9VKZ0", "--use-case", "large_cargo_integration"],
            )

        assert result.exit_code == 0
        mock_scorer.score_editor.assert_called_once_with(
            accession="A0A7C9VKZ0",
            use_case="large_cargo_integration",
            exclude_axes=None,
        )

    def test_score_editor_with_exclude_axes(self, runner, mock_scoring_result):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.score_editor.return_value = mock_scoring_result

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-editor", "A0A7C9VKZ0", "--exclude-axes", "S_Mature"],
            )

        assert result.exit_code == 0
        mock_scorer.score_editor.assert_called_once_with(
            accession="A0A7C9VKZ0",
            use_case="human_therapeutic_aav_insertion",
            exclude_axes=["S_Mature"],
        )

    def test_score_editor_multiple_exclude_axes(self, runner, mock_scoring_result):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.score_editor.return_value = mock_scoring_result

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-editor", "A0A7C9VKZ0", "--exclude-axes", "S_Mature,S_Immuno"],
            )

        assert result.exit_code == 0
        call_kwargs = mock_scorer.score_editor.call_args
        exclude = call_kwargs.kwargs.get("exclude_axes") or call_kwargs.args[2]
        assert "S_Mature" in exclude
        assert "S_Immuno" in exclude

    def test_score_editor_saves_parquet(self, runner, mock_scoring_result, tmp_path):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.score_editor.return_value = mock_scoring_result
        out_path = str(tmp_path / "out.parquet")

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-editor", "A0A7C9VKZ0", "--output", out_path],
            )

        assert result.exit_code == 0
        assert "Saved" in result.output

    def test_score_editor_no_accession_fails(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["score-editor"])
        assert result.exit_code != 0


class TestScoreAllEditorsCommand:
    """Tests for `pen-score score-all-editors` CLI command."""

    def test_score_all_editors_default_output(self, runner, mock_scorecard, tmp_path):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.get_scorecard.return_value = mock_scorecard

        with runner.isolated_filesystem(temp_dir=tmp_path):
            with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
                result = runner.invoke(main, ["score-all-editors"])

        assert result.exit_code == 0
        assert "scorecard" in result.output.lower() or "editors" in result.output.lower()

    def test_score_all_editors_custom_output(self, runner, mock_scorecard, tmp_path):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.get_scorecard.return_value = mock_scorecard
        out_path = str(tmp_path / "custom.parquet")

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-all-editors", "--output", out_path],
            )

        assert result.exit_code == 0

    def test_score_all_editors_custom_use_case(self, runner, mock_scorecard, tmp_path):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.get_scorecard.return_value = mock_scorecard
        out_path = str(tmp_path / "out.parquet")

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["score-all-editors", "--output", out_path, "--use-case", "large_cargo_integration"],
            )

        assert result.exit_code == 0
        mock_scorer.get_scorecard.assert_called_once_with(use_case="large_cargo_integration")


class TestSelectCommand:
    """Tests for `pen-score select` CLI command."""

    def test_select_requires_use_case(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["select"])
        assert result.exit_code != 0

    def test_select_basic(self, runner, mock_scorecard):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.select_editor.return_value = mock_scorecard

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["select", "--use-case", "human_therapeutic_aav_insertion"],
            )

        assert result.exit_code == 0

    def test_select_top_k_passed(self, runner, mock_scorecard):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.select_editor.return_value = mock_scorecard.head(3)

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                ["select", "--use-case", "human_therapeutic_aav_insertion", "--top-k", "3"],
            )

        assert result.exit_code == 0
        mock_scorer.select_editor.assert_called_once_with(
            use_case="human_therapeutic_aav_insertion",
            top_k=3,
            require_dsb_free=False,
        )

    def test_select_require_dsb_free_flag(self, runner, mock_scorecard):
        from pen_score.cli import main

        dsb_free_sc = mock_scorecard[mock_scorecard["S_DSB"] >= 0.85].copy()
        mock_scorer = MagicMock()
        mock_scorer.select_editor.return_value = dsb_free_sc

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                [
                    "select",
                    "--use-case",
                    "human_therapeutic_aav_insertion",
                    "--require-dsb-free",
                ],
            )

        assert result.exit_code == 0
        mock_scorer.select_editor.assert_called_once_with(
            use_case="human_therapeutic_aav_insertion",
            top_k=5,
            require_dsb_free=True,
        )

    def test_select_saves_parquet(self, runner, mock_scorecard, tmp_path):
        from pen_score.cli import main

        mock_scorer = MagicMock()
        mock_scorer.select_editor.return_value = mock_scorecard
        out_path = str(tmp_path / "ranked.parquet")

        with patch("pen_score.api.Scorer.load", return_value=mock_scorer):
            result = runner.invoke(
                main,
                [
                    "select",
                    "--use-case",
                    "human_therapeutic_aav_insertion",
                    "--output",
                    out_path,
                ],
            )

        assert result.exit_code == 0
        assert "Saved" in result.output


class TestCLIVersion:
    """Test --version flag."""

    def test_version_flag(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "pen-score" in result.output.lower() or "version" in result.output.lower()


class TestCLIHelp:
    """Test --help flags for all commands."""

    def test_main_help(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "score-editor" in result.output

    def test_score_editor_help(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["score-editor", "--help"])
        assert result.exit_code == 0
        assert "accession" in result.output.lower()

    def test_score_all_editors_help(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["score-all-editors", "--help"])
        assert result.exit_code == 0

    def test_select_help(self, runner):
        from pen_score.cli import main

        result = runner.invoke(main, ["select", "--help"])
        assert result.exit_code == 0
        assert "use-case" in result.output
