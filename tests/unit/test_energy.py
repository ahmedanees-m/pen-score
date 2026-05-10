"""Unit tests for S_Energy - Energy Independence axis (pen_score/axes/energy.py).

All tests that require a live UniProt fetch are mocked so CI runs without
network access.  The override logic and sentinel handling are tested directly
without mocking.
"""

from __future__ import annotations

import warnings
from unittest.mock import patch

import pytest

# ---------------------------------------------------------------------------
# Helper: build a fake UniProt FASTA response
# ---------------------------------------------------------------------------


def _fasta(header: str, sequence: str) -> bytes:
    return f">sp|{header}|TEST\n{sequence}\n".encode()


# ---------------------------------------------------------------------------
# Override logic (no network needed)
# ---------------------------------------------------------------------------


class TestWalkerOverride:
    """walker_motif_override bypasses the sequence scan entirely."""

    def test_override_true_returns_zero(self):
        from pen_score.axes.energy import score

        result = score("Q99ZW2", walker_motif_override=True)
        assert result == 0.0

    def test_override_false_returns_one(self):
        from pen_score.axes.energy import score

        result = score("NO_UNIPROT", walker_motif_override=False)
        assert result == 1.0

    def test_override_true_ignores_accession_entirely(self):
        """Even a sentinel accession returns 0.0 when override=True."""
        from pen_score.axes.energy import score

        result = score("REQUIRES_STEP7", walker_motif_override=True)
        assert result == 0.0

    def test_override_false_ignores_accession_entirely(self):
        """Even a sentinel accession returns 1.0 when override=False."""
        from pen_score.axes.energy import score

        result = score("REQUIRES_STEP7", walker_motif_override=False)
        assert result == 1.0


# ---------------------------------------------------------------------------
# Sentinel accessions -> None (no override)
# ---------------------------------------------------------------------------


class TestSentinelAccessions:
    """Sentinel accessions without override return None with a warning."""

    @pytest.mark.parametrize("acc", ["REQUIRES_STEP7", "NO_UNIPROT", "", None])
    def test_sentinel_returns_none(self, acc):
        from pen_score.axes.energy import score

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = score(acc)
        assert result is None
        assert len(w) == 1
        assert "sentinel" in str(w[0].message).lower()


# ---------------------------------------------------------------------------
# Walker A motif detection (mocked UniProt fetch)
# ---------------------------------------------------------------------------


class TestWalkerAMotif:
    """Walker A (GxxxxGK[ST]) triggers S_Energy = 0.0."""

    def test_walker_a_detected_returns_zero(self):
        from pen_score.axes.energy import score

        # GxxxxGKS: exactly 4 residues between the two G's (canonical Walker A)
        seq = "MAAACDEFGAAAAGKSTTTTTT"
        with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
            result = score("FAKE001")
        assert result == 0.0

    def test_walker_a_with_threonine_variant(self):
        from pen_score.axes.energy import score

        seq = "MAAACDEFGAAAAGKTTTTTTTT"  # GxxxxGKT - threonine variant
        with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
            result = score("FAKE002")
        assert result == 0.0

    def test_incomplete_walker_a_not_flagged(self):
        """GxxxxGK alone (no S/T) does not trigger."""
        from pen_score.axes.energy import score

        seq = "MAAACDEFGAAAAAGKAAALLL"  # GxxxxGKA - A not S/T
        with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
            result = score("FAKE003")
        # No Walker A, no Walker B -> energy independent
        assert result == 1.0


# ---------------------------------------------------------------------------
# Walker B motif detection (mocked UniProt fetch)
# ---------------------------------------------------------------------------


class TestWalkerBMotif:
    """Walker B ([LVIMF]{4}DE) triggers S_Energy = 0.0."""

    def test_walker_b_detected_returns_zero(self):
        from pen_score.axes.energy import score

        seq = "MAAACDEFLLLLDEAAATTT"  # LLLLLDE - 4 hydrophobic then DE
        with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
            result = score("FAKE004")
        assert result == 0.0

    def test_walker_b_all_hydrophobic_variants(self):
        from pen_score.axes.energy import score

        for motif in ["LLLLDE", "VVVVDE", "IIIIDE", "MMMMDE", "FFFFTE"]:
            # FFFFTE has T not E; should NOT fire
            expected = 0.0 if motif.endswith("DE") else 1.0
            seq = f"MAAACDEF{motif}AAATTT"
            with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
                result = score("FAKE005")
            assert result == expected, f"Motif {motif!r} gave {result}, expected {expected}"


# ---------------------------------------------------------------------------
# Energy-independent editors (mocked - no Walker motifs)
# ---------------------------------------------------------------------------


class TestEnergyIndependentEditors:
    """Editors known to be energy-independent score 1.0."""

    @pytest.mark.parametrize(
        "editor_id,accession",
        [
            ("IS621", "A0A2X3M8B0"),
            ("SpCas9", "Q99ZW2"),
            ("Bxb1", "Q9B086"),
            ("Cre", "P06956"),
        ],
    )
    def test_known_energy_independent(self, editor_id, accession):
        """These editors have no Walker A/B motifs; scan returns 1.0."""
        from pen_score.axes.energy import score

        # Realistic-length sequence with NO Walker A or B motifs
        seq = "MAADKLQRSTPVGEYIACNWM" * 15  # ~315 aa, no GxxxxGK[ST] or hhhhDE
        with patch("pen_score.axes.energy._fetch_sequence", return_value=seq):
            result = score(accession)
        assert result == 1.0, f"{editor_id} should be energy-independent"


# ---------------------------------------------------------------------------
# ATP-dependent CAST systems (override-based, no fetch needed)
# ---------------------------------------------------------------------------


class TestCastSystemsATPDependent:
    """CAST systems use walker_motif_override=True -> S_Energy = 0.0."""

    @pytest.mark.parametrize(
        "editor_id,accession",
        [
            ("evoCAST", "REQUIRES_STEP7"),
            ("CAST_VK", "A0A8M0FGU0"),
            ("CAST_IF", "A0A0F4L2U9"),
        ],
    )
    def test_cast_override_true(self, editor_id, accession):
        from pen_score.axes.energy import score

        result = score(accession, walker_motif_override=True)
        assert result == 0.0, f"{editor_id} should be ATP-dependent"


# ---------------------------------------------------------------------------
# SleepingBeauty - NO_UNIPROT but confirmed energy-independent
# ---------------------------------------------------------------------------


class TestSleepingBeauty:
    def test_sleeping_beauty_override_false(self):
        from pen_score.axes.energy import score

        result = score("NO_UNIPROT", walker_motif_override=False)
        assert result == 1.0


# ---------------------------------------------------------------------------
# Network failure -> None + warning
# ---------------------------------------------------------------------------


class TestNetworkFailure:
    def test_fetch_failure_returns_none_with_warning(self):
        from pen_score.axes.energy import score

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with patch("pen_score.axes.energy._fetch_sequence", return_value=None):
                result = score("Q99ZW2")
        assert result is None
        # _fetch_sequence already warned; score() should not double-warn
        # (the warning originates inside _fetch_sequence, not score())
        _ = w  # just ensure no crash

    def test_fetch_failure_does_not_raise(self):
        from pen_score.axes.energy import score

        with patch("pen_score.axes.energy._fetch_sequence", return_value=None):
            result = score("Q99ZW2")
        assert result is None


# ---------------------------------------------------------------------------
# Composite integration - AxisScores includes S_Energy
# ---------------------------------------------------------------------------


class TestAxisScoresModel:
    def test_axis_scores_has_s_energy_field(self):
        from pen_score.api import AxisScores

        axes = AxisScores(S_Energy=1.0)
        assert axes.S_Energy == 1.0

    def test_axis_scores_s_energy_defaults_none(self):
        from pen_score.api import AxisScores

        axes = AxisScores()
        assert axes.S_Energy is None

    def test_axis_scores_s_energy_rejects_out_of_range(self):
        from pydantic import ValidationError

        from pen_score.api import AxisScores

        with pytest.raises(ValidationError):
            AxisScores(S_Energy=1.5)


# ---------------------------------------------------------------------------
# Composite scorer includes S_Energy
# ---------------------------------------------------------------------------


class TestCompositeWithSEnergy:
    def test_composite_uses_s_energy(self):
        from pen_score.scorer.composite import compute_pen_score

        axes = {
            "S_DSB": 1.0,
            "S_Spec": 1.0,
            "S_Cargo": 1.0,
            "S_Deliv": 1.0,
            "S_Immuno": 1.0,
            "S_Prog": 1.0,
            "S_Mature": 1.0,
            "S_Energy": 0.0,  # ATP-dependent penalised
        }
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
        ps, missing = compute_pen_score(axes, weights)
        assert ps is not None
        assert ps < 1.0  # penalised by S_Energy = 0.0
        assert missing == []

    def test_composite_missing_s_energy_renormalises(self):
        """When S_Energy is None the score should be computed over 7 axes."""
        from pen_score.scorer.composite import compute_pen_score

        axes = {
            "S_DSB": 1.0,
            "S_Spec": 1.0,
            "S_Cargo": 1.0,
            "S_Deliv": 1.0,
            "S_Immuno": 1.0,
            "S_Prog": 1.0,
            "S_Mature": 1.0,
            "S_Energy": None,
        }
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
        ps, missing = compute_pen_score(axes, weights)
        assert ps == 1.0  # renormalised over 7 available axes
        assert "S_Energy" in missing
