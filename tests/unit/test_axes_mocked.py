"""Mock-based tests for external-dependency axis modules.

Covers the code paths in spec.py, immuno.py, mature.py, dsb.py, deliv.py, and
d7_homolumo.py that require BWA, netMHCpan, mhcflurry, mech-class, or xtb.

All external binaries and optional packages are mocked out, so these tests
run in vanilla CI without installing any extra system packages.

Mocking strategy
----------------
- ``pysam``           : injected as a MagicMock into sys.modules
- ``subprocess.run``  : patched via unittest.mock.patch
- ``requests.get``    : patched via unittest.mock.patch
- ``mech_class``      : injected as a MagicMock namespace into sys.modules
- ``mhcflurry``       : injected as a MagicMock namespace into sys.modules
- UniProt helpers     : patched at the call-site in pen_score.axes.deliv / pen_score.axes.mature
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_sam_line(nm: int) -> str:
    """Return a minimal SAM alignment line with the given NM tag."""
    fields = [
        "query",
        "0",
        "chr1",
        "100",
        "60",
        "20M",
        "*",
        "0",
        "0",
        "GGGTGGGGGGAGTTTGCTCC",
        "*",
        f"NM:i:{nm}",
        "MD:Z:20",
    ]
    return "\t".join(fields)


def _sam_output(*nm_values: int) -> str:
    """Build a fake BWA stdout string with header + alignment lines."""
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    for nm in nm_values:
        lines.append(_make_sam_line(nm))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# S_Spec  (spec.py)
# ---------------------------------------------------------------------------


class TestSSpecMocked:
    """Tests for S_Spec - mocks pysam + subprocess so BWA path is exercised."""

    def test_no_protospacer_returns_none(self):
        from pen_score.axes import spec

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = spec.score("Q99ZW2", protospacer=None)
        assert result is None
        assert any("protospacer" in str(warning.message) for warning in w)

    def test_pysam_not_installed_returns_none(self):
        """When pysam is absent from sys.modules the ImportError branch fires."""
        # Ensure pysam is absent
        pysam_backup = sys.modules.pop("pysam", None)
        try:
            from pen_score.axes import spec

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")
            assert result is None
            assert any("pysam" in str(warning.message) for warning in w)
        finally:
            if pysam_backup is not None:
                sys.modules["pysam"] = pysam_backup

    def test_bwa_two_hits_nm_le3_returns_valid_score(self):
        """Two NM<=3 hits -> off_target_count=2 -> sigmoid score in [0, 1].

        With only 2 hits in a 3.2 Gbp genome the ratio is ~6e-7, so sigmoid
        returns a value very close to (but possibly equal to) 1.0 after rounding
        to 4 dp.  We verify the return is a valid non-None float in [0, 1].
        """
        mock_pysam = MagicMock()
        mock_result = MagicMock()
        mock_result.stdout = _sam_output(0, 2)  # NM=0 and NM=2 both count

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", return_value=mock_result):
                from pen_score.axes import spec

                s = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        assert s is not None
        assert 0.0 <= s <= 1.0
        assert s > 0.9  # very few off-targets -> high specificity

    def test_bwa_zero_hits_scores_near_one(self):
        """Zero off-target hits -> log10(1e-10) floor -> sigmoid -> near 1.0."""
        mock_pysam = MagicMock()
        mock_result = MagicMock()
        mock_result.stdout = _sam_output()  # header only, no alignments

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", return_value=mock_result):
                from pen_score.axes import spec

                s = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        assert s is not None
        assert s > 0.9

    def test_bwa_many_hits_scores_lower_than_zero_hits(self):
        """More off-target hits -> lower specificity score (monotone property).

        10 000 NM=0 hits in a 3.2 Gbp genome gives ratio ~ 0.003, which still
        maps to sigmoid ~ 0.993 (the formula is calibrated for transposon-scale
        genomes).  We verify the score is strictly less than the zero-hit score
        rather than asserting an absolute threshold.
        """
        mock_pysam = MagicMock()

        mock_many = MagicMock()
        mock_many.stdout = _sam_output(*([0] * 10_000))

        mock_zero = MagicMock()
        mock_zero.stdout = _sam_output()  # no alignments

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", return_value=mock_many):
                from pen_score.axes import spec

                s_many = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")
            with patch("subprocess.run", return_value=mock_zero):
                s_zero = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        assert s_many is not None and s_zero is not None
        assert s_many < s_zero  # more hits -> lower score (monotone)

    def test_nm_tag_gt3_not_counted(self):
        """Alignments with NM>3 must not contribute to off_target_count."""
        mock_pysam = MagicMock()
        mock_result_low = MagicMock()
        mock_result_low.stdout = _sam_output(4, 5, 6)  # all NM>3 - none counted

        mock_result_none = MagicMock()
        mock_result_none.stdout = _sam_output()  # no alignments

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", return_value=mock_result_low):
                from pen_score.axes import spec

                s_no_nm = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

            with patch("subprocess.run", return_value=mock_result_none):
                s_zero = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        # Both should be near 1.0 (floor at 1e-10 ratio)
        assert s_no_nm is not None and s_zero is not None
        assert abs(s_no_nm - s_zero) < 0.001  # NM>3 lines contribute nothing

    def test_bwa_exception_returns_none(self):
        """If subprocess raises, the except block fires and None is returned."""
        mock_pysam = MagicMock()

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", side_effect=FileNotFoundError("bwa not found")):
                from pen_score.axes import spec

                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    result = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        assert result is None
        assert any("BWA" in str(warning.message) for warning in w)

    def test_score_rounded_to_4dp(self):
        """Return value must be rounded to 4 decimal places."""
        mock_pysam = MagicMock()
        mock_result = MagicMock()
        mock_result.stdout = _sam_output(1, 1, 1)  # 3 hits

        with patch.dict("sys.modules", {"pysam": mock_pysam}):
            with patch("subprocess.run", return_value=mock_result):
                from pen_score.axes import spec

                s = spec.score("Q99ZW2", protospacer="GGGTGGGGGGAGTTTGCTCC")

        assert s is not None
        assert s == round(s, 4)

    def test_sigmoid_helper(self):
        """Internal _sigmoid is a valid sigmoid (monotone, (0,1) range)."""
        from pen_score.axes.spec import _sigmoid

        assert abs(_sigmoid(0.0) - 0.5) < 1e-9
        assert _sigmoid(100.0) > 0.999
        assert _sigmoid(-100.0) < 0.001
        # Monotone
        prev = _sigmoid(-10.0)
        for x in [-5.0, 0.0, 5.0, 10.0]:
            cur = _sigmoid(x)
            assert cur > prev
            prev = cur


# ---------------------------------------------------------------------------
# S_Immuno  (immuno.py)
# ---------------------------------------------------------------------------


class TestSImmuno:
    """Tests for S_Immuno - uses cached parquet fast path and mocked mhcflurry."""

    def test_cached_parquet_returns_score(self, tmp_path):
        """Fast path: if parquet exists and contains the accession, return its score."""
        parquet_file = tmp_path / "immuno_scores.parquet"
        df = pd.DataFrame(
            {"canonical_accession": ["Q99ZW2", "A0A2X3M8B0"], "S_Immuno": [0.05, 0.82]}
        )
        df.to_parquet(parquet_file)

        from pen_score.axes import immuno

        s = immuno.score("Q99ZW2", cached_parquet=parquet_file)
        assert s is not None
        assert abs(s - 0.05) < 0.001

    def test_cached_parquet_second_row(self, tmp_path):
        parquet_file = tmp_path / "immuno_scores.parquet"
        df = pd.DataFrame(
            {"canonical_accession": ["Q99ZW2", "A0A2X3M8B0"], "S_Immuno": [0.05, 0.82]}
        )
        df.to_parquet(parquet_file)

        from pen_score.axes import immuno

        s = immuno.score("A0A2X3M8B0", cached_parquet=parquet_file)
        assert s is not None
        assert abs(s - 0.82) < 0.001

    def test_cached_parquet_accession_missing_falls_through_to_none(self, tmp_path):
        """If accession not in parquet and mhcflurry absent, returns None."""
        parquet_file = tmp_path / "immuno_scores.parquet"
        df = pd.DataFrame({"canonical_accession": ["OTHER"], "S_Immuno": [0.5]})
        df.to_parquet(parquet_file)

        mhcflurry_backup = sys.modules.pop("mhcflurry", None)
        try:
            from pen_score.axes import immuno

            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                result = immuno.score("Q99ZW2", cached_parquet=parquet_file)
            assert result is None
        finally:
            if mhcflurry_backup is not None:
                sys.modules["mhcflurry"] = mhcflurry_backup

    def test_no_parquet_no_mhcflurry_returns_none(self):
        """No parquet and mhcflurry absent -> ImportError branch -> None."""
        mhcflurry_backup = sys.modules.pop("mhcflurry", None)
        try:
            from pen_score.axes import immuno

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = immuno.score("Q99ZW2")
            assert result is None
            assert any("mhcflurry" in str(warning.message) for warning in w)
        finally:
            if mhcflurry_backup is not None:
                sys.modules["mhcflurry"] = mhcflurry_backup

    def test_mocked_mhcflurry_live_path(self):
        """With mhcflurry mocked, the live computation path runs and returns a float."""
        from pen_score.axes import immuno

        # Build a mock presentation predictor that returns a DataFrame
        mock_df = pd.DataFrame({"presentation_score": [0.8, 0.3, 0.6, 0.1, 0.9]})

        mock_predictor_instance = MagicMock()
        mock_predictor_instance.predict.return_value = mock_df

        mock_predictor_cls = MagicMock()
        mock_predictor_cls.load.return_value = mock_predictor_instance

        mock_mhcflurry = ModuleType("mhcflurry")
        mock_mhcflurry.Class1PresentationPredictor = mock_predictor_cls  # type: ignore[attr-defined]

        seq = "ACDEFGHIKLM" * 10  # 110 aa synthetic sequence

        with patch.dict("sys.modules", {"mhcflurry": mock_mhcflurry}):
            result = immuno.score("Q99ZW2", sequence=seq)

        assert result is not None
        assert 0.0 <= result <= 1.0

    def test_mhcflurry_exception_returns_none(self):
        """If the live predictor raises, the except block fires and returns None."""
        from pen_score.axes import immuno

        mock_predictor_cls = MagicMock()
        mock_predictor_cls.load.side_effect = RuntimeError("model weights missing")

        mock_mhcflurry = ModuleType("mhcflurry")
        mock_mhcflurry.Class1PresentationPredictor = mock_predictor_cls  # type: ignore[attr-defined]

        seq = "ACDEFGHIKLM" * 10

        with patch.dict("sys.modules", {"mhcflurry": mock_mhcflurry}):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = immuno.score("Q99ZW2", sequence=seq)

        assert result is None
        assert any("failed" in str(warning.message).lower() for warning in w)

    def test_cached_parquet_nonexistent_path(self):
        """If parquet path is given but does not exist, skip the fast path."""
        from pen_score.axes import immuno

        nonexistent = Path("/tmp/pen_score_nonexistent_12345.parquet")
        # mhcflurry absent -> returns None after falling through
        mhcflurry_backup = sys.modules.pop("mhcflurry", None)
        try:
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                result = immuno.score("Q99ZW2", cached_parquet=nonexistent)
            assert result is None
        finally:
            if mhcflurry_backup is not None:
                sys.modules["mhcflurry"] = mhcflurry_backup


# ---------------------------------------------------------------------------
# S_Mature  (mature.py)
# ---------------------------------------------------------------------------


class TestSMatureMocked:
    """Tests for S_Mature - mocks requests.get for PubMed E-utilities."""

    def _mock_pubmed_response(self, count: int) -> MagicMock:
        """Return a mock requests.Response with the given PubMed hit count."""
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"esearchresult": {"count": str(count)}}
        mock_resp.raise_for_status.return_value = None
        return mock_resp

    def test_query_pubmed_count_returns_int(self):
        """_query_pubmed_count returns the integer from esearchresult.count."""
        from pen_score.axes.mature import _query_pubmed_count

        with patch("requests.get", return_value=self._mock_pubmed_response(1234)):
            result = _query_pubmed_count(["SpCas9", "Cas9"])
        assert result == 1234

    def test_score_with_known_count(self):
        """score() normalises log10(count+1)/log10(max+1) correctly."""
        from pen_score.axes import mature

        max_count = 10_000
        pub_count = 1000
        expected = round(math.log10(pub_count + 1) / math.log10(max_count + 1), 4)

        with patch("requests.get", return_value=self._mock_pubmed_response(pub_count)):
            s = mature.score("Q99ZW2", search_terms=["SpCas9"], max_count=max_count)

        assert s is not None
        assert abs(s - expected) < 1e-4

    def test_zero_citations_returns_zero(self):
        """Zero PubMed hits -> log10(0+1)/log10(max+1) = 0/anything = 0."""
        from pen_score.axes import mature

        with patch("requests.get", return_value=self._mock_pubmed_response(0)):
            s = mature.score("SENTINEL", search_terms=["novel_unknown_editor_xyz"])
        assert s is not None
        assert s == 0.0

    def test_max_count_hits_returns_one(self):
        """If count == max_count, score -> 1.0."""
        from pen_score.axes import mature

        max_count = 5000
        with patch("requests.get", return_value=self._mock_pubmed_response(max_count)):
            s = mature.score("Q99ZW2", search_terms=["SpCas9"], max_count=max_count)
        assert s is not None
        assert abs(s - 1.0) < 1e-4

    def test_no_search_terms_returns_none(self):
        """Empty search_terms triggers warning and returns None."""
        from pen_score.axes import mature

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = mature.score("UNKNOWN_ACC_12345", search_terms=[])
        assert result is None
        assert any(
            "search_terms" in str(warning.message).lower()
            or "references_used_for_pubmed" in str(warning.message)
            for warning in w
        )

    def test_request_failure_returns_none(self):
        """If the HTTP request raises, the except block fires -> None."""
        from pen_score.axes import mature

        with patch("requests.get", side_effect=ConnectionError("network unreachable")):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = mature.score("Q99ZW2", search_terms=["SpCas9"])
        assert result is None
        assert any("PubMed" in str(warning.message) for warning in w)

    def test_score_from_cached_parquet(self, tmp_path):
        """Fast path: use pre-computed parquet rather than querying PubMed."""
        parquet_file = tmp_path / "mature_scores.parquet"
        df = pd.DataFrame({"canonical_accession": ["Q99ZW2"], "S_Mature": [0.97]})
        df.to_parquet(parquet_file)

        from pen_score.axes import mature

        # requests.get should NOT be called if parquet is used
        with patch("requests.get", side_effect=AssertionError("should not call network")):
            s = mature.score("Q99ZW2", cached_parquet=parquet_file)
        assert s is not None
        assert abs(s - 0.97) < 0.001

    def test_search_terms_from_universe(self):
        """If search_terms=None and editor is in universe, terms are auto-resolved."""
        from pen_score.axes import mature
        from pen_score.data.loader import load_editor_universe

        # Find an editor that has references_used_for_pubmed
        editors = [e for e in load_editor_universe() if e.references_used_for_pubmed]
        if not editors:
            pytest.skip("No editors with references_used_for_pubmed in universe")

        ed = editors[0]
        with patch("requests.get", return_value=self._mock_pubmed_response(500)):
            s = mature.score(ed.canonical_accession, search_terms=None)
        assert s is not None

    def test_score_bounded_between_0_and_1(self):
        """All output values must lie in [0, 1]."""
        from pen_score.axes import mature

        for count in [0, 1, 10, 100, 1000, 9999, 10000]:
            with patch("requests.get", return_value=self._mock_pubmed_response(count)):
                s = mature.score("Q99ZW2", search_terms=["SpCas9"], max_count=10000)
            assert s is not None
            assert 0.0 <= s <= 1.0, f"count={count} -> score={s} out of [0,1]"


# ---------------------------------------------------------------------------
# S_DSB  (dsb.py)
# ---------------------------------------------------------------------------


def _make_mech_class_mock(
    p_dsb: float = 0.95,
    p_dsb_free: float = 0.05,
    composite: bool = False,
    tier_a_gate_override: bool = False,
) -> tuple[ModuleType, ModuleType]:
    """Build mech_class and mech_class.api mock modules with a Predictor."""
    mock_pred_result = MagicMock()
    mock_pred_result.tier_a_probabilities = {
        "DSB_NUCLEASE": p_dsb,
        "DSB_FREE_TRANSEST_RECOMBINASE": p_dsb_free,
    }
    mock_pred_result.composite = composite
    mock_pred_result.tier_a_gate_override = tier_a_gate_override

    mock_predictor = MagicMock()
    mock_predictor.predict_from_sequence.return_value = mock_pred_result

    mock_predictor_cls = MagicMock()
    mock_predictor_cls.load.return_value = mock_predictor

    mock_api = ModuleType("mech_class.api")
    mock_api.Predictor = mock_predictor_cls  # type: ignore[attr-defined]

    mock_root = ModuleType("mech_class")

    return mock_root, mock_api


class TestSDSBMocked:
    """Tests for S_DSB - mocks mech-class Predictor."""

    def test_mech_class_not_installed_returns_none(self):
        """ImportError branch fires when mech_class absent -> None + warning."""
        mech_backup = sys.modules.pop("mech_class", None)
        mech_api_backup = sys.modules.pop("mech_class.api", None)
        try:
            from pen_score.axes import dsb

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = dsb.score("Q99ZW2")
            assert result is None
            assert any("mech-class" in str(warning.message) for warning in w)
        finally:
            if mech_backup is not None:
                sys.modules["mech_class"] = mech_backup
            if mech_api_backup is not None:
                sys.modules["mech_class.api"] = mech_api_backup

    def test_dsb_nuclease_scores_zero(self):
        """P(DSB_NUCLEASE)=1.0 -> base=0.0 -> S_DSB=0.0 (mech-class v0.5.1)."""
        mock_root, mock_api = _make_mech_class_mock(p_dsb=1.0, p_dsb_free=0.0)
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("Q99ZW2")
        assert s is not None
        assert abs(s - 0.0) < 1e-6

    def test_dsb_free_scores_high(self):
        """P(DSB_NUCLEASE)=0.0 -> base=1.0 -> S_DSB=1.0 (without composite)."""
        mock_root, mock_api = _make_mech_class_mock(p_dsb=0.0, p_dsb_free=1.0, composite=False)
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A2X3M8B0")
        assert s is not None
        assert abs(s - 1.0) < 1e-6

    def test_composite_flag_adds_bonus(self):
        """composite=True adds +0.1 to base; capped at 1.0."""
        mock_root, mock_api = _make_mech_class_mock(p_dsb=0.05, p_dsb_free=0.95, composite=True)
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A2X3M8B0")
        # base = 1.0 - 0.05 = 0.95, +0.1 = 1.0 (capped)
        assert s is not None
        assert abs(s - 1.0) < 1e-6

    def test_composite_flag_bonus_not_exceeding_one(self):
        """min(1.0, base + 0.1) must never exceed 1.0."""
        mock_root, mock_api = _make_mech_class_mock(p_dsb=0.0, p_dsb_free=1.0, composite=True)
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A2X3M8B0")
        assert s is not None
        assert s <= 1.0

    def test_intermediate_probability(self):
        """P(DSB_NUCLEASE)=0.3 -> base=0.7, no composite -> S_DSB~0.7."""
        mock_root, mock_api = _make_mech_class_mock(p_dsb=0.3, p_dsb_free=0.7, composite=False)
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A7X9Y9Z0")
        assert s is not None
        assert abs(s - 0.7) < 1e-6

    def test_is110_tier_a_gate_override_returns_one(self):
        """IS110 Tier-A gate (mech-class v0.5.2): tier_a_gate_override=True -> S_DSB=1.0.

        Simulates a novel IS110 protein (e.g. IS621) whose OOD ML features would
        have given P(DSB_NUCLEASE)=0.567-0.703.  The v0.5.2 hard gate detects
        PF01548 and PF02371 and sets tier_a_gate_override=True, bypassing the ML score.
        pen-score must return 1.0 regardless of the raw probabilities.
        """
        mock_root, mock_api = _make_mech_class_mock(
            p_dsb=0.567,  # buggy OOD ML output without gate
            p_dsb_free=0.433,
            composite=True,
            tier_a_gate_override=True,
        )
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A7C9VKZ0")  # IS621 accession
        assert s is not None
        assert abs(s - 1.0) < 1e-6, f"Expected S_DSB=1.0 for IS110 gate override, got {s}"

    def test_is110_gate_override_without_composite_still_returns_one(self):
        """tier_a_gate_override=True short-circuits before composite check -> always 1.0."""
        mock_root, mock_api = _make_mech_class_mock(
            p_dsb=0.703,
            p_dsb_free=0.297,
            composite=False,  # gate fires before composite is checked
            tier_a_gate_override=True,
        )
        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            s = dsb.score("A0A7C9VKZ0")
        assert s is not None
        assert abs(s - 1.0) < 1e-6

    def test_predictor_exception_returns_none(self):
        """If predict_from_sequence raises, S_DSB returns None with warning."""
        mock_root, mock_api = _make_mech_class_mock()
        # Override predict_from_sequence to raise
        mock_api.Predictor.load().predict_from_sequence.side_effect = RuntimeError(  # type: ignore[attr-defined]
            "model load failed"
        )

        with patch.dict("sys.modules", {"mech_class": mock_root, "mech_class.api": mock_api}):
            from pen_score.axes import dsb

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = dsb.score("Q99ZW2")

        assert result is None
        assert any("S_DSB" in str(warning.message) for warning in w)


# ---------------------------------------------------------------------------
# S_Deliv  (deliv.py) - UniProt fallback path
# ---------------------------------------------------------------------------


class TestSDelivUniProtPath:
    """Tests for S_Deliv when total_aa=None and UniProt must be queried."""

    def test_uniprot_returns_length_correctly(self):
        """When total_aa=None, fetch_sequence_length is called and sigmoid applied."""
        with patch("pen_score.axes.deliv.fetch_sequence_length", return_value=300):
            from pen_score.axes import deliv

            # Clear the _LENGTH_MAP sentinel for this accession
            deliv._LENGTH_MAP["FAKE_ACC"] = 0  # sentinel triggers fetch
            s = deliv.score("FAKE_ACC", total_aa=None)
        assert s is not None
        # 300 aa -> sigmoid(0.005*(300-900)) = sigmoid(-3) ~ 0.953
        assert s > 0.9

    def test_uniprot_failure_returns_none(self):
        """If UniProt fetch raises, S_Deliv returns None with warning."""
        with patch(
            "pen_score.axes.deliv.fetch_sequence_length",
            side_effect=ConnectionError("network error"),
        ):
            from pen_score.axes import deliv

            deliv._LENGTH_MAP["MISSING_ACC"] = 0
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = deliv.score("MISSING_ACC", total_aa=None)
        assert result is None
        assert any("S_Deliv" in str(warning.message) for warning in w)

    def test_total_aa_provided_bypasses_uniprot(self):
        """When total_aa is given directly, no network call is made."""
        with patch(
            "pen_score.axes.deliv.fetch_sequence_length",
            side_effect=AssertionError("should not call UniProt"),
        ):
            from pen_score.axes import deliv

            s = deliv.score("WHATEVER", total_aa=500)
        assert s is not None
        assert 0.0 < s < 1.0

    def test_accession_not_in_map_falls_through_to_uniprot(self):
        """Unknown accession not in _LENGTH_MAP -> fetches from UniProt."""
        from pen_score.axes import deliv

        # Ensure the key is absent from the map
        deliv._LENGTH_MAP.pop("COMPLETELY_NEW_ACC", None)
        with patch("pen_score.axes.deliv.fetch_sequence_length", return_value=600):
            s = deliv.score("COMPLETELY_NEW_ACC", total_aa=None)
        assert s is not None
        # 600 aa -> sigmoid(0.005*(600-900)) = sigmoid(-1.5) ~ 0.82
        assert 0.75 < s < 0.90


# ---------------------------------------------------------------------------
# d7_homolumo.py - ImportError branch (xtb not installed in CI)
# ---------------------------------------------------------------------------


class TestD7HomoLumo:
    """Tests for d7_homolumo - verifies graceful ImportError handling."""

    def test_xtb_not_installed_returns_none(self, tmp_path):
        """If xtb-python is absent, compute_homolumo_gap returns None + warning."""
        xtb_backup = sys.modules.pop("xtb", None)
        try:
            from pen_score.axes.d7_homolumo import compute_homolumo_gap

            dummy_pdb = tmp_path / "fake.pdb"
            dummy_pdb.write_text("")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = compute_homolumo_gap(dummy_pdb, [("A", 10, "ASP")])
            assert result is None
            assert any("xtb" in str(warning.message).lower() for warning in w)
        finally:
            if xtb_backup is not None:
                sys.modules["xtb"] = xtb_backup
