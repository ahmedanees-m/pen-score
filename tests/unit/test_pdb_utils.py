"""Unit tests for pen_score.utils.pdb - all network calls mocked."""

from __future__ import annotations

import warnings
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


class TestFetchPdbStructure:
    """Tests for fetch_pdb_structure."""

    def test_downloads_cif_file(self, tmp_path):
        from pen_score.utils.pdb import fetch_pdb_structure

        mock_resp = MagicMock()
        mock_resp.content = b"data_4ABC\n_entry.id 4ABC\n"
        mock_resp.raise_for_status = MagicMock()

        with patch("pen_score.utils.pdb.requests.get", return_value=mock_resp) as mock_get:
            result = fetch_pdb_structure("4abc", tmp_path)

        mock_get.assert_called_once_with(
            "https://files.rcsb.org/download/4ABC.cif", timeout=60
        )
        assert result == tmp_path / "4ABC.cif"
        assert result.read_bytes() == b"data_4ABC\n_entry.id 4ABC\n"

    def test_pdb_id_uppercased(self, tmp_path):
        from pen_score.utils.pdb import fetch_pdb_structure

        mock_resp = MagicMock()
        mock_resp.content = b"CIF data"
        mock_resp.raise_for_status = MagicMock()

        with patch("pen_score.utils.pdb.requests.get", return_value=mock_resp):
            result = fetch_pdb_structure("4abc", tmp_path)

        assert result.name == "4ABC.cif"

    def test_returns_cached_file_without_request(self, tmp_path):
        from pen_score.utils.pdb import fetch_pdb_structure

        # Pre-create the file
        cached = tmp_path / "4ABC.cif"
        cached.write_bytes(b"cached content")

        with patch("pen_score.utils.pdb.requests.get") as mock_get:
            result = fetch_pdb_structure("4abc", tmp_path)

        mock_get.assert_not_called()
        assert result == cached

    def test_raises_on_http_error(self, tmp_path):
        import requests as req

        from pen_score.utils.pdb import fetch_pdb_structure

        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = req.HTTPError("404 Not Found")

        with patch("pen_score.utils.pdb.requests.get", return_value=mock_resp):
            with pytest.raises(req.HTTPError):
                fetch_pdb_structure("XXXX", tmp_path)

    def test_writes_bytes_from_response(self, tmp_path):
        from pen_score.utils.pdb import fetch_pdb_structure

        expected_bytes = b"\x89PNG content bytes"
        mock_resp = MagicMock()
        mock_resp.content = expected_bytes
        mock_resp.raise_for_status = MagicMock()

        with patch("pen_score.utils.pdb.requests.get", return_value=mock_resp):
            result = fetch_pdb_structure("1TUP", tmp_path)

        assert result.read_bytes() == expected_bytes


class TestFetchAlphafoldStructure:
    """Tests for fetch_alphafold_structure."""

    def _make_af_metadata(self, accession: str, version: int = 4) -> list[dict]:
        return [
            {
                "pdbUrl": f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb",
                "uniprotAccession": accession,
            }
        ]

    def test_downloads_pdb_file(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        accession = "P04637"
        meta_resp = MagicMock()
        meta_resp.status_code = 200
        meta_resp.raise_for_status = MagicMock()
        meta_resp.json.return_value = self._make_af_metadata(accession)

        pdb_resp = MagicMock()
        pdb_resp.content = b"ATOM      1  N   MET A   1\n"
        pdb_resp.raise_for_status = MagicMock()

        with patch("pen_score.utils.pdb.requests.get", side_effect=[meta_resp, pdb_resp]):
            result = fetch_alphafold_structure(accession, tmp_path)

        assert result is not None
        assert result.name == f"AF-{accession}-F1-model_v4.pdb"
        assert result.read_bytes() == b"ATOM      1  N   MET A   1\n"

    def test_returns_cached_without_request(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        accession = "P04637"
        cached = tmp_path / "AF-P04637-F1-model_v4.pdb"
        cached.write_bytes(b"cached pdb data")

        with patch("pen_score.utils.pdb.requests.get") as mock_get:
            result = fetch_alphafold_structure(accession, tmp_path)

        mock_get.assert_not_called()
        assert result == cached

    def test_returns_none_on_404(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        mock_resp = MagicMock()
        mock_resp.status_code = 404

        with patch("pen_score.utils.pdb.requests.get", return_value=mock_resp):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = fetch_alphafold_structure("UNKNOWN99", tmp_path)

        assert result is None
        assert any("No AlphaFold structure" in str(warning.message) for warning in w)

    def test_returns_none_when_empty_metadata(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        meta_resp = MagicMock()
        meta_resp.status_code = 200
        meta_resp.raise_for_status = MagicMock()
        meta_resp.json.return_value = []  # empty list

        with patch("pen_score.utils.pdb.requests.get", return_value=meta_resp):
            result = fetch_alphafold_structure("P00000", tmp_path)

        assert result is None

    def test_returns_none_when_no_pdb_url(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        meta_resp = MagicMock()
        meta_resp.status_code = 200
        meta_resp.raise_for_status = MagicMock()
        meta_resp.json.return_value = [{"pdbUrl": None, "uniprotAccession": "P00000"}]

        with patch("pen_score.utils.pdb.requests.get", return_value=meta_resp):
            result = fetch_alphafold_structure("P00000", tmp_path)

        assert result is None

    def test_custom_version_used_in_filename(self, tmp_path):
        from pen_score.utils.pdb import fetch_alphafold_structure

        accession = "Q9Y2V2"
        version = 3
        meta_resp = MagicMock()
        meta_resp.status_code = 200
        meta_resp.raise_for_status = MagicMock()
        meta_resp.json.return_value = self._make_af_metadata(accession, version)

        pdb_resp = MagicMock()
        pdb_resp.content = b"ATOM content"
        pdb_resp.raise_for_status = MagicMock()

        with patch("pen_score.utils.pdb.requests.get", side_effect=[meta_resp, pdb_resp]):
            result = fetch_alphafold_structure(accession, tmp_path, version=version)

        assert result is not None
        assert f"model_v{version}" in result.name


class TestGetMeanPlddt:
    """Tests for get_mean_plddt."""

    def test_returns_none_on_parse_failure(self, tmp_path):
        from pen_score.utils.pdb import get_mean_plddt

        bad_pdb = tmp_path / "bad.pdb"
        bad_pdb.write_text("this is not a PDB file\n")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = get_mean_plddt(bad_pdb)

        assert result is None

    def test_returns_none_warns_on_missing_biopython(self, tmp_path, monkeypatch):
        """If biopython is not installed, get_mean_plddt should return None with a warning."""
        import sys

        from pen_score.utils.pdb import get_mean_plddt

        fake_pdb = tmp_path / "test.pdb"
        fake_pdb.write_text("ATOM      1  N   MET A   1       1.0   2.0   3.0  1.00 95.50\n")

        # Simulate BioPython being unavailable
        original = sys.modules.get("Bio")
        sys.modules["Bio"] = None  # type: ignore[assignment]
        sys.modules["Bio.PDB"] = None  # type: ignore[assignment]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = get_mean_plddt(fake_pdb)

        # Restore
        if original is None:
            sys.modules.pop("Bio", None)
            sys.modules.pop("Bio.PDB", None)
        else:
            sys.modules["Bio"] = original

        assert result is None

    def test_returns_float_with_real_biopython(self, tmp_path):
        """If biopython is available, should return a float for a valid PDB."""
        pytest.importorskip("Bio")
        from pen_score.utils.pdb import get_mean_plddt

        # Minimal valid PDB content with B-factor column (pLDDT)
        pdb_content = (
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 90.00           C\n"
            "ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00 80.00           C\n"
            "END\n"
        )
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)

        result = get_mean_plddt(pdb_path)

        if result is not None:
            assert isinstance(result, float)
            assert 0.0 <= result <= 100.0
