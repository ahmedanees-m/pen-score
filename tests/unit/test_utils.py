"""Unit tests for pen_score.utils modules.

Covers pen_score.utils.uniprot and pen_score.utils.pubmed using
unittest.mock.patch on requests.get - no real network calls.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resp(json_data: dict, status: int = 200) -> MagicMock:
    """Build a minimal mock requests.Response."""
    m = MagicMock()
    m.status_code = status
    m.json.return_value = json_data
    m.raise_for_status.return_value = None
    return m


def _uniprot_entry(accession: str = "Q99ZW2", length: int = 1368, sequence: str = "MAAA") -> dict:
    """Return a minimal UniProt REST JSON entry."""
    return {
        "primaryAccession": accession,
        "sequence": {"value": sequence, "length": length},
        "uniProtKBCrossReferences": [
            {"database": "Pfam", "id": "PF01548"},
            {"database": "Pfam", "id": "PF02371"},
            {"database": "PDB", "id": "4UN3"},  # should be ignored
        ],
    }


# ---------------------------------------------------------------------------
# utils.uniprot
# ---------------------------------------------------------------------------


class TestFetchEntry:
    """Tests for pen_score.utils.uniprot.fetch_entry."""

    def test_returns_dict_on_200(self):
        from pen_score.utils.uniprot import fetch_entry

        entry = _uniprot_entry()
        with patch("requests.get", return_value=_resp(entry)):
            result = fetch_entry("Q99ZW2")
        assert isinstance(result, dict)
        assert result["primaryAccession"] == "Q99ZW2"

    def test_raises_on_404(self):
        from pen_score.utils.uniprot import fetch_entry

        bad_resp = MagicMock()
        bad_resp.status_code = 404
        bad_resp.raise_for_status.side_effect = Exception("404 Not Found")
        with patch("requests.get", return_value=bad_resp):
            with pytest.raises(Exception, match="404"):
                fetch_entry("NONEXISTENT")

    def test_429_retried(self):
        """429 rate-limit responses are retried up to 4 times."""
        from pen_score.utils.uniprot import fetch_entry

        entry = _uniprot_entry()
        # First two calls return 429; third returns 200
        side_effects = [
            _resp({}, status=429),
            _resp({}, status=429),
            _resp(entry, status=200),
        ]
        with patch("requests.get", side_effect=side_effects):
            with patch("time.sleep"):  # skip actual sleeps in tests
                result = fetch_entry("Q99ZW2")
        assert result["primaryAccession"] == "Q99ZW2"


class TestFetchSequence:
    def test_returns_sequence_string(self):
        from pen_score.utils.uniprot import fetch_sequence

        entry = _uniprot_entry(sequence="MDIIQRTPKIQVYSRHPPEVGSSREQ")
        with patch("requests.get", return_value=_resp(entry)):
            seq = fetch_sequence("Q99ZW2")
        assert seq == "MDIIQRTPKIQVYSRHPPEVGSSREQ"

    def test_empty_sequence_returns_empty_string(self):
        from pen_score.utils.uniprot import fetch_sequence

        entry = _uniprot_entry(sequence="")
        with patch("requests.get", return_value=_resp(entry)):
            seq = fetch_sequence("Q99ZW2")
        assert seq == ""


class TestFetchSequenceLength:
    def test_returns_int_length(self):
        from pen_score.utils.uniprot import fetch_sequence_length

        entry = _uniprot_entry(length=1368)
        with patch("requests.get", return_value=_resp(entry)):
            length = fetch_sequence_length("Q99ZW2")
        assert length == 1368
        assert isinstance(length, int)

    def test_small_editor_length(self):
        from pen_score.utils.uniprot import fetch_sequence_length

        entry = _uniprot_entry(length=300)
        with patch("requests.get", return_value=_resp(entry)):
            length = fetch_sequence_length("SMALL_ACC")
        assert length == 300


class TestFetchPfamHits:
    def test_returns_only_pfam_ids(self):
        """Only cross-references with database == 'Pfam' should be returned."""
        from pen_score.utils.uniprot import fetch_pfam_hits

        entry = _uniprot_entry()  # has PF01548, PF02371, and one PDB ref
        with patch("requests.get", return_value=_resp(entry)):
            hits = fetch_pfam_hits("Q99ZW2")
        assert "PF01548" in hits
        assert "PF02371" in hits
        assert "4UN3" not in hits  # PDB ref must be excluded
        assert len(hits) == 2

    def test_no_pfam_hits_returns_empty_list(self):
        from pen_score.utils.uniprot import fetch_pfam_hits

        entry = {
            "primaryAccession": "NOPFAM",
            "sequence": {"value": "M", "length": 1},
            "uniProtKBCrossReferences": [],
        }
        with patch("requests.get", return_value=_resp(entry)):
            hits = fetch_pfam_hits("NOPFAM")
        assert hits == []


class TestSearchByPfam:
    def _search_resp(self, n: int, has_link: bool = False) -> MagicMock:
        results = [{"accession": f"ACC{i}"} for i in range(n)]
        data = {"results": results}
        if has_link:
            data["link"] = "https://rest.uniprot.org/uniprotkb/search?cursor=next"
        m = MagicMock()
        m.status_code = 200
        m.json.return_value = data
        m.raise_for_status.return_value = None
        return m

    def test_returns_list_of_entries(self):
        from pen_score.utils.uniprot import search_by_pfam

        with patch("requests.get", return_value=self._search_resp(5)):
            results = search_by_pfam("PF01548", max_results=10)
        assert len(results) == 5
        assert results[0]["accession"] == "ACC0"

    def test_pagination_followed(self):
        """When a link is returned, a second GET is issued."""
        from pen_score.utils.uniprot import search_by_pfam

        page1 = self._search_resp(3, has_link=True)
        page2 = self._search_resp(2, has_link=False)
        with patch("requests.get", side_effect=[page1, page2]):
            results = search_by_pfam("PF01548", max_results=100)
        assert len(results) == 5

    def test_max_results_respected(self):
        """Results are truncated to max_results."""
        from pen_score.utils.uniprot import search_by_pfam

        with patch("requests.get", return_value=self._search_resp(20)):
            results = search_by_pfam("PF01548", max_results=5)
        assert len(results) == 5

    def test_reviewed_filter_in_query(self):
        """reviewed=True must include 'reviewed:true' in the query."""
        from pen_score.utils.uniprot import search_by_pfam

        with patch("requests.get", return_value=self._search_resp(0)) as mock_get:
            search_by_pfam("PF01548", reviewed=True)
        call_kwargs = mock_get.call_args
        params = call_kwargs[1].get("params") or call_kwargs[0][1]
        assert "reviewed:true" in params.get("query", "")


# ---------------------------------------------------------------------------
# utils.pubmed
# ---------------------------------------------------------------------------


class TestFetchPubmedCount:
    def _pubmed_resp(self, count: int) -> MagicMock:
        m = MagicMock()
        m.status_code = 200
        m.json.return_value = {"esearchresult": {"count": str(count)}}
        m.raise_for_status.return_value = None
        return m

    def test_returns_integer_count(self):
        from pen_score.utils.pubmed import fetch_pubmed_count

        with patch("requests.get", return_value=self._pubmed_resp(987)):
            count = fetch_pubmed_count(["SpCas9", "Cas9"])
        assert count == 987
        assert isinstance(count, int)

    def test_zero_count(self):
        from pen_score.utils.pubmed import fetch_pubmed_count

        with patch("requests.get", return_value=self._pubmed_resp(0)):
            count = fetch_pubmed_count(["novel_editor_no_papers"])
        assert count == 0

    def test_single_term(self):
        from pen_score.utils.pubmed import fetch_pubmed_count

        with patch("requests.get", return_value=self._pubmed_resp(42)):
            count = fetch_pubmed_count(["IS621"])
        assert count == 42

    def test_query_contains_editor_term(self):
        """The HTTP request query string must contain the editor name."""
        from pen_score.utils.pubmed import fetch_pubmed_count

        with patch("requests.get", return_value=self._pubmed_resp(10)) as mock_get:
            fetch_pubmed_count(["IS621"])
        call_kwargs = mock_get.call_args
        params = call_kwargs[1].get("params") or call_kwargs[0][1]
        assert "IS621" in params.get("term", "")

    def test_extra_terms_included(self):
        """custom extra_terms are included in the query."""
        from pen_score.utils.pubmed import fetch_pubmed_count

        with patch("requests.get", return_value=self._pubmed_resp(10)) as mock_get:
            fetch_pubmed_count(["IS621"], extra_terms="AAV delivery")
        params = mock_get.call_args[1].get("params") or mock_get.call_args[0][1]
        assert "AAV delivery" in params.get("term", "")

    def test_429_triggers_retry(self):
        """Two 429 responses followed by a 200 should succeed."""
        from pen_score.utils.pubmed import fetch_pubmed_count

        responses = [
            MagicMock(status_code=429, raise_for_status=MagicMock()),
            MagicMock(status_code=429, raise_for_status=MagicMock()),
            self._pubmed_resp(55),
        ]
        with patch("requests.get", side_effect=responses):
            with patch("time.sleep"):
                count = fetch_pubmed_count(["SpCas9"])
        assert count == 55

    def test_http_error_raises(self):
        """Non-retryable HTTP error (e.g., 500) should propagate."""
        from pen_score.utils.pubmed import fetch_pubmed_count

        bad = MagicMock()
        bad.status_code = 500
        bad.raise_for_status.side_effect = Exception("Internal Server Error")
        with patch("requests.get", return_value=bad):
            with pytest.raises(Exception, match="Internal Server Error"):
                fetch_pubmed_count(["SpCas9"])


class TestBatchFetchCounts:
    def _pubmed_resp(self, count: int) -> MagicMock:
        m = MagicMock()
        m.status_code = 200
        m.json.return_value = {"esearchresult": {"count": str(count)}}
        m.raise_for_status.return_value = None
        return m

    def test_returns_dict_with_all_editors(self):
        """batch_fetch_counts returns one entry per input editor."""
        from pen_score.utils.pubmed import batch_fetch_counts

        editor_map = {
            "SpCas9": ["SpCas9", "Cas9"],
            "IS621": ["IS621", "bridge recombinase"],
        }
        responses = [self._pubmed_resp(1000), self._pubmed_resp(50)]
        with patch("requests.get", side_effect=responses):
            with patch("time.sleep"):  # skip throttle sleep
                counts = batch_fetch_counts(editor_map)
        assert set(counts.keys()) == {"SpCas9", "IS621"}
        assert counts["SpCas9"] == 1000
        assert counts["IS621"] == 50

    def test_empty_input_returns_empty_dict(self):
        from pen_score.utils.pubmed import batch_fetch_counts

        with patch("requests.get") as mock_get:
            counts = batch_fetch_counts({})
        assert counts == {}
        mock_get.assert_not_called()

    def test_sleep_called_between_requests(self):
        """batch_fetch_counts sleeps between API calls to respect rate limit."""
        from pen_score.utils.pubmed import batch_fetch_counts

        editor_map = {"A": ["termA"], "B": ["termB"]}
        responses = [self._pubmed_resp(10), self._pubmed_resp(20)]
        with patch("requests.get", side_effect=responses):
            with patch("time.sleep") as mock_sleep:
                batch_fetch_counts(editor_map, sleep_between=0.34)
        # sleep should be called once (between the two requests)
        assert mock_sleep.call_count >= 1
