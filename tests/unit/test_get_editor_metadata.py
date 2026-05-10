"""Tests for pen_score.get_editor_metadata() API (pen-score v0.1.3).

These tests verify the two boolean fields required by PEN-COMPARE v3.2:
  - intrinsic_cargo_mechanism  (Gate 3: native cargo vs HDR-template)
  - cell_based_evidence         (TRUE_WRITER tier: mammalian cell >1%)

The critical calibration values that PEN-COMPARE v3.2 depends on:
  ISCro4: intrinsic=True, cell_based=True   (TRUE_WRITER anchor)
  IS621:  intrinsic=True, cell_based=False  (PROBABLE_WRITER anchor)
  SpCas9: intrinsic=False, cell_based=True  (Gate 3 exclusion test)
"""

from __future__ import annotations

import warnings

import pytest

from pen_score import get_editor_metadata
from pen_score.api import EditorMetadata

# ---------------------------------------------------------------------------
# ISCro4 (canonical TRUE_WRITER anchor)
# ---------------------------------------------------------------------------

class TestISCro4:
    def test_canonical_name(self):
        md = get_editor_metadata("ISCro4")
        assert md.canonical_name == "ISCro4"

    def test_uniprot(self):
        md = get_editor_metadata("ISCro4")
        assert md.uniprot == "D2TGM5"

    def test_intrinsic_cargo_mechanism_true(self):
        """ISCro4 is an IS110 bridge recombinase - native cargo insertion."""
        md = get_editor_metadata("ISCro4")
        assert md.intrinsic_cargo_mechanism is True

    def test_cell_based_evidence_true(self):
        """ISCro4 has >6% insertion efficiency in human cells (Pelea 2026 Science)."""
        md = get_editor_metadata("ISCro4")
        assert md.cell_based_evidence is True, (
            "ISCro4 must have cell_based_evidence=True - "
            "this is the TRUE_WRITER anchor for PEN-COMPARE v3.2"
        )

    def test_cell_based_sources_populated(self):
        md = get_editor_metadata("ISCro4")
        assert len(md.cell_based_sources) >= 2
        sources_str = " ".join(md.cell_based_sources)
        assert "Pelea" in sources_str or "adz1884" in sources_str

    def test_is622_in_aliases(self):
        md = get_editor_metadata("ISCro4")
        assert "IS622" in md.aliases


# ---------------------------------------------------------------------------
# IS621 (PROBABLE_WRITER calibration anchor - no cell-based evidence)
# ---------------------------------------------------------------------------

class TestIS621:
    def test_intrinsic_cargo_true(self):
        """IS621 is also an IS110 bridge recombinase."""
        md = get_editor_metadata("IS621")
        assert md.intrinsic_cargo_mechanism is True

    def test_cell_based_evidence_false(self):
        """IS621 has only E. coli / cryo-EM evidence; no robust human-cell data.

        This is the CRITICAL PEN-COMPARE v3.2 keystone: IS621 must be False to
        be classified PROBABLE_WRITER (not TRUE_WRITER), distinguishing it from
        ISCro4 which has >6% human-cell insertion (Pelea 2026 Science).
        """
        md = get_editor_metadata("IS621")
        assert md.cell_based_evidence is False, (
            "IS621 must be False - this is what distinguishes it from ISCro4 "
            "in PEN-COMPARE v3.2 TRUE_WRITER tier"
        )

    def test_cell_based_sources_empty(self):
        md = get_editor_metadata("IS621")
        assert md.cell_based_sources == []


# ---------------------------------------------------------------------------
# SpCas9 (Gate 3 reformulation test: HDR-template-based)
# ---------------------------------------------------------------------------

class TestSpCas9:
    def test_intrinsic_cargo_mechanism_false(self):
        """SpCas9 uses an HDR donor template - not native cargo mechanism.

        If this were True, Gate 3 in PEN-COMPARE v3.2 would fail to exclude
        SpCas9+HDR from the native cargo class.
        """
        md = get_editor_metadata("SpCas9")
        assert md.intrinsic_cargo_mechanism is False, (
            "SpCas9 uses HDR template, not native cargo. "
            "intrinsic_cargo_mechanism must be False for Gate 3 to work."
        )

    def test_cell_based_evidence_true(self):
        """SpCas9 has extensive mammalian cell literature."""
        md = get_editor_metadata("SpCas9")
        assert md.cell_based_evidence is True


# ---------------------------------------------------------------------------
# IS622 alias resolution (backward-compatibility test)
# ---------------------------------------------------------------------------

class TestAliasResolution:
    def test_is622_resolves_to_iscro4(self):
        """IS622 is a deprecated alias that must resolve to ISCro4."""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            md = get_editor_metadata("IS622")
        assert md.canonical_name == "ISCro4"

    def test_is622_emits_deprecation_warning(self):
        """Alias lookup must emit a DeprecationWarning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            get_editor_metadata("IS622")
        assert any(
            issubclass(warning.category, DeprecationWarning) and
            "deprecated" in str(warning.message).lower()
            for warning in w
        ), "Expected DeprecationWarning for deprecated alias IS622"

    def test_is622_alias_returns_same_metadata_as_iscro4(self):
        """Alias and canonical ID must return identical metadata."""
        md_canonical = get_editor_metadata("ISCro4")
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            md_alias = get_editor_metadata("IS622")
        assert md_alias.editor_id == md_canonical.editor_id
        assert md_alias.intrinsic_cargo_mechanism == md_canonical.intrinsic_cargo_mechanism
        assert md_alias.cell_based_evidence == md_canonical.cell_based_evidence


# ---------------------------------------------------------------------------
# Return type and error handling
# ---------------------------------------------------------------------------

class TestAPIContract:
    def test_returns_editor_metadata_instance(self):
        md = get_editor_metadata("ISCro4")
        assert isinstance(md, EditorMetadata)

    def test_metadata_is_frozen(self):
        """EditorMetadata must be immutable."""
        md = get_editor_metadata("ISCro4")
        with pytest.raises((AttributeError, TypeError)):
            md.editor_id = "something_else"  # type: ignore[misc]

    def test_unknown_editor_raises_key_error(self):
        with pytest.raises(KeyError):
            get_editor_metadata("NotARealEditor_v99")

    def test_all_editors_have_required_fields(self):
        """Every editor in universe must have intrinsic_cargo_mechanism + cell_based_evidence."""
        from pen_score.data.loader import load_editor_universe
        for ed in load_editor_universe():
            md = get_editor_metadata(ed.id)
            assert isinstance(md.intrinsic_cargo_mechanism, bool), \
                f"{ed.id}: intrinsic_cargo_mechanism must be bool"
            assert isinstance(md.cell_based_evidence, bool), \
                f"{ed.id}: cell_based_evidence must be bool"
