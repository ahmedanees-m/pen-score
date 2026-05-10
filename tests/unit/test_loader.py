"""Unit tests for editor universe loader and data integrity."""

from __future__ import annotations


class TestEditorUniverse:
    def test_loads_without_error(self, editor_universe):
        assert len(editor_universe) > 0

    def test_n_editors_matches_plan(self, editor_universe):
        # Plan specifies ~30 editors; allow 25-35
        assert 25 <= len(editor_universe) <= 35, f"Expected ~30 editors, got {len(editor_universe)}"

    def test_all_have_canonical_accession(self, editor_universe):
        for ed in editor_universe:
            assert ed.canonical_accession, f"Missing accession for {ed.id}"

    def test_all_have_valid_mechanism_bucket(self, editor_universe):
        valid = {"DSB_NUCLEASE", "DSB_FREE_TRANSEST_RECOMBINASE", "TRANSPOSASE"}
        for ed in editor_universe:
            assert ed.mechanism_bucket in valid, (
                f"{ed.id} has unknown mechanism_bucket: {ed.mechanism_bucket}"
            )

    def test_cargo_capacity_positive(self, editor_universe):
        for ed in editor_universe:
            assert ed.cargo_capacity_bp > 0, f"{ed.id} has non-positive cargo_capacity_bp"

    def test_five_pre_registered_targets(self, editor_universe):
        targets = [ed for ed in editor_universe if ed.pre_registered_target]
        assert len(targets) == 5, f"Expected 5 pre-registered targets, got {len(targets)}"

    def test_pre_registered_target_ids(self, editor_universe):
        target_ids = {ed.id for ed in editor_universe if ed.pre_registered_target}
        assert "evoCAST" in target_ids
        assert "IS621" in target_ids
        assert "SpCas9" in target_ids
        assert "enNlovFz2" in target_ids
        assert "SpuFz1_V4" in target_ids

    def test_bxb1_corrected_accession(self, editor_universe):
        """MECH-CLASS corrected Bxb1 from Q8VVR2 -> Q9B086."""
        bxb1 = next((ed for ed in editor_universe if ed.id == "Bxb1"), None)
        assert bxb1 is not None
        assert bxb1.canonical_accession == "Q9B086", (
            f"Bxb1 accession should be Q9B086, got {bxb1.canonical_accession}"
        )

    def test_tn5_corrected_accession(self, editor_universe):
        """MECH-CLASS corrected Tn5 from P00509 -> Q46731."""
        tn5 = next((ed for ed in editor_universe if ed.id == "Tn5"), None)
        assert tn5 is not None
        assert tn5.canonical_accession == "Q46731", (
            f"Tn5 accession should be Q46731, got {tn5.canonical_accession}"
        )

    def test_is621_composite_architecture(self, editor_universe):
        """IS621 must have composite_architecture=True (IS110 family)."""
        is621 = next((ed for ed in editor_universe if ed.id == "IS621"), None)
        assert is621 is not None
        assert is621.composite_architecture is True

    def test_rna_guided_programmable_consistency(self, editor_universe):
        """Spot-check: Cre and Bxb1 are NOT RNA-guided."""
        for ed in editor_universe:
            if ed.id in ("Cre", "Bxb1", "Lambda_Int", "phiC31"):
                assert ed.rna_guided is False, f"{ed.id} should not be RNA-guided"
            if ed.id in ("SpCas9", "Cas12a", "IS621", "evoCAST"):
                assert ed.rna_guided is True, f"{ed.id} should be RNA-guided"


class TestUseCaseProfiles:
    def test_loads_without_error(self, use_case_profiles):
        assert len(use_case_profiles) > 0

    def test_default_profile_present(self, use_case_profiles):
        assert "human_therapeutic_aav_insertion" in use_case_profiles

    def test_weights_sum_to_one(self, use_case_profiles):
        for profile_name, weights in use_case_profiles.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 1e-6, (
                f"Profile '{profile_name}' weights sum to {total:.6f}, expected 1.0"
            )

    def test_all_axes_present_in_default(self, use_case_profiles):
        expected_axes = {
            "S_DSB",
            "S_Spec",
            "S_Cargo",
            "S_Deliv",
            "S_Immuno",
            "S_Prog",
            "S_Mature",
            "S_Energy",
        }
        default = use_case_profiles.get("human_therapeutic_aav_insertion", {})
        assert set(default.keys()) == expected_axes
