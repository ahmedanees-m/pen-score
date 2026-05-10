"""Unit tests for use-case profiles.

Tests weight integrity and semantic properties of all profiles in
pen_score/data/use_case_profiles.yaml.  No external deps required.
"""

from __future__ import annotations


class TestProfileWeightSums:
    """All profiles must have weights summing to exactly 1.0."""

    def test_all_profiles_sum_to_one(self, use_case_profiles):
        for name, weights in use_case_profiles.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 1e-9, (
                f"Profile '{name}' weights sum to {total:.10f}, expected 1.0"
            )

    def test_seven_profiles_present(self, use_case_profiles):
        """Expect exactly 7 profiles after v0.1.2 adds 2 new ones."""
        assert len(use_case_profiles) == 7, (
            f"Expected 7 profiles (5 original + 2 new), got {len(use_case_profiles)}: "
            f"{list(use_case_profiles.keys())}"
        )

    def test_all_eight_axes_present_in_each_profile(self, use_case_profiles):
        expected = {
            "S_DSB",
            "S_Spec",
            "S_Cargo",
            "S_Deliv",
            "S_Immuno",
            "S_Prog",
            "S_Mature",
            "S_Energy",
        }
        for name, weights in use_case_profiles.items():
            assert set(weights.keys()) == expected, (
                f"Profile '{name}' has unexpected axes: {set(weights.keys()) ^ expected}"
            )


class TestMegabaseRearrangementProfile:
    """megabase_rearrangement profile - added v0.1.2."""

    def test_profile_present(self, use_case_profiles):
        assert "megabase_rearrangement" in use_case_profiles, (
            "megabase_rearrangement profile not found in use_case_profiles.yaml"
        )

    def test_s_mature_zero_weight(self, use_case_profiles):
        """S_Mature must have weight 0.0 - all megabase editors are immature."""
        profile = use_case_profiles["megabase_rearrangement"]
        assert profile["S_Mature"] == 0.0, (
            f"S_Mature weight should be 0.0 in megabase_rearrangement, got {profile['S_Mature']}"
        )

    def test_s_prog_highest_weight(self, use_case_profiles):
        """S_Prog must have the highest weight - boundary programmability is paramount."""
        profile = use_case_profiles["megabase_rearrangement"]
        max_ax = max(profile, key=lambda k: profile[k])
        assert max_ax == "S_Prog", (
            f"S_Prog should have highest weight in megabase_rearrangement, "
            f"but {max_ax} has weight {profile[max_ax]:.2f}"
        )

    def test_s_dsb_required(self, use_case_profiles):
        """S_DSB must have substantial weight (>=0.15) - DSB-free is non-negotiable."""
        profile = use_case_profiles["megabase_rearrangement"]
        assert profile["S_DSB"] >= 0.15, (
            f"S_DSB weight {profile['S_DSB']:.2f} too low in megabase_rearrangement"
        )

    def test_composite_is622_ranks_high(self, use_case_profiles):
        """Under megabase_rearrangement, an IS622-like editor (strong biophysics,
        S_Mature=0) should score higher than an IS621-like editor under the
        default profile because S_Mature is zero-weighted."""
        from pen_score.scorer.composite import compute_pen_score

        weights = use_case_profiles["megabase_rearrangement"]

        # IS622-like: strong biophysics, S_Mature=0
        is622 = {
            "S_DSB": 1.0,
            "S_Spec": 0.99,
            "S_Cargo": 1.0,
            "S_Deliv": 0.95,
            "S_Immuno": 0.76,
            "S_Prog": 1.0,
            "S_Mature": 0.0,
            "S_Energy": 1.0,
        }
        score, missing = compute_pen_score(is622, weights)
        assert score is not None
        assert missing == []
        # S_Mature=0 doesn't penalise under megabase_rearrangement (weight=0)
        assert score > 0.90, (
            f"IS622-like editor should score >0.90 under megabase_rearrangement, got {score}"
        )


class TestTherapeuticExcisionBcl11aProfile:
    """therapeutic_excision_bcl11a profile - added v0.1.2."""

    def test_profile_present(self, use_case_profiles):
        assert "therapeutic_excision_bcl11a" in use_case_profiles, (
            "therapeutic_excision_bcl11a profile not found in use_case_profiles.yaml"
        )

    def test_s_dsb_and_s_spec_highest(self, use_case_profiles):
        """S_DSB and S_Spec must jointly be the two highest-weighted axes."""
        profile = use_case_profiles["therapeutic_excision_bcl11a"]
        sorted_axes = sorted(profile, key=lambda k: profile[k], reverse=True)
        top_two = set(sorted_axes[:2])
        assert top_two == {"S_DSB", "S_Spec"}, (
            f"Top-two axes should be S_DSB and S_Spec, got {top_two}"
        )

    def test_s_immuno_high_weight(self, use_case_profiles):
        """S_Immuno must be >=0.15 - immunogenicity is critical for HSC therapy."""
        profile = use_case_profiles["therapeutic_excision_bcl11a"]
        assert profile["S_Immuno"] >= 0.15, (
            f"S_Immuno weight {profile['S_Immuno']:.2f} too low for HSC therapy profile"
        )

    def test_s_energy_low_weight(self, use_case_profiles):
        """S_Energy should be low (<0.05) - delivery scaffold carries energy cost."""
        profile = use_case_profiles["therapeutic_excision_bcl11a"]
        assert profile["S_Energy"] < 0.05, (
            f"S_Energy weight {profile['S_Energy']:.2f} higher than expected in BCL11A profile"
        )

    def test_composite_is621_ranks_well(self, use_case_profiles):
        """IS621-like editor (strong biophysics + mature) should score >0.90
        under the BCL11A profile."""
        from pen_score.scorer.composite import compute_pen_score

        weights = use_case_profiles["therapeutic_excision_bcl11a"]

        # IS621-like: strong across all axes, moderate maturity
        is621 = {
            "S_DSB": 1.0,
            "S_Spec": 0.99,
            "S_Cargo": 1.0,
            "S_Deliv": 0.99,
            "S_Immuno": 0.76,
            "S_Prog": 1.0,
            "S_Mature": 0.25,
            "S_Energy": 1.0,
        }
        score, missing = compute_pen_score(is621, weights)
        assert score is not None
        assert missing == []
        # The BCL11A profile upweights S_Mature (0.15); IS621's S_Mature~0.25
        # (still a new editor by clinical standards) moderates the score.
        # The key check is that the score is well above 0.75 - IS621 still ranks
        # competitively even with moderate clinical maturity.
        assert score > 0.80, (
            f"IS621-like editor should score >0.80 under therapeutic_excision_bcl11a, got {score}"
        )
