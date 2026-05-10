"""Shared pytest fixtures for pen-score tests."""

from __future__ import annotations

from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(scope="session")
def editor_universe():
    """Return the full editor universe list."""
    from pen_score.data.loader import load_editor_universe

    return load_editor_universe()


@pytest.fixture(scope="session")
def use_case_profiles():
    """Return the use-case weight profiles dict."""
    from pen_score.data.loader import load_use_case_profiles

    return load_use_case_profiles()


@pytest.fixture
def mock_axis_scores():
    """Return a minimal axis scores dict for composite function tests (8 axes, v0.1.1)."""
    return {
        "S_DSB": 1.0,
        "S_Spec": 0.7,
        "S_Cargo": 1.0,
        "S_Deliv": 0.95,
        "S_Immuno": 0.8,
        "S_Prog": 1.0,
        "S_Mature": 0.3,
        "S_Energy": 1.0,
    }


@pytest.fixture
def default_weights():
    """Return default use-case weights for human_therapeutic_aav_insertion (8 axes, v0.1.1)."""
    from pen_score.data.loader import load_use_case_profiles

    profiles = load_use_case_profiles()
    return profiles.get(
        "human_therapeutic_aav_insertion",
        {
            "S_DSB": 0.24,
            "S_Spec": 0.14,
            "S_Cargo": 0.19,
            "S_Deliv": 0.19,
            "S_Immuno": 0.09,
            "S_Prog": 0.05,
            "S_Mature": 0.05,
            "S_Energy": 0.05,
        },
    )
