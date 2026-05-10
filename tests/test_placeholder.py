"""Package smoke test - verifies importability and correct version family."""

import re

import pen_score


def test_version():
    """Verify package has a version attribute in the 0.1.x family."""
    assert hasattr(pen_score, "__version__")
    v = pen_score.__version__
    # Matches 0.1.0, 0.1.0.post1, 0.1.1.dev0+g..., etc.
    assert re.match(r"^0\.1\.", v), (
        f"Expected version 0.1.x, got {v!r}. Run: pip3 install . after tagging v0.1.0 to fix."
    )


def test_scorer_importable():
    """Scorer class must be importable without optional dependencies."""
    from pen_score.api import Scorer

    scorer = Scorer.load()
    assert scorer is not None


def test_cli_importable():
    """CLI entry-point must be importable."""
    from pen_score.cli import main

    assert callable(main)
