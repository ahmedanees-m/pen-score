"""Tests for pen_score.axes.d7_homolumo - ImportError path only (xtb block is pragma: no cover)."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from unittest.mock import patch


class TestComputeHomolumoGapImportError:
    """Test the ImportError fallback - no xtb-python required."""

    def test_returns_none_when_xtb_not_installed(self, tmp_path):
        """compute_homolumo_gap returns None with a warning if xtb-python is absent."""
        from pen_score.axes.d7_homolumo import compute_homolumo_gap

        # Temporarily hide xtb from imports
        with patch.dict(sys.modules, {"xtb": None, "xtb.interface": None, "xtb.libxtb": None}):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = compute_homolumo_gap(
                    pdb_path=tmp_path / "dummy.pdb",
                    catalytic_residues=[("A", 100, "ASP")],
                )

        assert result is None
        warning_msgs = [str(warning.message) for warning in w]
        assert any("xtb-python" in msg or "HOMO-LUMO" in msg for msg in warning_msgs)

    def test_warning_category_is_user_warning(self, tmp_path):
        """ImportError path should emit a UserWarning (default for warnings.warn)."""
        from pen_score.axes.d7_homolumo import compute_homolumo_gap

        with patch.dict(sys.modules, {"xtb": None, "xtb.interface": None, "xtb.libxtb": None}):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                compute_homolumo_gap(
                    pdb_path=tmp_path / "dummy.pdb",
                    catalytic_residues=[("A", 100, "ASP")],
                )

        assert len(w) >= 1
        assert any(issubclass(warning.category, (UserWarning, ImportWarning)) for warning in w)

    def test_returns_none_type(self, tmp_path):
        from pen_score.axes.d7_homolumo import compute_homolumo_gap

        with patch.dict(sys.modules, {"xtb": None, "xtb.interface": None, "xtb.libxtb": None}):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                result = compute_homolumo_gap(
                    pdb_path=Path("/nonexistent/path.pdb"),
                    catalytic_residues=[("A", 1, "MET"), ("B", 50, "GLU")],
                )

        assert result is None

    def test_module_docstring_present(self):
        """Module has docstring documenting the Docker/VM requirement."""
        import pen_score.axes.d7_homolumo as mod

        assert mod.__doc__ is not None
        assert "xtb" in mod.__doc__.lower() or "biophysics" in mod.__doc__.lower()
