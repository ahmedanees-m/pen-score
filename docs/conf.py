"""Sphinx configuration for pen-score documentation."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    release = _pkg_version("pen-score")
except PackageNotFoundError:
    release = "0.1.0"

version = ".".join(release.split(".")[:2])
project = "pen-score"
copyright = "2026, Anees Ahmed Mahaboob Ali"
author = "Anees Ahmed Mahaboob Ali"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

html_theme = "furo"
html_title = f"pen-score {release}"
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_extra_path = ["scorecards"]

autodoc_typehints = "description"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

html_theme_options = {
    "source_repository": "https://github.com/ahmedanees-m/pen-score/",
    "source_branch": "main",
    "source_directory": "docs/",
}
