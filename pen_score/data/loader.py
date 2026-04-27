"""Load editor universe and use-case profiles with Pydantic validation."""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel

DATA_DIR = Path(__file__).parent


class EditorEntry(BaseModel):
    id: str
    aliases: list[str] = []
    canonical_accession: str
    organism: str
    organism_taxid: int = 0
    canonical_pdb: str | None = None
    mechanism_bucket: str
    rna_guided: bool
    year_discovered: int
    primary_reference: str
    primary_doi: str | None = None
    cargo_capacity_bp: int
    cargo_capacity_note: str = ""
    composite_architecture: bool = False
    pre_registered_target: bool = False
    parent_editor: str | None = None
    notes: str = ""
    references_used_for_pubmed: list[str] = []
    walker_motif_override: bool | None = None
    # v1.0.7 - PEN-COMPARE v3.2 fields
    intrinsic_cargo_mechanism: bool = False
    cell_based_evidence: bool = False
    cell_based_sources: list[str] = []


class EditorUniverse(BaseModel):
    version: str
    created: str
    curator: str
    n_editors: int
    editors: list[EditorEntry]


def load_editor_universe() -> list[EditorEntry]:
    """Load and validate the curated editor universe from YAML."""
    raw = yaml.safe_load((DATA_DIR / "editor_universe.yaml").read_text(encoding="utf-8"))
    universe = EditorUniverse(**raw)
    return universe.editors


def load_use_case_profiles() -> dict[str, dict[str, float]]:
    """Load axis weight profiles per use case."""
    path = DATA_DIR / "use_case_profiles.yaml"
    if not path.exists():
        return {}
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    return raw.get("profiles", {})


def load_axis_definitions() -> dict:
    """Load score axis mathematical definitions."""
    path = DATA_DIR / "axis_definitions.yaml"
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text(encoding="utf-8"))
