"""PDB structure helpers for active-site extraction and structure fetching."""

from __future__ import annotations

import warnings
from pathlib import Path

import requests

_RCSB_API = "https://data.rcsb.org/rest/v1/core/entry"
_RCSB_FILE = "https://files.rcsb.org/download"
_AF_API = "https://alphafold.ebi.ac.uk/api/prediction"


def fetch_pdb_structure(pdb_id: str, out_dir: Path) -> Path:
    """Download a PDB file from RCSB.

    Parameters
    ----------
    pdb_id:
        4-character PDB ID (case-insensitive).
    out_dir:
        Directory to write the .cif file.

    Returns
    -------
    Path to the downloaded file.
    """
    pdb_id = pdb_id.upper()
    out_path = out_dir / f"{pdb_id}.cif"
    if out_path.exists():
        return out_path
    url = f"{_RCSB_FILE}/{pdb_id}.cif"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)
    return out_path


def fetch_alphafold_structure(accession: str, out_dir: Path, version: int = 4) -> Path | None:
    """Download AlphaFold predicted structure for a UniProt accession.

    Parameters
    ----------
    accession:
        UniProt accession.
    out_dir:
        Directory to write the PDB file.
    version:
        AlphaFold DB version (default: 4).

    Returns
    -------
    Path to the PDB file, or None if unavailable.
    """
    out_path = out_dir / f"AF-{accession}-F1-model_v{version}.pdb"
    if out_path.exists():
        return out_path
    url = f"{_AF_API}/{accession}"
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        warnings.warn(f"No AlphaFold structure for {accession}.", stacklevel=2)
        return None
    resp.raise_for_status()
    entries = resp.json()
    if not entries:
        return None
    pdb_url = entries[0].get("pdbUrl")
    if not pdb_url:
        return None
    pdb_resp = requests.get(pdb_url, timeout=60)
    pdb_resp.raise_for_status()
    out_path.write_bytes(pdb_resp.content)
    return out_path


def get_mean_plddt(pdb_path: Path) -> float | None:
    """Compute mean pLDDT from an AlphaFold PDB file (B-factor column).

    Returns
    -------
    Mean pLDDT (0-100) or None if the file cannot be parsed.
    """
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("s", str(pdb_path))
        b_factors = [
            atom.bfactor
            for model in structure
            for chain in model
            for residue in chain
            for atom in residue
            if atom.bfactor > 0
        ]
        if not b_factors:
            return None
        return float(sum(b_factors) / len(b_factors))
    except Exception as exc:
        warnings.warn(f"pLDDT extraction failed for {pdb_path}: {exc}", stacklevel=2)
        return None
