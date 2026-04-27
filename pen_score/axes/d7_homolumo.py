"""d₇ HOMO-LUMO gap - supplementary biophysical descriptor.

NOT a primary PenScore axis.  Computed as supplementary biophysical evidence
corroborating the S_DSB / MECH-CLASS Tier A mechanism distinction.

Predicted ordering (pre-registered before any computation):
    Transesterases (DSB-free) > Nucleases (DSB) in HOMO-LUMO gap

Method: GFN2-xTB semiempirical tight-binding applied to truncated active-site
clusters (~20-30 atoms around catalytic metal centre), extracted from PDB or
AlphaFold structures.

Runtime: ~30 s per active-site cluster on CPU.
Requires [biophysics] extra (xtb-python) running inside pen-stack/biophysics
Docker image on the VM.  NOT installable on a CPU-only laptop.

Run script 17_compute_d7_homolumo.py inside pen-stack/biophysics:0.1.0.
"""

from __future__ import annotations

import warnings
from pathlib import Path


def compute_homolumo_gap(
    pdb_path: Path, catalytic_residues: list[tuple[str, int, str]]
) -> float | None:
    """Compute HOMO-LUMO gap (eV) for a truncated active-site cluster.

    Parameters
    ----------
    pdb_path:
        Path to PDB file (experimental or AlphaFold).
    catalytic_residues:
        List of (chain_id, residue_number, residue_name) tuples defining the
        active-site cluster centre.

    Returns
    -------
    HOMO-LUMO gap in eV, or None if computation fails.
    """
    try:
        import numpy as np
        from xtb.interface import Calculator, Param
        from xtb.libxtb import VERBOSITY_MUTED
    except ImportError:
        warnings.warn(
            "xtb-python is not installed. d₇ HOMO-LUMO computation unavailable. "
            "Run inside pen-stack/biophysics:0.1.0 Docker image on the VM.",
            stacklevel=2,
        )
        return None

    try:  # pragma: no cover
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("target", str(pdb_path))

        # Extract active-site atom coordinates and atomic numbers
        atoms = []
        numbers = []
        _ELEMENT_TO_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "MG": 12, "ZN": 30}
        for model in structure:
            for chain in model:
                for residue in chain:
                    chain_id = chain.id
                    res_num = residue.id[1]
                    if any(
                        chain_id == cr[0] and abs(res_num - cr[1]) <= 5 for cr in catalytic_residues
                    ):
                        for atom in residue:
                            el = atom.element.strip().upper() if atom.element else "C"
                            z = _ELEMENT_TO_Z.get(el)
                            if z is not None:
                                atoms.append(atom.coord)
                                numbers.append(z)

        if len(atoms) < 5:
            warnings.warn("Active-site cluster too small (<5 atoms); skipping.", stacklevel=2)
            return None

        coords = np.array(atoms) * 1.8897259886  # Å -> Bohr
        numbers_arr = np.array(numbers)

        calc = Calculator(Param.GFN2xTB, numbers_arr, coords)
        calc.set_verbosity(VERBOSITY_MUTED)
        res = calc.singlepoint()
        orbital_e = res.get_orbital_energies()  # Hartree
        n_occ = res.get_number_of_electrons() // 2
        homo = orbital_e[n_occ - 1]
        lumo = orbital_e[n_occ]
        gap_ev = (lumo - homo) * 27.2114  # Hartree -> eV
        return round(float(gap_ev), 4)

    except Exception as exc:  # pragma: no cover
        warnings.warn(f"d₇ HOMO-LUMO computation failed: {exc}", stacklevel=2)
        return None
