# pen-stack/biophysics:0.1.0
# Extends pen-stack/structure:0.1.0 with:
#   - xtb-python (GFN2-xTB semiempirical quantum chemistry) for d_7 HOMO-LUMO gap
#   - APBS (Adaptive Poisson-Boltzmann Solver) for electrostatic surfaces
#   - PDB2PQR (protonation state assignment for APBS input)
#   - fpocket (pocket detection for active-site geometry)
#   - rdkit (SMILES <-> 3D conversion utilities)
#
# Build:
#   docker build -t pen-stack/biophysics:0.1.0 -t pen-stack/biophysics:latest .
#
# Run (on VM, with data and repo mounts):
#   docker run --rm \
#       -v ~/pen-stack/data:/data \
#       -v ~/pen-stack/code/repos/pen-score:/pkg \
#       -w /pkg pen-stack/biophysics:0.1.0 \
#       python scripts/00_smoke_import.py

FROM pen-stack/structure:0.1.0

# â”€â”€ System packages â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# fpocket is not in Ubuntu 22.04 apt repos and requires custom compilation;
# it is supplementary (d7 HOMO-LUMO active-site cluster) and does not
# gate any headline claim - install separately on the host if needed.
RUN apt-get update && apt-get install -y --no-install-recommends \
    apbs \
    pdb2pqr \
    && rm -rf /var/lib/apt/lists/*

# â”€â”€ Python packages â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
# xtb-python has no PyPI wheel for Ubuntu 22.04 + Python 3.x from CUDA base;
# install via conda-forge separately when running the d7 HOMO-LUMO step.
# rdkit provides SMILES <-> 3D conversion utilities for cluster extraction.
RUN pip install --no-cache-dir "rdkit==2024.3.1"

# â”€â”€ Smoke test â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
RUN python -c "from rdkit import Chem; print('rdkit:', Chem.__name__, 'OK'); print('biophysics smoke OK')"

CMD ["/bin/bash"]
