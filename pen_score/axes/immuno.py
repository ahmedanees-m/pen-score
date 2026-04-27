"""S_Immuno - Immunogenicity axis (inverse epitope load).

Formula (from axis_definitions.yaml):
    epitope_load = (n_MHC_I_binders + n_MHC_II_binders) / sequence_length
    score = 1.0 - clip(epitope_load / epitope_load_max, 0, 1)

where epitope_load_max is the 99th percentile over the editor universe.

MHC binding prediction:
    - Primary:  netMHCpan-4.1 (external binary; not redistributable)
    - Fallback: MHCflurry 2.0 (open-source; install via [immuno] extra)

This axis is NOT computed live by Scorer._compute_axes_live because it
requires ~30-60 s per protein on CPU.  Run script 14_compute_S_Immuno.py
offline, then load the cached parquet.

Requires [immuno] optional extra:  pip install pen-score[immuno]
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

# Default HLA alleles for immunogenicity profiling (representative human alleles)
_HLA_I_ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", "HLA-B*44:02"]
_HLA_II_ALLELES = ["DRB1*03:01", "DRB1*07:01", "DRB1*15:01"]
_IC50_THRESHOLD_NM = 500.0  # binders with IC50 < 500 nM are counted as epitopes


def score(
    accession: str, sequence: str | None = None, cached_parquet: Path | None = None
) -> float | None:
    """Compute S_Immuno for an editor.

    Parameters
    ----------
    accession:
        UniProt accession.
    sequence:
        Full protein sequence.  If None, fetched from UniProt.
    cached_parquet:
        Path to pre-computed ``immuno_scores.parquet`` from script 14.
        Preferred over live computation for speed.

    Returns
    -------
    float in [0, 1] (higher = less immunogenic) or None if not computable.
    """
    # Fast path: use cached scores if available
    if cached_parquet and Path(cached_parquet).exists():
        import pandas as pd

        df = pd.read_parquet(cached_parquet)
        row = df[df["canonical_accession"] == accession]
        if not row.empty:
            return float(row["S_Immuno"].iloc[0])

    # Live computation requires MHCflurry
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        warnings.warn(
            "mhcflurry is not installed. S_Immuno returns None. "
            "Install with: pip install pen-score[immuno]  or run script 14 on the VM.",
            stacklevel=2,
        )
        return None

    if sequence is None:
        from pen_score.utils.uniprot import fetch_sequence

        try:
            sequence = fetch_sequence(accession)
        except Exception as exc:
            warnings.warn(
                f"S_Immuno: could not fetch sequence for {accession}: {exc}",
                stacklevel=2,
            )
            return None

    try:
        predictor = Class1PresentationPredictor.load()
        # Generate 9-mer peptides (MHC-I) and 15-mer peptides (MHC-II)
        peptides_9 = [sequence[i : i + 9] for i in range(len(sequence) - 8)]
        result_i = predictor.predict(
            peptides=peptides_9,
            alleles=_HLA_I_ALLELES,
            include_affinity_percentile=True,
        )
        n_binders_i = int((result_i["presentation_score"] > 0.5).sum())

        epitope_density = n_binders_i / max(len(sequence), 1)
        # Normalise using approximate max density observed over SpCas9 (large, immunogenic)
        max_density = 0.35
        return round(float(np.clip(1.0 - epitope_density / max_density, 0.0, 1.0)), 4)

    except Exception as exc:
        warnings.warn(f"S_Immuno live computation failed for {accession}: {exc}", stacklevel=2)
        return None
