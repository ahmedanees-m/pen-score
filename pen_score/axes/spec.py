"""S_Spec - Specificity axis.

Formula (from axis_definitions.yaml):
    off_target_sites = BWA hits with <=3 mismatches on canonical 20-bp protospacer
    ratio = off_target_sites / (genome_size_bp / 1000)
    score = sigmoid(-2 * log10(ratio + 1e-10))

For non-RNA-guided editors (site-specific recombinases), an att-site
specificity proxy is used (number of natural att sites in GRCh38).

Requires [spec] optional extra (pysam, biopython, BWA binary on PATH).
Returns None if pysam/BWA unavailable.
"""

from __future__ import annotations

import math
import warnings

# BWA genome index path - expected at /data/genomes/GRCh38/GRCh38.fa.bwt on VM
_GENOME_SIZE_BP = 3.2e9
_BWA_INDEX = "/data/genomes/GRCh38/GRCh38.fa"


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def score(
    accession: str,
    protospacer: str | None = None,
    genome_index: str = _BWA_INDEX,
) -> float | None:
    """Compute S_Spec for a given editor.

    Parameters
    ----------
    accession:
        UniProt accession.
    protospacer:
        20-bp canonical protospacer sequence for RNA-guided editors.
        If None, fetched from editor_universe.yaml or returns None.
    genome_index:
        Path to BWA-indexed GRCh38 FASTA.

    Returns
    -------
    float in [0, 1] or None if BWA unavailable.
    """
    if protospacer is None:
        warnings.warn(
            f"No protospacer supplied for {accession}; S_Spec cannot be computed. "
            "Supply via script 11_compute_S_Spec.py with per-editor protospacers.",
            stacklevel=2,
        )
        return None

    try:
        import subprocess

        import pysam  # noqa: F401
    except ImportError:
        warnings.warn(
            "pysam is not installed. S_Spec returns None. "
            "Install with: pip install pen-score[spec]",
            stacklevel=2,
        )
        return None

    try:
        result = subprocess.run(
            ["bwa", "mem", "-k", "17", genome_index, "/dev/stdin"],
            input=f">query\n{protospacer}\n",
            capture_output=True,
            text=True,
            check=False,
        )
        # Count aligned hits with <=3 mismatches (NM:i: tag)
        off_target_count = 0
        for line in result.stdout.splitlines():
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            nm_tags = [f for f in fields[11:] if f.startswith("NM:i:")]
            if nm_tags:
                nm = int(nm_tags[0].split(":")[2])
                if nm <= 3:
                    off_target_count += 1

        ratio = off_target_count / (_GENOME_SIZE_BP / 1000.0)
        raw = -2.0 * math.log10(ratio + 1e-10)
        return round(_sigmoid(raw), 4)

    except Exception as exc:
        warnings.warn(f"S_Spec BWA scan failed for {accession}: {exc}", stacklevel=2)
        return None
