"""S_Energy - Energy Independence axis.

Scores whether an editor requires ATP hydrolysis for its core editing reaction.
Energy-independent editors (S_Energy = 1.0) act by site-specific chemistry
(transesterification, strand exchange) without a dedicated ATPase subunit.
Energy-dependent editors (S_Energy = 0.0) require a Walker A/B ATPase for
productive integration or translocation.

Method
------
1. Fetch primary protein sequence from UniProt REST API.
2. Scan for Walker A motif  : ``G[A-Z]{4}GK[ST]``   (P-loop consensus)
3. Scan for Walker B motif  : ``[LVIMF]{4}DE``        (β-strand; hydrophobic then DE)
4. If either motif is found -> 0.0 (ATP-dependent)
5. If neither is found      -> 1.0 (energy-independent)
6. Sequence unavailable     -> None

Override
--------
``walker_motif_override`` in ``editor_universe.yaml`` bypasses the sequence scan:

- ``true``  -> force S_Energy = 0.0  (ATP-dependent)
- ``false`` -> force S_Energy = 1.0  (energy-independent; use for NO_UNIPROT editors
              whose mechanism is confirmed energy-independent, e.g. SleepingBeauty)

Use the override for **multi-subunit systems** where the primary accession is
NOT the ATPase but the complex requires ATP to function.  The canonical examples
are Type V-K (ShCAST Cas12k) and Type I-F (VchCAST TnsB) CAST systems, where
TnsC (AAA+ ATPase, a *separate* subunit) drives target-site duplication; the
primary Cas12k or TnsB accession lacks Walker motifs on its own.

Limitations
-----------
- Walker A/B regex scan is a sequence-level proxy; it is not a structural
  pocket analysis (fpocket) or a biochemical assay.  Degenerate Walker motifs
  may occur by chance in non-ATPase proteins; true false-positive rate in the
  PEN-SCORE universe is expected to be low (1-2/28 editors) given the small,
  curated universe and the requirement that BOTH motifs co-occur for a Walker
  NTPase fold.
- The scan covers only the primary UniProt accession.  Multi-subunit systems
  where the ATPase is a separate polypeptide **must** use ``walker_motif_override``.
- REQUIRES_STEP7 sentinels return None; they are excluded from the composite.
"""

from __future__ import annotations

import re
import urllib.error
import urllib.request
import warnings

# Walker A (P-loop): G followed by any 4 residues, then GK, then S or T.
# Consensus: GxxxxGK[ST]
_WALKER_A: re.Pattern[str] = re.compile(r"G[A-Z]{4}GK[ST]")

# Walker B: four hydrophobic residues (LVIMF) followed by DE.
# Consensus: hhhhDE (h = hydrophobic)
_WALKER_B: re.Pattern[str] = re.compile(r"[LVIMF]{4}DE")


def _fetch_sequence(accession: str) -> str | None:
    """Return the plain amino-acid sequence for *accession* from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            fasta: str = resp.read().decode()
        lines = fasta.strip().split("\n")
        # First line is the FASTA header; join remaining lines as sequence.
        return "".join(lines[1:]).upper()
    except (urllib.error.URLError, OSError) as exc:
        warnings.warn(
            f"S_Energy: UniProt sequence fetch failed for {accession}: {exc}",
            stacklevel=3,
        )
        return None


def score(
    accession: str,
    walker_motif_override: bool | None = None,
) -> float | None:
    """Compute S_Energy for *accession*.

    Parameters
    ----------
    accession:
        UniProt accession of the primary editor protein.
        Special sentinels (``REQUIRES_STEP7``, ``NO_UNIPROT``) return ``None``
        unless *walker_motif_override* is provided.
    walker_motif_override:
        ``True``  -> return 0.0 immediately (ATP-dependent, no sequence scan).
        ``False`` -> return 1.0 immediately (energy-independent, no scan).
        ``None``  -> run the Walker A/B motif scan (default).

    Returns
    -------
    float in {0.0, 1.0}, or ``None`` if the sequence could not be fetched
    and no override was provided.
    """
    # Override takes priority - used for multi-subunit ATPase systems and
    # editors without a UniProt accession.
    if walker_motif_override is True:
        return 0.0
    if walker_motif_override is False:
        return 1.0

    # Sentinel accessions cannot be fetched from UniProt.
    if accession in ("REQUIRES_STEP7", "NO_UNIPROT", "", None):
        warnings.warn(
            f"S_Energy: sentinel accession '{accession}' - returning None. "
            "Set walker_motif_override in editor_universe.yaml if the mechanism "
            "class is known.",
            stacklevel=2,
        )
        return None

    seq = _fetch_sequence(accession)
    if seq is None:
        return None

    if _WALKER_A.search(seq) or _WALKER_B.search(seq):
        return 0.0
    return 1.0
