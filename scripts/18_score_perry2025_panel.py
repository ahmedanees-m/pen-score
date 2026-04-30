"""Score Perry-2025 IS110 ortholog panel.

Perry et al. 2026 Science (doi: 10.1126/science.adz0276) screened ~72 IS622
orthologs in human HEK293 cells and reported insertion efficiencies.  This
script scores the top-performing orthologs (>=5% insertion efficiency) through
the 8-axis pen-score pipeline and appends them to editor_universe.yaml v1.0.6.

Data requirement
----------------
Place the Perry 2025 ortholog CSV at::

    data/perry2025_supplementary_orthologs.csv

Expected columns:
    name                      : short editor identifier (e.g. "IS622_Ecl")
    accession                 : UniProt accession (or REQUIRES_STEP7 if pending)
    organism                  : source organism name
    organism_taxid            : NCBI taxid (int)
    insertion_efficiency_pct  : % insertion in HEK293 cells (primary screen)
    cargo_capacity_bp         : max cargo in bp (use 930000 as IS622-class default)
    primary_doi               : primary reference DOI
    notes                     : free text

Output
------
- ``results/perry2025_penscore_panel.csv``  - axis scores + PenScore for each
  ortholog that passes the >=5% threshold.
- Appended YAML entries printed to stdout for review before adding to
  pen_score/data/editor_universe.yaml.

Usage
-----
::

    python scripts/18_score_perry2025_panel.py

"""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).parent.parent
DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
UNIVERSE_PATH = REPO_ROOT / "pen_score" / "data" / "editor_universe.yaml"
PERRY_CSV = DATA_DIR / "perry2025_supplementary_orthologs.csv"

THRESHOLD_PCT = 5.0  # minimum insertion efficiency to include

# Default values for IS110-class editors (all are bridge recombinases)
IS110_DEFAULTS: dict[str, object] = {
    "mechanism_bucket": "DSB_FREE_TRANSEST_RECOMBINASE",
    "rna_guided": True,
    "composite_architecture": True,
    "cargo_capacity_bp": 930000,  # IS622 megabase class default
    "primary_reference": "Perry et al. 2026 Science",
    "primary_doi": "10.1126/science.adz0276",
    "references_used_for_pubmed": [
        "IS622 ortholog",
        "bridge recombinase",
        "IS110 recombinase human cells",
    ],
}

# IS621 proxy values for axes that require live computation
IS621_PROXY = {
    "S_DSB": 1.0,       # IS110 Tier-A gate (PF01548+PF02371)
    "S_Energy": 1.0,    # No Walker motifs; IS110 is energy-independent
    "S_Mature": 0.0,    # New proteins; 0 clinical PubMed hits at scoring time
    # S_Spec / S_Immuno use IS621 proxy values when proxy_source="IS621"
    "S_Spec_proxy": 0.9891,    # IS621 BWA-MEM S_Spec (from v0.1.1 scorecard)
    "S_Immuno_proxy": 0.7594,  # IS621 MHCflurry S_Immuno (from v0.1.1 scorecard)
    "S_Prog": 1.0,      # RNA-guided -> programmable
}


def _s_cargo(cargo_bp: float) -> float:
    """S_Cargo = log10(cargo_bp) / log10(1e6), clipped [0, 1]."""
    import math

    score = math.log10(max(cargo_bp, 1)) / math.log10(1e6)
    return float(min(1.0, max(0.0, score)))


def _s_deliv_from_uniprot(accession: str) -> float | None:
    """S_Deliv = sigmoid on protein length.  Returns None if unavailable."""
    import math

    from pen_score.utils.uniprot import fetch_sequence_length

    try:
        length = fetch_sequence_length(accession)
        if length is None:
            return None
        return 1.0 / (1.0 + math.exp(0.005 * (length - 900)))
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"S_Deliv unavailable for {accession}: {exc}", stacklevel=2)
        return None


def _compute_pen_score(axes: dict[str, float | None], weights: dict[str, float]) -> float | None:
    from pen_score.scorer.composite import compute_pen_score

    score, _ = compute_pen_score(axes, weights)
    return score


def score_perry_panel(
    csv_path: Path = PERRY_CSV,
    threshold_pct: float = THRESHOLD_PCT,
) -> pd.DataFrame:
    """Score all orthologs in csv_path with insertion_efficiency_pct >= threshold_pct.

    Parameters
    ----------
    csv_path:
        Path to Perry supplementary CSV (see module docstring for columns).
    threshold_pct:
        Minimum insertion efficiency (%) to include.

    Returns
    -------
    DataFrame with columns: name, accession, insertion_efficiency_pct, PenScore,
    S_DSB, S_Spec, S_Cargo, S_Deliv, S_Immuno, S_Prog, S_Mature, S_Energy.
    """
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Perry supplementary CSV not found at {csv_path}.\n"
            "Download the supplementary table from:\n"
            "  doi: 10.1126/science.adz0276 (Perry et al. 2026 Science)\n"
            "and place it at data/perry2025_supplementary_orthologs.csv"
        )

    perry = pd.read_csv(csv_path)
    candidates = perry[perry["insertion_efficiency_pct"] >= threshold_pct].copy()
    print(f"Found {len(candidates)} orthologs with insertion_efficiency >= {threshold_pct}%")

    from pen_score.data.loader import load_use_case_profiles

    weights = load_use_case_profiles()["human_therapeutic_aav_insertion"]
    rows: list[dict] = []

    for _, row in candidates.iterrows():
        acc = str(row["accession"])
        name = str(row["name"])
        cargo_bp = float(row.get("cargo_capacity_bp", 930000))

        s_deliv = _s_deliv_from_uniprot(acc) if not acc.startswith(("REQUIRES_", "NO_")) else None

        axes: dict[str, float | None] = {
            "S_DSB": IS621_PROXY["S_DSB"],
            "S_Spec": IS621_PROXY["S_Spec_proxy"],
            "S_Cargo": _s_cargo(cargo_bp),
            "S_Deliv": s_deliv,
            "S_Immuno": IS621_PROXY["S_Immuno_proxy"],
            "S_Prog": IS621_PROXY["S_Prog"],
            "S_Mature": IS621_PROXY["S_Mature"],
            "S_Energy": IS621_PROXY["S_Energy"],
        }
        pen_score = _compute_pen_score(axes, weights)

        rows.append(
            {
                "name": name,
                "accession": acc,
                "organism": row.get("organism", ""),
                "insertion_efficiency_pct": row["insertion_efficiency_pct"],
                "PenScore": pen_score,
                **axes,
            }
        )

    df = pd.DataFrame(rows).sort_values("PenScore", ascending=False).reset_index(drop=True)
    print(
        df[["name", "PenScore", "S_DSB", "S_Cargo", "S_Deliv", "S_Energy", "S_Mature"]].to_string(
            index=False
        )
    )
    return df


def print_yaml_entries(df: pd.DataFrame) -> None:
    """Print editor_universe.yaml entries for review."""
    print("\n# --- Perry-2025 ortholog entries for editor_universe.yaml ---\n")
    for _, row in df.iterrows():
        entry = {
            "id": str(row["name"]),
            "canonical_accession": str(row["accession"]),
            "organism": str(row.get("organism", "")),
            "mechanism_bucket": IS110_DEFAULTS["mechanism_bucket"],
            "rna_guided": IS110_DEFAULTS["rna_guided"],
            "composite_architecture": IS110_DEFAULTS["composite_architecture"],
            "cargo_capacity_bp": int(row.get("cargo_capacity_bp", 930000)),
            "year_discovered": 2025,
            "primary_reference": IS110_DEFAULTS["primary_reference"],
            "primary_doi": IS110_DEFAULTS["primary_doi"],
            "notes": (
                f"IS622 ortholog; insertion_efficiency={row['insertion_efficiency_pct']:.1f}% "
                f"in HEK293 (Perry 2026); S_Spec/S_Immuno use IS621 proxy"
            ),
            "references_used_for_pubmed": IS110_DEFAULTS["references_used_for_pubmed"],
        }
        print(yaml.dump([entry], default_flow_style=False, sort_keys=False))


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    df = score_perry_panel()
    output_path = RESULTS_DIR / "perry2025_penscore_panel.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved to {output_path}")

    print_yaml_entries(df)
    print(
        "\nReview the YAML entries above, then append to "
        "pen_score/data/editor_universe.yaml manually after verification."
    )


if __name__ == "__main__":
    main()
