"""Compute S_Immuno (Immunogenicity axis).

Formula (axis_definitions.yaml v1.0.0):
    n_i  = unique 9-mer positions with netMHCpan-4.1 %Rank_EL < 0.5 (strong binders)
              over HLA-A*02:01, HLA-A*01:01, HLA-B*07:02, HLA-B*44:02
    n_ii = unique 15-mer positions with netMHCIIpan-4.0 %Rank_EL <= 10.0 (weak binders)
              over DRB1*01:01, DRB1*03:01, DRB1*04:01
    total = n_i + 0.5 * n_ii
    score = 1.0 - clip(total / max_total_over_universe, 0.0, 1.0)
        where max_total_over_universe = 95th-percentile of raw combined counts

NOTE: axis_definitions uses raw total (NOT divided by seq_length).
      Density normalisation was incorrect in v1/v2 runs - this is v3 (corrected).
      HLA-II alleles corrected: DRB1_0101/0301/0401 per axis_definitions (not 0301/0701/1501).

n_binders = unique positions where ANY allele binds (union per class).
max_density = 95th-percentile over all computed editors.
score = 1.0 - clip(epitope_load / max_density, 0, 1)   higher = less immunogenic.

Tools (DTU, academic license):
    netMHCpan-4.1:    ~/netmhc/netMHCpan-4.1/
    netMHCIIpan-4.0:  ~/netmhc/netMHCIIpan-4.0/

Run directly on the VM (NOT inside Docker):
    python3 ~/pen-stack/code/repos/pen-score/scripts/14_compute_S_Immuno.py \\
    2>&1 | tee ~/pen-stack/logs/pen-score/S_Immuno_$(date +%Y%m%d).log
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd
import requests
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"
OUT = Path.home() / "pen-stack" / "data" / "pen-score" / "axes" / "immuno"
OUT.mkdir(parents=True, exist_ok=True)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
REQ_DELAY_S = 0.3

# MHC-I (netMHCpan-4.1)
NMHOME_I = Path.home() / "netmhc" / "netMHCpan-4.1"
NMPLAT_I = NMHOME_I / "Linux_x86_64"
BIN_I = NMPLAT_I / "bin" / "netMHCpan"

# Alleles: remove '*' for the tool's format (HLA-A*02:01 -> HLA-A02:01)
HLA_I_ALLELES = ["HLA-A02:01", "HLA-A01:01", "HLA-B07:02", "HLA-B44:02"]
RANK_I = 0.5          # %Rank_EL threshold for Class I strong binder (axis_definitions v1.0.0)

# MHC-II (netMHCIIpan-4.0)
NMHOME_II = Path.home() / "netmhc" / "netMHCIIpan-4.0"
NMPLAT_II = NMHOME_II / "Linux_x86_64"
PERL_II = NMHOME_II / "NetMHCIIpan-4.0.pl"

HLA_II_ALLELES = ["DRB1_0101", "DRB1_0301", "DRB1_0401"]   # per axis_definitions v1.0.0
RANK_II = 10.0        # %Rank_EL threshold for Class II weak binder


def fetch_sequence(accession: str) -> str:
    base = (accession.split("_")[0]
            if "_" in accession and not accession.startswith("A0A")
            else accession)
    r = requests.get(f"{UNIPROT_BASE}/{base}.json", timeout=20)
    if r.status_code != 200:
        raise ValueError(f"UniProt HTTP {r.status_code} for {base}")
    return r.json()["sequence"]["value"]


def _write_fasta(sequence: str) -> str:
    """Write sequence to a temp FASTA file; return path (caller must delete)."""
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".fsa", delete=False)
    fh.write(f">seq\n{sequence}\n")
    fh.close()
    return fh.name


def run_netmhcpan_I(sequence: str) -> int:
    """Run netMHCpan-4.1 (Class I); return count of unique binding positions."""
    env = os.environ.copy()
    env["NETMHCpan"] = str(NMPLAT_I)
    env.setdefault("TMPDIR", "/tmp")

    fsa = _write_fasta(sequence)
    binding_positions: set[int] = set()
    try:
        for allele in HLA_I_ALLELES:
            cmd = [str(BIN_I), "-f", fsa, "-a", allele, "-l", "9"]
            result = subprocess.run(
                cmd, env=env, capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"netMHCpan failed for {allele}: {result.stderr[:200]}"
                )
            for line in result.stdout.splitlines():
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("-") \
                        or line.startswith("HLA"):
                    continue
                parts = line.split()
                # Columns: Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity
                #          Score_EL %Rank_EL [BindLevel]
                if len(parts) < 13:
                    continue
                try:
                    pos = int(parts[0])
                    rank = float(parts[12])
                except ValueError:
                    continue
                if rank <= RANK_I:
                    binding_positions.add(pos)
    finally:
        Path(fsa).unlink(missing_ok=True)

    return len(binding_positions)


def run_netmhcpan_II(sequence: str) -> int:
    """Run netMHCIIpan-4.0 (Class II); return count of unique binding positions."""
    env = os.environ.copy()
    env["NETMHCIIpan"] = str(NMHOME_II)
    env["NetMHCIIpanPLAT"] = str(NMPLAT_II)
    env.setdefault("TMPDIR", "/tmp")

    fsa = _write_fasta(sequence)
    binding_positions: set[int] = set()
    try:
        for allele in HLA_II_ALLELES:
            cmd = ["perl", str(PERL_II), "-a", allele, "-f", fsa]
            result = subprocess.run(
                cmd, env=env, capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"netMHCIIpan failed for {allele}: {result.stderr[:200]}"
                )
            for line in result.stdout.splitlines():
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("-"):
                    continue
                parts = line.split()
                # Columns: Pos MHC Peptide Of Core Core_Rel Identity
                #          Score_EL %Rank_EL Exp_Bind [BindLevel]
                if len(parts) < 9:
                    continue
                try:
                    pos = int(parts[0])
                    rank = float(parts[8])
                except ValueError:
                    continue
                if rank <= RANK_II:
                    binding_positions.add(pos)
    finally:
        Path(fsa).unlink(missing_ok=True)

    return len(binding_positions)


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Immuno computation v3 (MHC-I + MHC-II, raw total, corrected alleles)")
    print(f"  MHC-I  binary : {BIN_I}")
    print(f"  MHC-II Perl   : {PERL_II}")
    print(f"  HLA-I  alleles: {HLA_I_ALLELES}  (rank < {RANK_I}%  STRONG binders)")
    print(f"  HLA-II alleles: {HLA_II_ALLELES}  (rank <= {RANK_II}%  weak binders)")
    print(f"  Combined formula: n_I + 0.5*n_II  raw total (axis_definitions v1.0.0)")
    print("=" * 60)

    if not BIN_I.exists():
        sys.exit(f"ERROR: netMHCpan binary not found at {BIN_I}")
    if not PERL_II.exists():
        sys.exit(f"ERROR: netMHCIIpan Perl script not found at {PERL_II}")

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]

        if acc.startswith("REQUIRES") or acc == "NO_UNIPROT":
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "seq_length": None,
                "n_I_binders": None, "n_II_binders": None,
                "epitope_load": None, "S_Immuno": None, "note": "sentinel",
            })
            print(f"  {eid:<22} SKIP ({acc})")
            continue

        # Fetch sequence
        try:
            seq = fetch_sequence(acc)
            time.sleep(REQ_DELAY_S)
        except Exception as exc:
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "seq_length": None,
                "n_I_binders": None, "n_II_binders": None,
                "epitope_load": None, "S_Immuno": None,
                "note": f"uniprot_err: {str(exc)[:60]}",
            })
            print(f"  {eid:<22} UniProt ERROR: {exc}")
            continue

        # Run both tools
        try:
            t0 = time.time()
            n_I = run_netmhcpan_I(seq)
            n_II = run_netmhcpan_II(seq)
            elapsed = time.time() - t0
            # axis_definitions v1.0.0: total = n_i + 0.5 * n_ii (RAW, not per-residue density)
            combined = n_I + 0.5 * n_II
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "seq_length": len(seq),
                "n_I_binders": n_I, "n_II_binders": n_II,
                "epitope_load": round(combined, 4), "S_Immuno": None,
                "note": f"MHC-I+II; {elapsed:.0f}s",
            })
            print(f"  {eid:<22} len={len(seq):>5}  n_I={n_I:>4}  n_II={n_II:>4}"
                  f"  combined={combined:.1f}  [{elapsed:.0f}s]")
        except Exception as exc:
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "seq_length": len(seq),
                "n_I_binders": None, "n_II_binders": None,
                "epitope_load": None, "S_Immuno": None,
                "note": f"mhc_err: {str(exc)[:60]}",
            })
            print(f"  {eid:<22} MHC ERROR: {exc}")

        sys.stdout.flush()

    # Normalise: axis_definitions says max_total_over_universe; use 95th-pct for robustness
    computed = [r for r in rows if r["epitope_load"] is not None]
    if computed:
        totals = sorted(r["epitope_load"] for r in computed)
        p95_idx = max(0, int(math.ceil(0.95 * len(totals))) - 1)
        max_total = max(totals[p95_idx], 1.0)   # floor at 1 to avoid div/0
        print(f"\n95th-pct raw combined total = {max_total:.2f}  (normalisation denominator)")

        for row in rows:
            if row["epitope_load"] is not None:
                row["S_Immuno"] = round(
                    float(max(0.0, min(1.0, 1.0 - row["epitope_load"] / max_total))),
                    4,
                )
    else:
        max_total = None
        print("\nWARNING: no computed rows - cannot normalise")

    print("\nFinal scores:")
    for row in rows:
        sc = row["S_Immuno"]
        print(f"  {row['editor_id']:<22} S_Immuno={str(sc):<6}  [{row['note']}]")

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "immuno_scores.parquet", index=False)
    df.to_csv(OUT / "immuno_scores.csv", index=False)
    n_ok = df["S_Immuno"].notna().sum()
    n_skip = (df["note"] == "sentinel").sum()
    print(f"\nWritten -> {OUT}/immuno_scores.parquet (.csv)")
    print(f"  Computed: {n_ok}  |  Skipped (sentinel): {n_skip}  |  "
          f"Errors: {df['S_Immuno'].isna().sum() - n_skip}")
    print(f"  max_total (95th-pct) used: {max_total}")


if __name__ == "__main__":
    main()
