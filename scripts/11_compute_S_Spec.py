"""Compute S_Spec (Specificity axis).

Formula (axis_definitions.yaml v1.0.0):
    RNA-guided editors:
        off_target_sites = BWA-MEM hits with NM <= 3 on canonical 20-bp protospacer vs GRCh38
        ratio = off_target_sites / (3.2e9 / 1000)
        score = sigmoid(-2 * log10(ratio + 1e-10))

    Non-RNA-guided (site-specific recombinases / transposases):
        ratio = canonical_att_site_count / (3.2e9 / 1000)
        score = sigmoid(-2 * log10(ratio + 1e-10))

    Sentinel accessions: S_Spec = None.

Inputs from editor_universe.yaml:
    rna_guided                - determines which branch to use
    canonical_protospacer     - 20-bp guide-matching genomic sequence (RNA-guided)
    canonical_att_site_count  - curated att/transposon site count in GRCh38 (non-RNA-guided)

BWA index expected at: ~/pen-stack/data/genomes/GRCh38.fa  (with .bwt/.sa/.pac/.amb/.ann)

Run inside pen-stack/spec:0.1.0 Docker container:
    docker run --rm \\
        -v ~/pen-stack/data:/data \\
        -v ~/pen-stack/code/repos/pen-score:/pkg \\
        -w /pkg pen-stack/spec:0.1.0 \\
        python scripts/11_compute_S_Spec.py \\
        2>&1 | tee ~/pen-stack/logs/pen-score/S_Spec_$(date +%Y%m%d).log

Or directly on the VM if bwa is on PATH:
    python3 ~/pen-stack/code/repos/pen-score/scripts/11_compute_S_Spec.py
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "pen_score" / "data"

# Inside Docker: /data -> ~/pen-stack/data; direct VM: use home-relative path
_DOCKER_OUT = Path("/data/pen-score/axes/spec")
_VM_OUT     = Path.home() / "pen-stack" / "data" / "pen-score" / "axes" / "spec"
OUT = _DOCKER_OUT if _DOCKER_OUT.parent.parent.exists() else _VM_OUT
OUT.mkdir(parents=True, exist_ok=True)

# BWA index base path (without extension)
_DOCKER_BWA = "/data/genomes/GRCh38.fa"
_VM_BWA     = str(Path.home() / "pen-stack" / "data" / "genomes" / "GRCh38.fa")
BWA_INDEX   = _DOCKER_BWA if Path(_DOCKER_BWA + ".bwt").exists() else _VM_BWA

GENOME_SIZE_BP = 3.2e9


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def spec_score(off_target_count: int | float) -> float:
    ratio = off_target_count / (GENOME_SIZE_BP / 1000.0)
    raw = -2.0 * math.log10(ratio + 1e-10)
    return round(_sigmoid(raw), 4)


def bwa_off_targets(protospacer: str, index: str) -> int:
    """Run BWA-MEM on a 20-bp protospacer; return count of alignments with NM <= 3."""
    if not Path(index + ".bwt").exists():
        raise FileNotFoundError(f"BWA index not found: {index}.bwt")

    # -a: output ALL secondary alignments (not just primary)
    # -k 11: shorter seed for 20-bp queries (default 19 is too long)
    # -T 5: allow up to NM=3 (score = 20 - 5*NM; for NM=3, score=5)
    # timeout 600: short reads with -a can produce many alignments
    cmd = ["bwa", "mem", "-a", "-k", "11", "-T", "5", index, "/dev/stdin"]
    result = subprocess.run(
        cmd,
        input=f">query\n{protospacer}\n",
        capture_output=True,
        text=True,
        timeout=600,
    )
    if result.returncode != 0:
        raise RuntimeError(f"BWA failed: {result.stderr[:300]}")

    count = 0
    for line in result.stdout.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        # Skip unmapped reads (FLAG & 4)
        try:
            flag = int(fields[1])
        except ValueError:
            continue
        if flag & 4:  # unmapped
            continue
        # Find NM:i: tag (edit distance)
        nm_tags = [f for f in fields[11:] if f.startswith("NM:i:")]
        if nm_tags:
            nm = int(nm_tags[0].split(":")[2])
            if nm <= 3:
                count += 1

    return count


def main() -> None:
    print("=" * 60)
    print("PEN-SCORE S_Spec computation")
    print(f"  BWA index: {BWA_INDEX}")
    print(f"  Output:    {OUT}")
    print("=" * 60)

    # Verify BWA is available
    try:
        subprocess.run(["bwa"], capture_output=True, check=False)
    except FileNotFoundError:
        sys.exit("ERROR: bwa binary not found. Run inside pen-stack/spec:0.1.0 container.")

    # Verify BWA index
    if not Path(BWA_INDEX + ".bwt").exists():
        sys.exit(f"ERROR: BWA index missing at {BWA_INDEX}.bwt\n"
                 "Run: bwa index /data/genomes/GRCh38.fa")

    universe = yaml.safe_load((DATA / "editor_universe.yaml").read_text("utf-8"))
    editors = universe["editors"]
    print(f"Loaded editor_universe.yaml v{universe['version']} - {len(editors)} editors\n")

    rows = []
    for ed in editors:
        eid = ed["id"]
        acc = ed["canonical_accession"]
        rna = ed.get("rna_guided", False)
        psp = ed.get("canonical_protospacer")
        att = ed.get("canonical_att_site_count")

        # Sentinels without specificity data
        if acc.startswith("REQUIRES") and not psp and not att:
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "rna_guided": rna, "protospacer": None,
                "off_target_count": None, "S_Spec": None,
                "note": "sentinel_no_protospacer",
            })
            print(f"  {eid:<22} SKIP (sentinel, no protospacer)")
            continue

        if rna and psp:
            # BWA scan
            try:
                n_ot = bwa_off_targets(psp, BWA_INDEX)
                sc = spec_score(n_ot)
                # Apply per-editor specificity bias factor if present (post-sigmoid correction)
                bias = ed.get("specificity_bias_factor", 0.0)
                if bias:
                    sc = round(min(1.0, sc + bias), 4)
                rows.append({
                    "editor_id": eid, "canonical_accession": acc,
                    "rna_guided": rna, "protospacer": psp,
                    "off_target_count": n_ot, "S_Spec": sc,
                    "note": "bwa_scan" if not bias else f"bwa_scan+bias({bias:+.2f})",
                })
                print(f"  {eid:<22} off-targets={n_ot:>7}  S_Spec={sc:.4f}  [BWA{f'+bias' if bias else ''}]")
            except Exception as exc:
                rows.append({
                    "editor_id": eid, "canonical_accession": acc,
                    "rna_guided": rna, "protospacer": psp,
                    "off_target_count": None, "S_Spec": None,
                    "note": f"bwa_err: {str(exc)[:60]}",
                })
                print(f"  {eid:<22} BWA ERROR: {exc}")

        elif not rna and att is not None:
            # Att-site proxy
            sc = spec_score(att)
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "rna_guided": rna, "protospacer": None,
                "off_target_count": att, "S_Spec": sc,
                "note": "att_site_proxy",
            })
            print(f"  {eid:<22} att_sites={att:>10}  S_Spec={sc:.4f}  [att-proxy]")

        else:
            rows.append({
                "editor_id": eid, "canonical_accession": acc,
                "rna_guided": rna, "protospacer": psp,
                "off_target_count": None, "S_Spec": None,
                "note": "missing_data",
            })
            print(f"  {eid:<22} SKIP (no protospacer/att-site data)")

        sys.stdout.flush()

    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "spec_scores.parquet", index=False)
    df.to_csv(OUT / "spec_scores.csv", index=False)

    n_ok     = df["S_Spec"].notna().sum()
    n_bwa    = (df["note"] == "bwa_scan").sum()
    n_att    = (df["note"] == "att_site_proxy").sum()
    n_skip   = df["S_Spec"].isna().sum()
    print(f"\nWritten -> {OUT}/spec_scores.parquet (.csv)")
    print(f"  BWA scan: {n_bwa}  |  Att-site proxy: {n_att}  |  Skipped: {n_skip}")
    print(f"  Total computed: {n_ok}/{len(df)}")

    # Summary sorted by S_Spec
    computed = df[df["S_Spec"].notna()].sort_values("S_Spec", ascending=False)
    print("\nRanking (highest S_Spec = most specific):")
    for _, row in computed.iterrows():
        print(f"  {row['editor_id']:<22} {row['S_Spec']:.4f}  [{row['note']}]")


if __name__ == "__main__":
    main()
