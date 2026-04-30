"""Accession Validation Gate.

Validates all entries in pen_score/data/editor_universe.yaml against live databases
(UniProt REST API and NCBI Entrez).  Must PASS before any score-axis script is run.

Three special sentinel values in canonical_accession are handled explicitly:

  NO_UNIPROT        Lab-engineered or reconstructed protein with no UniProt entry.
                    Validated via ncbi_protein_accession (if present) and
                    protein_length_aa (hardcoded).  Does NOT block the gate.

  REQUIRES_STEP7    Accession not yet resolved from primary paper SI.
                    BLOCKS the gate (exit code 2) until replaced with a real accession.

  REQUIRES_STEP7_EN Same semantics; used for engineered variants.

Exit codes:
  0  All resolvable accessions valid; any NO_UNIPROT or WARN entries flagged in report.
  1  One or more resolvable accessions fail UniProt/NCBI lookup (wrong protein, 404, etc.).
  2  One or more REQUIRES_STEP7 entries are still unresolved.

Run:
    docker run --rm \\
        -v ~/pen-stack/data:/data \\
        -v ~/pen-stack/code/repos/pen-score:/pkg \\
        -w /pkg pen-stack/biophysics:0.1.0 \\
        bash -c "pip install -e '.' --quiet && python scripts/02_validate_accessions.py" \\
        2>&1 | tee ~/pen-stack/logs/pen-score/step7_accession_validation_$(date +%Y%m%d).log
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import requests
import yaml

# Paths
REPO_ROOT  = Path(__file__).resolve().parent.parent
DATA_DIR   = REPO_ROOT / "pen_score" / "data"
REPORT_DIR = Path("/data/pen-score/validation")

# Constants
REQ_DELAY_S       = 0.25      # conservative: 4 req/s to stay within UniProt limits
UNIPROT_BASE      = "https://rest.uniprot.org/uniprotkb"
NCBI_ESUMMARY     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

SENTINEL_NO_UNIPROT = "NO_UNIPROT"
SENTINEL_REQUIRES   = "REQUIRES_STEP7"


# Result dataclass
@dataclass
class ValResult:
    editor_id:  str
    accession:  str
    status:     str           # PASS | WARN | FAIL | SKIP_NO_UNIPROT | UNRESOLVED
    length_aa:  Optional[int] = None
    organism:   Optional[str] = None
    taxid_got:  Optional[int] = None
    taxid_exp:  Optional[int] = None
    error:      Optional[str] = None
    notes:      str           = ""


# Validation helpers
def _get(url: str, **kw) -> requests.Response:
    """GET with timeout; raises on network error."""
    return requests.get(url, timeout=15, **kw)


def validate_uniprot(accession: str, expected_taxid: int) -> ValResult:
    """
    Validate a UniProt accession (TrEMBL or SwissProt).
    Strips engineered-variant suffixes (_PE2, _V4, _BE3 ...) before lookup.
    Returns PASS, WARN (taxid mismatch), or FAIL.
    """
    # Strip variant suffixes: Q99ZW2_PE2 -> Q99ZW2; A0A0L0H5U9_V4 -> A0A0L0H5U9
    # but NOT A0A2X3M8B0 (TrEMBL IDs contain digits after the second letter block)
    base = accession.split("_")[0] if "_" in accession and not accession.startswith("A0A") else accession
    if "_" in base:
        base = base.split("_")[0]

    url = f"{UNIPROT_BASE}/{base}.json"
    try:
        r = _get(url)
        if r.status_code == 404:
            return ValResult(
                editor_id="", accession=accession, status="FAIL",
                error=f"UniProt 404: accession {base} not found",
            )
        if r.status_code != 200:
            return ValResult(
                editor_id="", accession=accession, status="FAIL",
                error=f"UniProt HTTP {r.status_code}",
            )
        data      = r.json()
        length    = data["sequence"]["length"]
        sci_name  = data["organism"]["scientificName"]
        taxid_got = data["organism"]["taxonId"]

        if taxid_got != expected_taxid:
            return ValResult(
                editor_id="", accession=accession, status="WARN",
                length_aa=length, organism=sci_name,
                taxid_got=taxid_got, taxid_exp=expected_taxid,
                error=(
                    f"taxid mismatch: YAML={expected_taxid}, "
                    f"UniProt={taxid_got} ({sci_name})"
                ),
            )
        return ValResult(
            editor_id="", accession=accession, status="PASS",
            length_aa=length, organism=sci_name,
            taxid_got=taxid_got, taxid_exp=expected_taxid,
        )
    except Exception as exc:
        return ValResult(
            editor_id="", accession=accession, status="FAIL",
            error=f"Exception: {exc}",
        )


def validate_ncbi_protein(ncbi_acc: str, expected_length: int) -> ValResult:
    """
    Validate via NCBI protein accession (e.g. YCX28314.1).
    Returns PASS/WARN/FAIL with length.
    """
    try:
        r = _get(NCBI_ESUMMARY, params={"db": "protein", "id": ncbi_acc, "retmode": "json"})
        if r.status_code != 200:
            return ValResult(
                editor_id="", accession=ncbi_acc, status="FAIL",
                error=f"NCBI HTTP {r.status_code}",
            )
        data    = r.json()
        entries = [v for v in data["result"].values() if isinstance(v, dict)]
        if not entries:
            return ValResult(
                editor_id="", accession=ncbi_acc, status="FAIL",
                error="NCBI: no record returned",
            )
        slen    = int(entries[0].get("slen", 0))
        title   = entries[0].get("title", "")
        if expected_length and abs(slen - expected_length) > 5:
            return ValResult(
                editor_id="", accession=ncbi_acc, status="WARN",
                length_aa=slen,
                error=f"length mismatch: hardcoded={expected_length}, NCBI={slen} ({title})",
            )
        return ValResult(
            editor_id="", accession=ncbi_acc, status="PASS",
            length_aa=slen,
            notes=f"NCBI protein: {title[:80]}",
        )
    except Exception as exc:
        return ValResult(
            editor_id="", accession=ncbi_acc, status="FAIL",
            error=f"Exception: {exc}",
        )


# Main
def main() -> None:
    print("=" * 70)
    print("PEN-SCORE: Mandatory Accession Validation Gate")
    print("=" * 70)

    yaml_path = DATA_DIR / "editor_universe.yaml"
    data      = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    editors   = data["editors"]
    version   = data.get("version", "unknown")
    print(f"\nLoaded {yaml_path.name} v{version} - {len(editors)} entries\n")

    results:    list[ValResult] = []
    unresolved: list[str]       = []

    for entry in editors:
        eid    = entry["id"]
        acc    = entry.get("canonical_accession", "")
        taxid  = entry.get("organism_taxid", 0)

        # 1. REQUIRES_STEP7 - gate-blocking unresolved accession
        if acc.startswith("REQUIRES_STEP7") or acc.startswith("REQUIRES_"):
            res = ValResult(
                editor_id=eid, accession=acc, status="UNRESOLVED",
                notes=(
                    "Must retrieve from primary paper SI and replace sentinel "
                    "before accession-validation gate can pass."
                ),
            )
            results.append(res)
            unresolved.append(eid)
            doi = entry.get("primary_doi", "no DOI")
            print(f"  [UNRESOLVED] {eid:<22} {acc}  ->  DOI {doi}")
            continue

        # 2. NO_UNIPROT - engineered/reconstructed; use NCBI if available
        if acc == SENTINEL_NO_UNIPROT:
            ncbi_acc    = entry.get("ncbi_protein_accession")
            hard_len    = entry.get("protein_length_aa")
            if ncbi_acc:
                res = validate_ncbi_protein(ncbi_acc, hard_len or 0)
                res.notes = f"NO_UNIPROT; validated via NCBI {ncbi_acc}"
                time.sleep(REQ_DELAY_S)
            else:
                res = ValResult(
                    editor_id="", accession=acc, status="SKIP_NO_UNIPROT",
                    length_aa=hard_len,
                    notes=f"NO_UNIPROT; using hardcoded protein_length_aa={hard_len}",
                )
            res.editor_id = eid
            results.append(res)
            sym = "ok" if res.status in ("PASS", "WARN", "SKIP_NO_UNIPROT") else "xx"
            extra = f"{res.length_aa} aa" if res.length_aa else ""
            errmsg = f"  [{sym} {res.status:<16}] {eid:<22} {extra}  {res.notes}"
            print(errmsg)
            continue

        # 3. Standard UniProt accession (with optional variant suffix)
        res           = validate_uniprot(acc, taxid)
        res.editor_id = eid
        results.append(res)
        time.sleep(REQ_DELAY_S)

        sym = {"PASS": "ok", "WARN": "!!", "FAIL": "xx"}.get(res.status, "?")
        detail = (
            f"{res.length_aa} aa, {res.organism}"
            if res.status in ("PASS", "WARN") and res.length_aa
            else (res.error or "")
        )
        print(f"  [{sym} {res.status:<16}] {eid:<22} {acc:<24} {detail}")

    # Summary
    n_pass       = sum(1 for r in results if r.status in ("PASS", "SKIP_NO_UNIPROT"))
    n_warn       = sum(1 for r in results if r.status == "WARN")
    n_fail       = sum(1 for r in results if r.status == "FAIL")
    n_unresolved = len(unresolved)

    print()
    print("-" * 70)
    print(f"  Total editors  : {len(results)}")
    print(f"  PASS           : {n_pass}")
    print(f"  WARN (taxid)   : {n_warn}")
    print(f"  FAIL           : {n_fail}")
    print(f"  UNRESOLVED: {n_unresolved}  <- must be resolved before axis computation")

    if unresolved:
        print()
        print("  Editors with unresolved accessions:")
        for eid in unresolved:
            entry = next(e for e in editors if e["id"] == eid)
            doi   = entry.get("primary_doi", "no DOI")
            ref   = entry.get("primary_reference", "")
            print(f"    {eid:<22} {ref}  (DOI: {doi})")

    if n_warn:
        print()
        print("  WARN - taxid mismatches (review but do not block gate):")
        for r in results:
            if r.status == "WARN":
                print(f"    {r.editor_id:<22} {r.error}")

    # Write JSON report
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report = {
        "editor_universe_version": version,
        "n_editors":    len(results),
        "n_pass":       n_pass,
        "n_warn":       n_warn,
        "n_fail":       n_fail,
        "n_unresolved": n_unresolved,
        "unresolved_editors": unresolved,
        "results": [asdict(r) for r in results],
    }
    out = REPORT_DIR / "accession_validation_report.json"
    out.write_text(json.dumps(report, indent=2))
    print()
    print(f"  Report written -> {out}")
    print("=" * 70)

    # Exit codes
    if n_fail > 0:
        print("\nGATE: FAIL  - accession errors must be corrected before axis computation.")
        sys.exit(1)
    if n_unresolved > 0:
        print(
            "\nGATE: BLOCKED  - unresolved REQUIRES_STEP7 entries must be replaced "
            "with verified accessions.\nRe-run this script after each fix."
        )
        sys.exit(2)
    print("\nGATE: PASS  - all resolvable accessions verified. Safe to proceed to axis computation.")
    sys.exit(0)


if __name__ == "__main__":
    main()
