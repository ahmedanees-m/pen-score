"""Pre-Registration Hash Lock.

Computes SHA-256 hashes of the four specification files that define
the PEN-SCORE scoring framework and embeds them in pre_registration.yaml.

This script must be run AFTER accession validation and BEFORE any
score-axis computation.  The output file is committed and tagged
`pre-registration-v1.0.0` - providing an immutable, timestamped record of
the exact specification against which all five retrospective predictions
will be evaluated.

Usage:
    python scripts/compute_prereg_hash.py [--verify]

    --verify   Re-compute hashes and check they match pre_registration.yaml.
               Exits 0 if match, 1 if mismatch.  Run before every axis-computation step.

Output:
    pen_score/data/pre_registration.yaml  (created or updated)
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR  = REPO_ROOT / "pen_score" / "data"
PREREG    = DATA_DIR / "pre_registration.yaml"

# Ordered list of files covered by the pre-registration hash.
# Do NOT reorder - ordering is part of the specification.
TARGETS: list[tuple[str, Path]] = [
    ("editor_universe.yaml",        DATA_DIR / "editor_universe.yaml"),
    ("axis_definitions.yaml",       DATA_DIR / "axis_definitions.yaml"),
    ("use_case_profiles.yaml",      DATA_DIR / "use_case_profiles.yaml"),
    ("scripts/02_validate_accessions.py",
                                    REPO_ROOT / "scripts" / "02_validate_accessions.py"),
]


def sha256_file(path: Path) -> str:
    """Return hex SHA-256 of file bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def combined_hash(file_hashes: dict[str, str]) -> str:
    """
    Stable combined hash: SHA-256 of the sorted (filename, hash) pairs
    joined as 'filename:hash\\n' lines.
    """
    h = hashlib.sha256()
    for fname in sorted(file_hashes):
        h.update(f"{fname}:{file_hashes[fname]}\n".encode())
    return h.hexdigest()


def compute() -> dict[str, str]:
    hashes: dict[str, str] = {}
    for label, path in TARGETS:
        if not path.exists():
            print(f"  ERROR: target file not found: {path}", file=sys.stderr)
            sys.exit(1)
        hashes[label] = sha256_file(path)
    return hashes


def verify_mode() -> None:
    """Re-compute and compare against committed pre_registration.yaml."""
    if not PREREG.exists():
        print("ERROR: pre_registration.yaml not found - run without --verify first.")
        sys.exit(1)

    committed = yaml.safe_load(PREREG.read_text(encoding="utf-8"))
    committed_hashes = committed.get("file_hashes", {})
    committed_combined = committed.get("combined_sha256", "")

    current_hashes  = compute()
    current_combined = combined_hash(current_hashes)

    ok = True
    print("Verifying pre-registration hashes ...")
    for label, _ in TARGETS:
        exp = committed_hashes.get(label, "<missing>")
        got = current_hashes.get(label, "<missing>")
        match = "OK" if exp == got else "MISMATCH"
        print(f"  {match}  {label}")
        if exp != got:
            print(f"         expected: {exp}")
            print(f"         current:  {got}")
            ok = False

    if committed_combined != current_combined:
        print(f"\n  COMBINED hash mismatch")
        print(f"    expected: {committed_combined}")
        print(f"    current:  {current_combined}")
        ok = False

    if ok:
        print("\nPRE-REGISTRATION VERIFIED: specification files unchanged.")
        sys.exit(0)
    else:
        print(
            "\nPRE-REGISTRATION MISMATCH: one or more specification files have been "
            "modified after pre-registration lock.\n"
            "If this change is intentional (e.g. correcting a sentinel accession),\n"
            "you MUST re-run without --verify to issue a new pre_registration.yaml\n"
            "and create a new tag pre-registration-v1.0.1."
        )
        sys.exit(1)


def write_mode() -> None:
    """Compute hashes and write pre_registration.yaml."""
    hashes   = compute()
    combined = combined_hash(hashes)
    ts       = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    doc = {
        "# Pre-registration lock for PEN-SCORE": None,
        "locked_at": ts,
        "pen_score_version": "0.0.1",
        "editor_universe_version": "1.0.2",
        "axis_definitions_version": "1.0.0",
        "n_editors_verified": 24,   # 28 total minus 4 REQUIRES_STEP7 sentinels
        "n_editors_sentinel": 4,    # evoCAST, NlovFz2, enNlovFz2, MmeFz2 - sequences in paper SIs
        "pre_registered_predictions": [
            {
                "id": 1,
                "editor": "evoCAST",
                "prediction": "Ranks in top-5 of AAV-deliverable DSB-free integrases",
                "threshold": "Top 5 of ~10-15 editors in that subset",
                "use_case": "human_therapeutic_aav_insertion",
            },
            {
                "id": 2,
                "editor": "IS621",
                "prediction": "Ranks in top-3 of programmable DSB-free systems",
                "threshold": "Top 3 of ~5-10 editors in that subset",
                "use_case": "human_therapeutic_aav_insertion",
            },
            {
                "id": 3,
                "editor": "SpCas9",
                "prediction": "Ranks in bottom 30% of PenScore (human therapeutic + AAV)",
                "threshold": "Bottom 9 of 30 editors",
                "use_case": "human_therapeutic_aav_insertion",
            },
            {
                "id": 4,
                "editor": "enNlovFz2",
                "prediction": "S_Deliv strictly greater than NlovFz2 WT",
                "threshold": "enNlovFz2.S_Deliv > NlovFz2.S_Deliv",
                "use_case": "any",
            },
            {
                "id": 5,
                "editor": "SpuFz1_V4",
                "prediction": "S_Spec strictly greater than SpuFz1 WT",
                "threshold": "SpuFz1_V4.S_Spec > SpuFz1.S_Spec",
                "use_case": "any",
            },
        ],
        "outcome_policy": {
            "5/5": "Strong claim supported",
            "4/5": "Report which prediction failed",
            "3/5": "reframe as scoring framework + lessons",
            "<=2/5": "Framework needs structural rework",
        },
        "sentinel_resolution_paths": {
            "evoCAST":   "Witte et al. 2025 Science SI Table S1 - PACE-evolved PseCAST TnsB",
            "NlovFz2":   "Wei et al. 2025 Nat Chem Biol SI - BLAST vs XP_04454xxxx (Naegleria lovaniensis ATCC 30569)",
            "enNlovFz2": "Same as NlovFz2",
            "MmeFz2":    "Saito et al. 2023 Nature SI Extended Data - MmeFanzor2 from Mercenaria mercenaria; A0A803GVX0 is 404 in UniProt; retrieve protein_id from paper SI Table",
        },
        "file_hashes": hashes,
        "combined_sha256": combined,
        "notes": (
            "combined_sha256 = SHA-256 of sorted 'filename:hash\\n' pairs. "
            "Run `python scripts/compute_prereg_hash.py --verify` before each axis-computation step "
            "to confirm specification files are unchanged."
        ),
    }

    # Write clean YAML (skip None-value comment keys)
    lines = ["# Pre-registration lock for PEN-SCORE",
             "# DO NOT EDIT - generated by scripts/compute_prereg_hash.py",
             "# Tag: pre-registration-v1.0.1", ""]
    real_doc = {k: v for k, v in doc.items() if v is not None and not k.startswith("#")}
    lines.append(yaml.dump(real_doc, default_flow_style=False, sort_keys=False, allow_unicode=True))
    PREREG.write_text("\n".join(lines), encoding="utf-8")

    print("Pre-registration hashes computed:")
    for label, h in hashes.items():
        print(f"  {label:<45} {h}")
    print(f"\n  combined_sha256: {combined}")
    print(f"\nWritten -> {PREREG}")
    print("\nNext steps:")
    print("  git add pen_score/data/pre_registration.yaml scripts/compute_prereg_hash.py")
    print("  git commit -m 'feat(prereg): re-lock pre-registration hashes v1.0.1 (IscB+MmeFz2 fixes)'")
    print("  git tag pre-registration-v1.0.1")
    print("  git push && git push --tags")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="Verify existing pre_registration.yaml; do not update.")
    args = parser.parse_args()
    if args.verify:
        verify_mode()
    else:
        write_mode()


if __name__ == "__main__":
    main()
