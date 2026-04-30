"""Smoke test: verify GENOME-ATLAS and MECH-CLASS are importable and artifacts load.

Run via:
    docker run --rm \\
        -v ~/pen-stack/data:/data \\
        -v ~/pen-stack/code/repos/pen-score:/pkg \\
        -w /pkg pen-stack/biophysics:0.1.0 \\
        bash -c "pip install -e '.[mech-class,atlas]' --quiet && python scripts/00_smoke_import.py"
"""

from __future__ import annotations

import sys
from pathlib import Path

DATA_DIR = Path("/data")

print("=" * 60)
print("PEN-SCORE Smoke Import Test")
print("=" * 60)

# -----------------------------------------------------------------------
# 1. pen-score itself
# -----------------------------------------------------------------------
print("\n[1] Importing pen_score ...")
import pen_score

print(f"    pen_score {pen_score.__version__}")

from pen_score.data.loader import load_editor_universe, load_use_case_profiles

editors = load_editor_universe()
profiles = load_use_case_profiles()
print(f"    Editor universe: {len(editors)} editors loaded")
print(f"    Use-case profiles: {list(profiles.keys())}")

# Quick axis smoke tests (no external deps)
from pen_score.axes import cargo, deliv, prog

is621_cargo = cargo.score("IS621")
is621_deliv = deliv.score("IS621", total_aa=300)
is621_prog = prog.score("IS621")
spcas9_prog = prog.score("SpCas9")
cre_prog = prog.score("Cre")
assert is621_cargo is not None and is621_cargo > 0.9, f"IS621 S_Cargo unexpected: {is621_cargo}"
assert is621_deliv is not None and is621_deliv > 0.9, f"IS621 S_Deliv unexpected: {is621_deliv}"
assert is621_prog == 1.0, f"IS621 S_Prog unexpected: {is621_prog}"
assert spcas9_prog == 1.0, f"SpCas9 S_Prog unexpected: {spcas9_prog}"
assert cre_prog == 0.0, f"Cre S_Prog unexpected: {cre_prog}"
print("    Axis smoke checks PASS (S_Cargo, S_Deliv, S_Prog)")

# -----------------------------------------------------------------------
# 2. GENOME-ATLAS - genome_atlas
# -----------------------------------------------------------------------
print("\n[2] Importing genome_atlas ...")
try:
    import genome_atlas

    print(f"    genome_atlas {genome_atlas.__version__}")

    atlas_db = DATA_DIR / "graphs" / "atlas.duckdb"
    esm2_parquet = DATA_DIR / "embeddings" / "esm2_150M_v6.parquet"
    graphsage_parquet = DATA_DIR / "embeddings" / "graphsage_v6.parquet"

    assert atlas_db.exists(), f"MISSING: {atlas_db}"
    assert esm2_parquet.exists(), f"MISSING: {esm2_parquet}"
    assert graphsage_parquet.exists(), f"MISSING: {graphsage_parquet}"
    print(f"    atlas.duckdb:          {atlas_db.stat().st_size / 1e6:.1f} MB")
    print(f"    esm2_150M_v6.parquet:  {esm2_parquet.stat().st_size / 1e6:.1f} MB")
    print(f"    graphsage_v6.parquet:  {graphsage_parquet.stat().st_size / 1e6:.1f} MB")

except ImportError:
    print("    genome_atlas not installed (install with: pip install pen-score[atlas])")
    print("    Skipping GENOME-ATLAS artifact check.")

# -----------------------------------------------------------------------
# 3. MECH-CLASS - mech_class
# -----------------------------------------------------------------------
print("\n[3] Importing mech_class ...")
try:
    import mech_class

    print(f"    mech_class {mech_class.__version__}")

    tier_a_pkl = DATA_DIR / "models" / "tier_a" / "model.pkl"
    composite_pkl = DATA_DIR / "models" / "composite_head" / "model.pkl"

    assert tier_a_pkl.exists(), f"MISSING: {tier_a_pkl}"
    assert composite_pkl.exists(), f"MISSING: {composite_pkl}"
    print(f"    tier_a/model.pkl:       {tier_a_pkl.stat().st_size / 1e3:.0f} KB")
    print(f"    composite_head/model.pkl: {composite_pkl.stat().st_size / 1e3:.0f} KB")

    from mech_class.api import Predictor

    predictor = Predictor.load(model_dir=DATA_DIR / "models")
    # Smoke prediction: verify API runs; full-sequence accuracy checked at full validation.
    # Stub sequence - too short for correct ESM2 embedding; class not asserted.
    _SPCAS9_SEQ_STUB = "MDKKYSIGLDIGTNSVGWAV"
    pred = predictor.predict_from_sequence(
        "Q99ZW2",
        _SPCAS9_SEQ_STUB,
        pfam_hits=["PF18541", "PF16595", "PF18516"],
    )
    assert hasattr(pred, "tier_a") and pred.tier_a is not None, f"predict_from_sequence returned invalid result: {pred}"
    print(f"    SpCas9 smoke prediction: {pred.tier_a} (conf {pred.tier_a_confidence:.3f}) API OK (full-seq accuracy at full validation)")

except ImportError:
    print("    mech_class not installed (install with: pip install pen-score[mech-class])")
    print("    Skipping MECH-CLASS artifact check.")
    print("    NOTE: S_DSB and S_Prog axes will return None without mech-class.")

# -----------------------------------------------------------------------
# 4. biophysics image check
# -----------------------------------------------------------------------
print("\n[4] Checking biophysics dependencies (xtb, APBS, fpocket) ...")
try:
    import xtb

    print(f"    xtb-python: {xtb.__version__}")
except ImportError:
    print("    xtb-python not available (expected outside biophysics container)")

import subprocess

for tool in ["apbs", "fpocket"]:
    try:
        out = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=5)
        v = (out.stdout + out.stderr).strip().split("\n")[0]
        print(f"    {tool}: {v}")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print(f"    {tool}: not found (expected outside biophysics container)")

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
print("\n" + "=" * 60)
print("Smoke import test complete.")
print("Prerequisite check: PASS")
print("=" * 60)
