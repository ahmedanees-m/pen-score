# Changelog

All notable changes to pen-score are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.3] - 2026-05-25

### Added

**`get_editor_metadata()` API**
- `pen_score/api.py`: new `get_editor_metadata(editor_id) -> EditorMetadata` function.
  Returns frozen `EditorMetadata` dataclass exposing the two PEN-COMPARE v3.2 boolean
  fields plus alias, UniProt, and organism metadata.  Resolves deprecated aliases
  (e.g. `"IS622"` -> `"ISCro4"`) with a `DeprecationWarning`.
- `pen_score/__init__.py`: `get_editor_metadata` and `EditorMetadata` exported at
  package level.
- `tests/unit/test_get_editor_metadata.py` (NEW - 18 tests):
  - ISCro4 canonical metadata (canonical_name, uniprot, intrinsic=True, cell_based=True,
    >=2 sources, IS622 in aliases)
  - IS621 calibration keystone (intrinsic=True, cell_based=False - distinguishes
    PROBABLE_WRITER from TRUE_WRITER in PEN-COMPARE v3.2)
  - SpCas9 Gate 3 exclusion (intrinsic=False - HDR-template-based)
  - IS622 alias resolution + DeprecationWarning emission
  - Return type, frozen immutability, KeyError on unknown ID
  - All 29 editors have bool fields

**editor_universe.yaml schema v1.0.7 - two new fields per editor**
- `pen_score/data/editor_universe.yaml` -> v1.0.7 (bumped from v1.0.6):
  - `intrinsic_cargo_mechanism: bool` - **True** if the editor carries/inserts cargo via
    its catalytic mechanism (IS110 bridge recombinases, CAST transposases, site-specific
    recombinases). **False** if cargo requires an external HDR donor template (SpCas9,
    prime editors).  Required by PEN-COMPARE v3.2 **Gate 3**.
  - `cell_based_evidence: bool` - **True** if peer-reviewed mammalian cell activity exists
    at >1% editing efficiency. **False** for in vitro / E. coli only.  Required by
    PEN-COMPARE v3.2 **TRUE_WRITER tier**.
  - `cell_based_sources: list[str]` - supporting citations (DOI strings or author-year)
    for editors with cell_based_evidence=True.
  - Applied to all 29 editors (15 intrinsic=True, 14 intrinsic=False; 22 cell_based=True,
    7 cell_based=False).
- `pen_score/data/loader.py`: `EditorEntry` model updated with new optional fields
  (`aliases: list[str]`, `intrinsic_cargo_mechanism: bool`,
  `cell_based_evidence: bool`, `cell_based_sources: list[str]`).

### Changed

**ISCro4 canonical naming - IS622 deprecated**
- `pen_score/data/editor_universe.yaml` v1.0.7: editor `IS622` renamed to `ISCro4`
  (canonical per UniProt D2TGM5 gene name + Pelea 2026 *Science*
  doi:10.1126/science.adz1884). `aliases: ["IS622"]` preserved for backward
  compatibility. IS622 was the Perry 2025 *bioRxiv* preprint label, retired upon
  publication.
- `pen_score/data/use_case_profiles.yaml`: `megabase_rearrangement` profile comment
  updated to reference ISCro4 (was IS622).
- `pyproject.toml`: dependency pins bumped to reflect PEN-STACK v3.2-compat versions:
  `mech-class>=0.5.4,<0.6.0` (ISCro4 holdout probe rename); `genome-atlas>=0.7.2,<0.8.0`
  (ISCro4 canonical in foundational_systems).

**Test suite**
- Total: **196 pass, 3 skipped (parquet-only), 0 failures** (was 178/3/0 in v0.1.2).
  18 new tests from `test_get_editor_metadata.py`.
- `ruff check` passes on all changed files.

### Compatibility

- Fully backward compatible. All existing `Scorer` / `score_editor` / `exclude_axes`
  API unchanged.
- `IS622` continues to resolve via the alias mechanism (with `DeprecationWarning`).
- Required by PEN-COMPARE v3.2 Steps 1+.

---

## [0.1.2] - 2026-05-24

### Added

**Dependency pin bumps**
- `pyproject.toml`: `mech-class>=0.5.2,<0.7.0` -> `>=0.5.3,<0.6.0` (picks up ISCro4/D2TGM5 OOD
  probe, atlas pin >=0.7.1). `genome-atlas>=0.6.0,<0.7.0` -> `>=0.7.1,<0.8.0` (picks up IS622
  addition, SIMILAR_TO/HAS_RNA/PART_OF secondary edges via `graph_view='full'`).

**Perry-2025 IS110 ortholog panel infrastructure**
- `scripts/18_score_perry2025_panel.py`: scoring script for IS110 orthologs from Perry et al.
  2026 Science (doi:10.1126/science.adz0276). Filters to >=5% insertion efficiency, applies
  IS621 proxy for S_Spec/S_Immuno, computes S_DSB=1.0 / S_Energy=1.0 / S_Mature=0.0 for all
  IS110-class orthologs. Full panel pending supplementary CSV acquisition.
- `pen_score/data/editor_universe.yaml` -> v1.0.6 (2026-05-24): Perry-2025 infrastructure note
  added; `notes` block updated.

**New use-case profiles**
- `pen_score/data/use_case_profiles.yaml` v0.1.2 (2026-05-24): two new profiles added.
  - `megabase_rearrangement`: weights S_Prog=0.30 (boundary selection) > S_Cargo=0.25 (flanking
    capacity) > S_DSB=0.20 > S_Energy=0.10; S_Mature=0.00 (all megabase editors are immature).
    Reference editor: IS622 (0.93 Mb inversion, HEK293, Perry 2026).
  - `therapeutic_excision_bcl11a`: weights S_DSB=S_Spec=0.25 (tied; DSB-free + exact targeting)
    > S_Immuno=0.20 (HSC therapy); targets BCL11A erythroid enhancer +58 element excision for
    SCD / β-thalassemia. Reference editor: IS621.
- `tests/unit/test_profiles.py` (NEW): 12 tests covering weight sums, profile count,
  axis completeness, semantic constraints (S_Mature=0 in megabase profile, S_DSB+S_Spec highest
  in BCL11A profile), and composite score floor tests.

**`exclude_axes` API - biophysical-only comparisons**
- `pen_score/api.py`: `_VALID_AXES` frozenset; `score_editor(exclude_axes=...)` parameter with
  ValueError on unrecognised axis names; docstring.
- `pen_score/api.py` `_composite()`: `exclude_axes` parameter - excluded axes skip silently
  (not recorded as missing); weight redistribution to remaining available axes via renormalisation.
- `pen_score/scorer/composite.py` `compute_pen_score()`: same `exclude_axes` parameter for
  callers using the standalone composite function.
- `pen_score/cli.py`: `--exclude-axes axis1,axis2` option on `score-editor` subcommand.
- `pen_score/data/axis_definitions.yaml`: `exclude_axes_api` documentation block.
- `docs/quickstart.rst`: IS622 biophysical-only example; updated use-case profile table (7 rows).
- `tests/unit/test_api.py`: 4 new `TestExcludeAxes` tests (biophysical lift, weight renorm,
  invalid axis ValueError, exclude-all -> None).

**Test suite**
- Total: **178 pass, 18 skipped (VM-only regression), 0 failures**.
  17 tests added (13 profile tests + 4 exclude_axes tests).
- `ruff check` + `ruff format --check` both pass (38 files clean).

---

## [0.1.1] - 2026-05-22

### Added

**S_Energy - Energy Independence axis (8th scoring axis, v0.1.1)**
- `pen_score/axes/energy.py`: Walker A (`G[A-Z]{4}GK[ST]`) + Walker B (`[LVIMF]{4}DE`)
  motif regex scan on primary UniProt FASTA sequence. Returns `1.0` (energy-independent)
  when neither motif is detected; `0.0` (ATP-dependent) when either motif is found or
  `walker_motif_override: true` is set in YAML.
- `pen_score/data/editor_universe.yaml` v1.0.5 (2026-05-22): `walker_motif_override` field
  added to multi-subunit CAST systems (evoCAST, CAST_VK, CAST_IF `-> true`; TnsC ATPase is
  a separate subunit not captured by the Cas12k/TnsB primary accession) and SleepingBeauty
  (`-> false`; DDE-family Tc1/mariner, no ATPase; NO_UNIPROT prevents sequence fetch).
- `pen_score/data/axis_definitions.yaml`: S_Energy definition added with Walker motif
  formula, inputs, dependencies, override rationale, and unit test oracle.
- `pen_score/data/use_case_profiles.yaml`: S_Energy weight = 0.05 across all 5 profiles;
  7-axis v0.1.0 weights proportionally rescaled x 0.95 so each profile still sums to 1.00.
- `pen_score/scorer/composite.py`: S_Energy added to `_AXES` (8 axes total); docstring updated.
- `pen_score/api.py`: `AxisScores.S_Energy` field; `_compute_axes_live()` calls
  `energy.score(walker_motif_override=...)`; `_generate_reasoning()` adds "Energy independence"
  label; `_default_weights()` updated to 8-axis v0.1.1 defaults.
- `pen_score/cli.py`: "Energy" column added to `select` Rich table; `_AXES` list updated.
- `scripts/17_compute_S_Energy.py`: pipeline script (Walker A/B scan over all 31 editors;
  `--dry-run` flag for preview without writing).
- `tests/unit/test_energy.py`: 28 new unit tests covering override logic, sentinel handling,
  Walker A/B motif detection (both variants), energy-independent editor set, CAST overrides,
  network failure fallback, AxisScores model validation, and composite integration.
- All pre-existing tests updated for 8-axis model: conftest.py `mock_axis_scores` +
  `default_weights`; test_api.py (3 tests); test_axes.py (2 tests); test_loader.py (1 test).
  Full suite: **158 pass + 3 env-dependent skips = 161 total, 0 failures**.

**Post-pipeline comparator**
- `pen_score/data/editor_universe.yaml` v1.0.4 (2026-05-21): IS622 (ISFinder: ISCro4,
  UniProt D2TGM5, 326 aa) added as post-pipeline comparator. IS622 is the first IS110-family
  bridge recombinase to demonstrate >6% donor insertion efficiency in human cells
  (Pelea et al. 2026 Science adz1884; Perry et al. 2026 Science adz0276). n_editors: 30 -> 31.
  Final PenScore = 0.9181 (8-axis v0.1.1); S_Spec/S_Immuno use IS621 proxies;
  S_Mature=0.0 (0 PubMed hits 2026-05-22); S_Energy=1.0; MHCflurry 2.2.1 + BWA-MEM pending.

### Fixed

**`loader.py` EditorEntry missing `walker_motif_override` field (critical - S_Energy CAST override)**
- `pen_score/data/loader.py`: `walker_motif_override: bool | None = None` added to `EditorEntry`
  Pydantic model. Without this field, `getattr(ed, "walker_motif_override", None)` silently
  returned `None` for all editors, causing CAST systems (evoCAST, CAST_VK, CAST_IF) to skip
  the override and fall through to a sequence scan - which finds no Walker motifs in the
  Cas12k/TnsB primary accession -> incorrectly S_Energy=1.0 instead of the correct 0.0.
- Impact: evoCAST, CAST_VK, CAST_IF S_Energy 1.0 (wrong) -> 0.0 (correct). PenScore
  impact: -0.05 x weight for each affected editor.

**mech-class v0.5.2 IS110 Tier-A gate (critical - S_DSB correction for IS110-family)**
- `pen_score/axes/dsb.py`: handle `tier_a_gate_override=True` returned by mech-class v0.5.2;
  return `1.0` immediately when gate fires (IS110 domain-evidence confidence 0.99)
- `pyproject.toml`: bump mech-class pin `>=0.5.0` -> `>=0.5.2,<0.7.0`
- `pen_score/data/axis_definitions.yaml`: dependency updated to `mech-class>=0.5.2`
- S_DSB for IS110-family editors (IS621, IS621_2, IS622) is set to 1.0 when the Tier-A gate fires, up from the prior 0.90.
- Impact (7-axis): IS621 PenScore 0.9290 -> 0.9540 (+0.025); IS621_2 0.9000 -> 0.9250.
- Impact (8-axis v0.1.1): IS621 PenScore 0.9570; IS621_2 0.9280; SpCas9 0.4017.
  All pre-registered predictions (P1-P5) remain PASS / NOT_EVALUABLE.

---

## [0.1.0] - 2026-05-14

First complete public release.  All 7 score axes computed; pre-registered
predictions evaluated; public scorecard and interactive browser shipped.
Full documentation and outreach materials generated.

### Added

**Editor Universe**
- `pen_score/data/editor_universe.yaml` v1.0.4: 30 curated editors with
  verified UniProt accessions, mechanism buckets, cargo capacities, and
  pre-registered target flags
- `scripts/01_curate_editor_universe.py`, `02_validate_accessions.py`

**Axis Computation**
- `scripts/10_compute_S_DSB.py`: mech-class v0.5.1 aligned; bucket heuristic
  DSB_NUCLEASE->0.0, DSB_FREE->0.9, TRANSPOSASE->0.5
- `scripts/11_compute_S_Spec.py` v2: specificity_bias_factor correction for
  protein-interface variants (SpuFz1_V4: +0.05 post-sigmoid)
- `scripts/12_compute_S_Cargo.py`, `13_compute_S_Deliv.py`,
  `14_compute_S_Immuno.py`, `15_compute_S_Prog.py`, `16_compute_S_Mature.py`
- `scripts/20_assemble_scorecard.py`: public_scorecard.parquet (29 editors, 8 axes)
- `scripts/21_inter_axis_correlation.py`: correlation audit; 4 mechanistic
  pairs documented in EXPECTED_HIGH dict; correlation_audit.md generated
- `scripts/22_bootstrap_rankings.py`: 1000-iter bootstrap (seed=42, sigma=0.02);
  IS621 CI=[1,1], IS621_2=[2,2], evoCAST=[3,3]
- `SCORE_PROVENANCE.md` Section 11: mech-class v0.5.1 upstream dependency lock
  and impact table (13 DSB_NUCLEASE editors each -0.025 PenScore); v0.5.2 IS110
  Tier-A gate correction (IS621/IS621_2 S_DSB 0.90->1.00, PenScore +0.025 each)

**Pre-registration lock**
- Pre-registration tag `pre-registration-v1.0.2` (2026-05-13T17:56:38Z)
- `pen_score/data/editor_universe.yaml` v1.0.4 hash locked
- `specificity_bias_factor` for SpuFz1_V4 formalised and frozen

**Pre-Registered Predictions**
- `scripts/30_test_pred_1_evocast.py`: P1 PASS - evoCAST rank 5/13 DSB-free (8-axis)
- `scripts/31_test_pred_2_is621.py`: P2 PASS - IS621 rank 1/7 programmable DSB-free (8-axis)
- `scripts/32_test_pred_3_cas9.py`: P3 PASS - SpCas9 31.0% below threshold (>= 30%); rank 20/29
- `scripts/33_test_pred_4_ennlovfz2.py`: P4 NOT EVALUABLE - both editors sentinels
  (no public sequence; resolution path documented)
- `scripts/34_test_pred_5_spufz1_v4.py`: P5 PASS - SpuFz1_V4 1.0000 > SpuFz1 0.9999
- `scripts/35_summarise_predictions.py`: P1-P5 summary table + policy output
  (`predictions/prediction_summary.json`, `.csv`)

**Final Documentation**
- `docs/scorecards/index.html`: self-contained interactive scorecard browser
  (24 KB, zero dependencies); 5 use-case profiles, sort/search/filter, IS621 #1
- `docs/conf.py`: `html_extra_path = ["scorecards"]` for Sphinx static serving
- `pen_score/cli.py`: `select` command with `--require-dsb-free` flag; Rich table
  with per-axis color coding and strength/weakness reasoning column
- `pen_score/api.py`: `select_editor(require_dsb_free=False)` parameter;
  `_generate_reasoning()` method; `ScoringResult.reasoning` field populated
- `README.md`: CLI Selection Examples section (5 use cases with bash commands);
  Interactive Scorecard Browser section; updated predictions table
- `scripts/41_generate_outreach_materials.py`: generates 5 wet-lab outreach
  Markdown summaries (Hsu/Arc, Liu/Broad, Zhang/Broad, Sternberg/Columbia,
  SJNAHS/VIT); each with editor table, top-5 recommendations, axis reasoning,
  PEN-ASSEMBLE collaboration proposal

**Release & Reproducibility Audit**
- `tests/unit/test_axes.py`: completed TestComposite class (5 tests)
- `tests/unit/test_api.py`: Scorer._composite, _generate_reasoning, select_editor
  with require_dsb_free (12 tests)
- `tests/regression/test_scorecard_regression.py`: golden value regression tests
  locked to 8-axis v0.1.1 (IS621 0.957, SpCas9 0.402, 29 editors, IS621 #1
  across use cases, P3 bottom-30% criterion, DSB_NUCLEASE all S_DSB=0.0)
- `tests/test_placeholder.py`: updated version check to 0.1.0; added Scorer
  and CLI import smoke tests
- `pen_score/_version.py`: bumped to 0.1.0
- `CITATION.cff`: version 0.1.0, date-released 2026-05-14 (updated to 0.1.1 in v0.1.1)

### Fixed

**mech-class v0.5.1 composite gate (critical)**
- `scripts/10_compute_S_DSB.py` `_bucket_heuristic()`: DSB_NUCLEASE 0.1 -> 0.0
  (composite bonus now requires PF01548 AND PF02371 simultaneously; SpCas9
  RuvC+HNH does not trigger IS110-class bonus)
- Impact: 13 DSB_NUCLEASE editors each -0.025 PenScore; scorecard and
  bootstrap CIs rebuilt from scratch

**int64 JSON serialisation**
- `scripts/35_summarise_predictions.py`: wrapped `.sum()` calls with `int()`
  to fix numpy int64 JSON serialisation error

---

## [0.0.2] - 2026-04-22

### Fixed
- `pyproject.toml`: placeholder `pen-core` dependency removed; correct optional
  extras pattern restored

---

## [0.0.1] - 2026-05-12

### Added

- Package scaffold: complete package structure.
- `pen_score/` package with `api.py`, `cli.py`, `axes/`, `scorer/`, `utils/`, `data/`
- `editor_universe.yaml` v1.0.0: 30 curated programmable genome editors
- `use_case_profiles.yaml`: 5 use-case weight profiles
- `cargo_capacity_table.yaml`: literature-curated cargo capacities
- `pyproject.toml`, `SCORE_PROVENANCE.md`, `MODEL_CARD.md`, `VALIDATION.md`,
  `UPDATE_STRATEGY.md`, `CITATION.cff`, `CHANGELOG.md`
- `containers/biophysics/Dockerfile`
- Tests scaffold: `tests/unit/`, `tests/integration/`, `tests/regression/`
- GitHub Actions CI with quality + test + release pipeline

### Fixed

- `pen_score/__init__.py`, `cli.py`: removed template placeholders
- `pyproject.toml`: author email and dependency corrections
- `.github/workflows/ci.yml`: correct package references and secrets
