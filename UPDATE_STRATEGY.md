# UPDATE_STRATEGY.md - PEN-SCORE

**Package:** pen-score
**Author:** Anees Ahmed Mahaboob Ali
**Analogous to:** UPDATE_STRATEGY.md in genome-atlas and mech-class
**Last updated:** 2026-05-12

---

## Version Roadmap

| Version | Milestone | Key changes |
|---|---|---|
| v0.0.1 | Scaffold | Package skeleton; PyPI name reservation |
| v0.5.0 | PEN-SCORE public release | Full 7-axis pipeline; public scorecard; 30 editors |
| v0.6.0 | Post-review | Fusion editor S_Deliv fix; SaProt S_Immuno; Tier B S_Prog |
| v1.0.0 | Stable | PyPI public release |

---

## When to Update the Editor Universe

**v1.x minor updates** (patch, backward-compatible):
- New editors added to `editor_universe.yaml` (n_editors++, version bump)
- Accession corrections (follow MECH-CLASS discipline: document in SCORE_PROVENANCE.md)
- PubMed count refresh (S_Mature values stale within ~6 months)

**v0.6 changes** (breaking for S_Deliv scores):
- Fusion editor lengths: replace parent-accession length with full fusion protein size
  for PE2, ABE, BE3, TwinPE, PE5max once a standard fusion-length lookup is available

**Not updated without re-computation:**
- Score axis formulas (any formula change invalidates all cached scores -> requires
  full re-run of affected axis scripts)
- Axis weight profiles (use-case profiles may be updated but re-ranking is needed)

---

## Deprecation Policy

- Scores computed with v0.5 axis formulas are tagged with `pen_score_version: "0.5.0"`
  in the scorecard Parquet.
- Downstream tools (PEN-COMPARE) must check this tag and warn if stale.
- No backward-compatibility shims: if the formula changes, re-run the scripts.

---

## Raw Data

The raw data files that accompany a release are listed below.
- `editor_universe.yaml`, all axis parquets, `scorecard.parquet`,
  `bootstrap_rankings.parquet`, `holdout_results.json`, `MODEL_CARD.md`,
  `SCORE_PROVENANCE.md`, `VALIDATION.md`
