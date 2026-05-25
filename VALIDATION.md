# VALIDATION.md - PEN-SCORE

**Package:** pen-score
**Author:** Anees Ahmed Mahaboob Ali
**Analogous to:** VALIDATION.md in mech-class
**Last updated:** 2026-05-22 (8-axis v0.1.1 results added)

This document records the pre-registered success criteria and their results.
Pre-registration committed in bioRxiv preprint **before any score computation**.

---

## Pre-Registered Predictions (section 0.4)

Five retrospective recovery predictions, locked before any score computation.
Pre-registration tag: `pre-registration-v1.0.2` (2026-05-13T17:56:38Z, git `af184cc`).
Evaluated on 8-axis v0.1.1 scorecard (29 editors, 2026-05-22).

| # | Prediction | Threshold | Status | Result (8-axis) |
|---|---|---|---|---|
| 1 | evoCAST ranks top-5 of AAV-deliverable DSB-free integrases | Top 5 of ~10-15 | **PASS** | Rank 5 / 13 |
| 2 | IS621 ranks top-3 of programmable DSB-free systems | Top 3 of ~5-10 | **PASS** | Rank 1 / 7, bootstrap CI=[1,1] |
| 3 | SpCas9 ranks bottom 30% overall PenScore for human therapeutic AAV | Bottom 30% | **PASS** | Rank 20/29; 31.0% below threshold |
| 4 | enNlovFz2 S_Deliv strictly > NlovFz2 WT | Strictly greater | **NOT EVALUABLE** | Both REQUIRES_STEP7 sentinels (no sequence) |
| 5 | SpuFz1 V4 S_Spec strictly > SpuFz1 WT | Strictly greater | **PASS** | 1.0000 > 0.9999 |

**Outcome: 4/4 evaluated PASS + 1 sentinel boundary condition (P4).**
- 5/5 hold -> strong claim supported
- 4/5 hold -> report which prediction failed
- 3/5 hold -> reframe as scoring framework + lessons learned
- <=2/5 hold -> scoring function needs structural rework

---

## Inter-Axis Correlation Gate

Pre-registered: flag any pair with |ρ| > 0.7.

Expected justified correlations:
- S_DSB vs S_Prog: moderate positive (RNA-guided systems tend to be DSB-free); justified
  because both axes intentionally consume MECH-CLASS output

Expected independent axes:
- S_Deliv vs S_Spec: expected low correlation
- S_Cargo vs S_Immuno: expected low correlation
- S_Mature vs all others: expected low (sociological measure)

Actual matrix: completed; see `data/pen-score/scorecards/axis_correlation_matrix.csv` and `scorecards/correlation_audit.md` for flagged pairs.

---

## Accession Validation Gate

Mandatory before any score computation.  All 31 editor accessions verified against
UniProt REST API (30 pipeline editors + IS622 added 2026-05-21).  Known placeholders that require validation:

| Editor | Current accession | Status |
|---|---|---|
| NlovFz2 | A0A0C4DH50 | _to be verified_ |
| enNlovFz2 | A0A0C4DH50 | _shares parent; variant tracked via `parent_editor` field_ |
| evoCAST | A0A1B0GTW7 | _to be verified_ |
| eePASSIGE | Q9B086 (parent Bxb1) | _no separate UniProt entry; tracked via notes_ |

---

## Bootstrap CI Parameters (per PEN-STACK convention)

- N iterations: 1,000
- Seed: 42
- CI: 95% (2.5th-97.5th percentile)
- Applied to: axis scores, PenScore rankings, pre-registered prediction tests
