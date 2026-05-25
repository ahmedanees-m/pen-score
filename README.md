<div align="center">

# PEN-SCORE

### *Programmable Enzyme Networks - Systematic Comparative Output for Ranking Editors*

[![CI](https://github.com/ahmedanees-m/pen-score/actions/workflows/ci.yml/badge.svg)](https://github.com/ahmedanees-m/pen-score/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ahmedanees-m/pen-score/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/pen-score)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://ahmedanees-m.github.io/pen-score/)
[![Release](https://img.shields.io/github/v/release/ahmedanees-m/pen-score)](https://github.com/ahmedanees-m/pen-score/releases)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.3-green)](CHANGELOG.md)

**Upstream:**
[![GENOME-ATLAS](https://img.shields.io/badge/GENOME--ATLAS-v0.7.2-blue)](https://github.com/ahmedanees-m/genome-atlas)
[![MECH-CLASS](https://img.shields.io/badge/MECH--CLASS-v0.5.4-blue)](https://github.com/ahmedanees-m/mech-class)

**Part of [PEN-STACK](https://github.com/ahmedanees-m)** - Programmable Enzyme Networks, Systematic Tool for Atlas and Knowledge

Multi-axis scoring framework for ranking programmable genome editors by therapeutic potential

</div>

---

## The Problem PEN-SCORE Solves

The genome editing field has exploded. As of 2026, wet-lab researchers face more than **29 distinct editor classes** - CRISPR-Cas nucleases, prime editors, base editors, CAST transposases, IS110 bridge recombinases, Fanzors, site-specific recombinases, and more. Each has a different safety profile, cargo limit, delivery route, immunogenicity, and maturity level.

**Selecting the right editor for a clinical application is not obvious, and mistakes are costly.** A team that chooses SpCas9 for AAV gene therapy because it is the most published editor may overlook that its 997-amino-acid size makes single-AAV delivery borderline, it creates double-strand breaks that can cause chromosomal translocations, and its immunogenicity is among the highest of any editor tested. A team that reaches for a newer IS110 bridge recombinase because of impressive recent results may not realise that it lacks clinical safety data entirely.

No principled multi-axis comparison framework existed. Publications compare two or three editors on one or two properties. Reviews are qualitative. Without a structured, reproducible way to weigh all relevant properties simultaneously - and to tune those weights to a *specific clinical scenario* - researchers were making consequential decisions on incomplete information.

**PEN-SCORE was built to fill that gap.** It is a computational triage tool: score all 29 editors on eight biologically motivated axes, combine them into a single weighted composite score tuned to your clinical application, and produce a ranked shortlist in seconds. Use it to narrow 29 candidates to 3-5 worth characterising in the lab.

---

## Key Result

> **IS621** (an IS110-family bridge recombinase) ranks **#1 out of 29 editors** with a PenScore of **0.9570**. SpCas9 - the most widely used gene editor in the world - ranks **20th** with a PenScore of **0.4017**. The difference is driven almost entirely by mechanism class: IS621 edits DNA without any strand breaks (S_DSB = 1.0); SpCas9 cuts both strands (S_DSB = 0.0), which carries the highest axis weight (0.24) in the AAV therapeutic profile.

The top 13 editors are all DSB-free. This clean biological stratification was not imposed by the weights - it emerged from the data.

---

## How PEN-SCORE Works

PEN-SCORE is a three-stage pipeline that runs on top of two upstream tools in the PEN-STACK ecosystem.

### Stage 1 - Curate the editor universe

A YAML-defined universe of 29 scored editors (plus 2 sentinels) specifies each editor's UniProt canonical accession, mechanism override flags, literature cargo values, and PubMed query terms. All 29 accessions were verified against the GENOME-ATLAS knowledge graph; 27 corrections were made during curation (wrong organism, isoform confusion, obsolete IDs). The curated universe is the single source of truth that all downstream scripts consume.

### Stage 2 - Compute eight axes independently

Each axis is computed from first principles by a dedicated script and a corresponding module in `pen_score/axes/`. Axes are **completely independent** - one axis failing or being unavailable does not affect the others.

```
Editor accession (UniProt)
         │
         ├──→ S_DSB   ──── MECH-CLASS v0.5.4 Tier-A gate
         │                  └─ DSB_NUCLEASE → 0.0  |  DSB_FREE → 1.0  |  TRANSPOSASE → 0.5
         │
         ├──→ S_Spec  ──── CRISPOR + BWA-MEM off-target scan (GRCh38, 20-bp protospacer)
         │                  └─ sigmoid(off-target count) → [0, 1]
         │
         ├──→ S_Cargo ──── Literature-curated maximum payload capacity
         │                  └─ log-sigmoid(cargo_kb) → [0, 1]
         │
         ├──→ S_Deliv ──── UniProt sequence length lookup
         │                  └─ sigmoid(length_aa, centre=900) → [0, 1]
         │
         ├──→ S_Immuno ─── netMHCpan-4.1 + netMHCIIpan-4.0 epitope load
         │                  └─ combined_epitopes / 95th_percentile → [0, 1]  (inverted: fewer = better)
         │
         ├──→ S_Prog  ──── MECH-CLASS v0.5.4 Tier-B programmability class
         │                  └─ RNA-guided → 1.0  |  fixed-specificity → 0.0
         │
         ├──→ S_Mature ─── NCBI E-utilities PubMed citation count
         │                  └─ log(citations + 1) / log(max_citations + 1) → [0, 1]
         │
         └──→ S_Energy ─── Walker A (G[A-Z]{4}GK[ST]) + Walker B ([LVIMF]{4}DE) regex scan
                            └─ motif absent → 1.0 (ATP-free)  |  motif present → 0.0 (ATP-dependent)
```

### Stage 3 - Compose a weighted PenScore

The eight axis scores are combined into a single **PenScore** using a weighted average, where the weights come from a clinical use-case profile. Seven profiles are pre-registered - each assigns different importance to different axes depending on the therapeutic context.

```
PenScore = Σ(weight_i x score_i) / Σ(weight_i for axes with data)
```

If an axis score is unavailable (e.g. a sentinel editor with no public sequence), that axis's weight is redistributed proportionally to the axes that *do* have data. This means a missing axis deflates the score by absence-of-evidence rather than treating it as zero.

The `exclude_axes` parameter allows you to drop an axis entirely from the composite - useful for newly characterised editors where `S_Mature = 0.0` reflects publication recency, not a biophysical limitation.

---

## Architecture - PEN-SCORE in PEN-STACK

PEN-SCORE is the third tool in the PEN-STACK pipeline. It consumes outputs from GENOME-ATLAS and MECH-CLASS, and feeds its scorecard into PEN-COMPARE (the web interface).

```
PEN-STACK pipeline
══════════════════════════════════════════════════════════════════════════════

  ┌──────────────────────────────────────────────────────────────────────────┐
  │  GENOME-ATLAS                                                            │
  │  github.com/ahmedanees-m/genome-atlas · v0.7.2                           │
  │                                                                          │
  │  Heterogeneous GNN knowledge graph                                       │
  │  13,401 nodes · 11,817 edges · 28 foundational systems                   │
  │  UniProt + PDB + Pfam + ESM-2 embeddings                                 │
  │                                                                          │
  │  Provides PEN-SCORE with:                                                │
  │    • Verified UniProt accessions for all 29 editors (27 corrections)     │
  │    • ISCro4 / D2TGM5 canonical node (formerly IS622, Pelea 2026 Science) │
  │    • SIMILAR_TO / HAS_RNA / PART_OF edges via graph_view='full'          │
  └────────────────────────────┬─────────────────────────────────────────────┘
                               │  accession lookup · knowledge graph edges
                               ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │  MECH-CLASS                                                              │
  │  github.com/ahmedanees-m/mech-class · v0.5.4                             │
  │                                                                          │
  │  Two-tier LightGBM classifier · 572-protein gold set                     │
  │  ESM-2 (640-dim) + Pfam domain flags (26-dim) + SaProt structural feats  │
  │  Tier-A macro-F1 = 0.9862                                                │
  │                                                                          │
  │  Provides PEN-SCORE with:                                                │
  │    • Mechanism class → S_DSB score                                       │
  │      v0.5.2 IS110 Tier-A gate: IS621 = 1.0 (corrects OOD mis-scoring)    │
  │    • Programmability tier (Tier-B) → S_Prog score                        │
  └────────────────────────────┬─────────────────────────────────────────────┘
                               │  mechanism class · programmability tier
                               ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │  PEN-SCORE  ← YOU ARE HERE · v0.1.3                                      │
  │  github.com/ahmedanees-m/pen-score                                       │
  │                                                                          │
  │  8-axis scoring framework · 29 editors · 7 clinical use-case profiles    │
  │  Pre-registered predictions: 4 / 4 evaluable PASS                        │
  │  Coverage: 93 % · 267 tests pass                                         │
  │                                                                          │
  │  Outputs:                                                                │
  │    • public_scorecard.parquet  (29 editors x 8 axes + PenScore)          │
  │    • bootstrap_rankings.parquet  (95 % CIs, 1000 iterations, seed = 42)  │
  │    • Interactive scorecard browser  (docs/scorecards/index.html)         │
  └────────────────────────────┬─────────────────────────────────────────────┘
                               │  ranked scorecard · PenScore per editor
                               ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │  PEN-ASSEMBLE                                                            │
  │  Chimeric editor design from top-ranked components                       │
  └────────────────────────────┬─────────────────────────────────────────────┘
                               │
                               ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │  PEN-COMPARE                                                             │
  │  Web interface - browse, compare, and filter all editors by axis         │
  │  PEN-COMPARE v3.2 uses get_editor_metadata() API from pen-score v0.1.3   │
  └──────────────────────────────────────────────────────────────────────────┘
```

---

## The Eight Score Axes

Each axis measures one biologically important property on a **[0, 1]** scale (1 = better for therapeutic use). The composite PenScore is a weighted average across whichever axes have data - missing axes are handled by weight redistribution, not zero-imputation.

| Axis | Symbol | Score of 1.0 means... | Score of 0.0 means... | Source |
|------|--------|---------------------|---------------------|--------|
| **DSB Safety** | `S_DSB` | Edits without cutting both DNA strands - zero double-strand break risk | Cuts both strands (SpCas9, Cas12a) - DSBs can cause chromosomal translocations | MECH-CLASS v0.5.4 Tier-A gate |
| **Specificity** | `S_Spec` | Hits only the intended target - no detectable off-target edits across the genome | High off-target activity at many unintended genomic sites | CRISPOR + BWA-MEM (GRCh38, 20 bp) |
| **Cargo Capacity** | `S_Cargo` | Can carry a very large genetic payload (> 10 kb) | Tiny payload capacity (< 100 bp) | Literature-curated, log-sigmoid |
| **Deliverability** | `S_Deliv` | Small enough for single-AAV delivery (< 4.5 kb) | Too large for single-AAV - needs split-intein or LNP | Sigmoid centred at 900 aa |
| **Immunogenicity** | `S_Immuno` | Few MHC epitopes - unlikely to trigger T-cell or antibody responses | Highly immunogenic - many predicted MHC-I and MHC-II epitopes | netMHCpan-4.1 + netMHCIIpan-4.0 |
| **Programmability** | `S_Prog` | RNA-guided - retarget to any sequence by changing a guide RNA | Fixed specificity - only works at its natural att-site or recognition sequence | MECH-CLASS Tier-B |
| **Tech Maturity** | `S_Mature` | Many published clinical studies - established safety and efficacy data | Characterised in 2025-2026 with limited published data | PubMed citation count, log-normalised |
| **Energy Independence** | `S_Energy` | No ATP hydrolysis required - simpler biochemistry, no energy co-factor | Requires ATP - adds a dependency on cellular energy | Walker A/B motif regex scan (v0.1.1) |

### How MECH-CLASS feeds S_DSB and S_Prog

MECH-CLASS classifies each protein into one of three mechanism tiers. This is the most important axis - it drives the clean stratification between the top 13 (all DSB-free) and the rest of the universe.

```
MECH-CLASS Tier-A output          →  S_DSB
──────────────────────────────────────────
DSB_NUCLEASE                      →  0.0   (SpCas9, Cas12a, ZFN, TALEN)
DSB_FREE_TRANSEST_RECOMBINASE     →  1.0   (IS621, IS621_2, ISCro4, IscB, Lambda-Int)
TRANSPOSASE                       →  0.5   (Sleeping Beauty, Tc1/mariner)
IS110 Tier-A gate override        →  1.0   (bypasses LightGBM for novel IS110 proteins)
```

The **v0.5.2 IS110 Tier-A hard gate** was critical for correctness. Novel IS110 proteins (IS621, ISCro4) had no pre-computed ESM-2 embeddings at inference time, causing the LightGBM model to return `DSB_NUCLEASE` with P = 0.567-0.703 due to out-of-distribution feature vectors. The gate detects the dual PF01548 + PF02371 Pfam domain signature and returns `DSB_FREE` with confidence >= 0.90 regardless of ML output. Without this fix, IS621 S_DSB would be 0.90 instead of 1.0 and its PenScore 0.9540 instead of 0.9570.

---

## Public Scorecard - Top 10 (Human Therapeutic AAV Insertion)

*IS621 ranks #1. All top-13 editors are DSB-free. SpCas9 ranks 20th out of 29.*

| Rank | Editor | PenScore | Class | S_DSB | S_Spec | S_Cargo | S_Deliv | S_Prog | S_Energy |
|------|--------|----------|-------|-------|--------|---------|---------|--------|----------|
| 1 | **IS621** | **0.9570** | Bridge recombinase | 1.00 | 1.0000 | 1.0000 | 0.9421 | 1.0 | 1.0 |
| 2 | IS621_2 | 0.9280 | Bridge recombinase | 1.00 | 1.0000 | 1.0000 | 0.9421 | 1.0 | 1.0 |
| 3 | ISCro4 † | 0.9181 | Bridge recombinase | 1.00 | 1.0000 | 1.0000 | 0.9463 | 1.0 | 1.0 |
| 4 | Lambda_Int | 0.8295 | Site-specific recombinase | 0.90 | 0.9880 | 0.7832 | 0.9382 | 0.0 | 1.0 |
| 5 | evoCAST * | 0.8128 | CAST transposase | 0.90 | 1.0000 | 0.7832 | - | 1.0 | 0.0 |
| 6 | eePASSIGE_v2 | 0.8045 | Serine integrase | 0.90 | 0.9999 | 0.8835 | 0.8808 | 0.0 | 1.0 |
| 7 | Cre | 0.8058 | Tyrosine recombinase | 0.90 | 0.8866 | 0.6667 | 0.9419 | 0.0 | 1.0 |
| 8 | Bxb1 | 0.7960 | Serine integrase | 0.90 | 0.9999 | 0.7832 | 0.8808 | 0.0 | 1.0 |
| 9 | IscB | 0.7939 | OMEGA/Fanzor nuclease | 0.90 | 1.0000 | 0.3835 | 0.9241 | 1.0 | 1.0 |
| 10 | phiC31 | 0.7940 | Serine integrase | 0.90 | 0.9995 | 0.8333 | 0.8138 | 0.0 | 1.0 |

*\* Provisional - sentinel editor with no public sequence for S_Deliv/S_Immuno*
*† Post-pipeline comparator (Pelea et al. 2026 Science); S_Mature = 0.0 (no prior clinical citations)*

Full 29-editor scorecard: [`scorecards/public_scorecard.parquet`](scorecards/)
Interactive browser: [`docs/scorecards/index.html`](docs/scorecards/index.html) - open locally in any browser, no server needed.

---

## Install

```bash
pip install pen-score
```

For full axis computation, install the relevant extras:

```bash
pip install "pen-score[mech-class]"   # S_DSB, S_Prog  (requires mech-class >= 0.5.4)
pip install "pen-score[spec]"         # S_Spec  (requires BWA + pysam)
pip install "pen-score[immuno]"       # S_Immuno  (requires MHCflurry or netMHCpan)
pip install "pen-score[ml]"           # LightGBM + scikit-learn for ML axes
pip install "pen-score[dev]"          # All of the above + pytest, ruff, mypy, sphinx
```

---

## Quick Start

### Score a single editor

```python
from pen_score.api import Scorer

scorer = Scorer.load()

# Score IS621 bridge recombinase for AAV gene therapy
result = scorer.score_editor(
    accession="A0A2X3M8B0",          # IS621 UniProt accession
    use_case="human_therapeutic_aav_insertion",
)

print(f"PenScore:  {result.pen_score}")      # 0.9570 - rank 1/29
print(f"S_DSB:     {result.axes.S_DSB}")    # 1.0 - no double-strand breaks
print(f"S_Cargo:   {result.axes.S_Cargo}")  # 1.0 - megabase-scale capacity
print(f"S_Energy:  {result.axes.S_Energy}") # 1.0 - no Walker A/B motifs (ATP-free)
print(f"S_Prog:    {result.axes.S_Prog}")   # 1.0 - RNA-guided, reprogrammable

# Axis-level strengths and weaknesses
for bullet in result.reasoning:
    print(bullet)
# [strength] DSB safety:       1.0000 (weight=0.24) - IS110 Tier-A gate; no double-strand breaks
# [strength] Cargo capacity:   1.0000 (weight=0.19) - megabase-scale payload
# [strength] Deliverability:   0.9421 (weight=0.19) - 326 aa; single-AAV compatible
# [weakness] Tech maturity:    0.0000 (weight=0.05) - new 2025 enzyme; limited clinical data
```

### Compare SpCas9 (the most widely used editor)

```python
# SpCas9 - included to illustrate what PEN-SCORE penalises
cas9 = scorer.score_editor("Q99ZW2", use_case="human_therapeutic_aav_insertion")
print(f"SpCas9 PenScore: {cas9.pen_score}")   # 0.4017 - rank 20/29
print(f"SpCas9 S_DSB:    {cas9.axes.S_DSB}")  # 0.0    - cuts both DNA strands
# SpCas9 is penalised for DSB risk (S_DSB=0.0), large size (997 aa), and high immunogenicity
```

### Query PEN-COMPARE v3.2 metadata (new in v0.1.3)

```python
from pen_score import get_editor_metadata

# ISCro4 (formerly IS622) - TRUE_WRITER tier in PEN-COMPARE v3.2
meta = get_editor_metadata("ISCro4")
print(meta.intrinsic_cargo_mechanism)  # True  - integrates cargo via its own catalytic mechanism
print(meta.cell_based_evidence)        # True  - peer-reviewed mammalian cell data > 1% efficiency
print(meta.cell_based_sources)         # ['Perry 2026 Science doi:10.1126/science.adz0276', ...]

# IS621 - PROBABLE_WRITER tier (intrinsic mechanism, but E. coli data only so far)
meta2 = get_editor_metadata("IS621")
print(meta2.cell_based_evidence)       # False - cryo-EM + E. coli; no robust human-cell data yet

# Backward-compatible alias resolution (IS622 was the preprint label for ISCro4)
meta3 = get_editor_metadata("IS622")   # DeprecationWarning: use 'ISCro4' (canonical name)
```

### Select the best editors for your use case

```python
# Top-5 DSB-free editors for AAV gene therapy
top5 = scorer.select_editor(
    use_case="human_therapeutic_aav_insertion",
    top_k=5,
    require_dsb_free=True,
)
print(top5[["editor_id", "PenScore", "S_DSB", "S_Cargo", "S_Prog"]])

# Top editors for megabase-scale rearrangements
mega = scorer.select_editor(use_case="megabase_rearrangement", top_k=5)

# BCL11A enhancer excision for sickle cell disease
scd = scorer.select_editor(use_case="therapeutic_excision_bcl11a", top_k=5)
```

### Biophysical-only comparison (exclude S_Mature)

```python
# ISCro4 (formerly IS622) was characterised in 2026 - S_Mature = 0.0 because
# there are no prior PubMed hits. This is publication recency, not a limitation.
# exclude_axes redistributes the weight to the remaining 7 axes automatically.
bio_result = scorer.score_editor(
    accession="D2TGM5",                # ISCro4 UniProt accession
    use_case="megabase_rearrangement",
    exclude_axes=["S_Mature"],         # remove recency penalty for new enzymes
)
print(f"ISCro4 biophysical PenScore: {bio_result.pen_score}")  # ~0.94
```

### CLI

```bash
# Score a single editor and display all axes
pen-score score-editor A0A2X3M8B0 --use-case human_therapeutic_aav_insertion

# Score ISCro4 excluding S_Mature (biophysical comparison for new enzymes)
pen-score score-editor D2TGM5 --use-case megabase_rearrangement --exclude-axes S_Mature

# Build the full 29-editor scorecard
pen-score score-all-editors --output scorecard.parquet

# Select top editors for each clinical profile
pen-score select --use-case human_therapeutic_aav_insertion   --top-k 5 --require-dsb-free
pen-score select --use-case human_therapeutic_electroporation --top-k 5
pen-score select --use-case large_cargo_integration           --top-k 5 --require-dsb-free
pen-score select --use-case base_editing_small_correction     --top-k 5
pen-score select --use-case research_discovery                --top-k 10
pen-score select --use-case megabase_rearrangement            --top-k 5 --require-dsb-free
pen-score select --use-case therapeutic_excision_bcl11a       --top-k 5 --require-dsb-free
```

---

## Use-Case Profiles

Each clinical context weights the eight axes differently. Switching profiles re-ranks the editors without changing any axis score - only the weights change.

| Profile key | Clinical scenario | Highest-weighted axes | Reference editor |
|-------------|-------------------|-----------------------|-----------------|
| `human_therapeutic_aav_insertion` | AAV-compatible gene addition for rare disease | S_DSB (0.24), S_Cargo (0.19), S_Deliv (0.19) | IS621 |
| `human_therapeutic_electroporation` | Ex vivo HSC / T-cell editing | S_DSB (0.24), S_Spec (0.20), S_Immuno (0.15) | IS621 |
| `large_cargo_integration` | Payload > 4.5 kb (lentiviral / LNP) | S_Cargo (0.30), S_DSB (0.25) | IS621 |
| `base_editing_small_correction` | Single-base correction (point mutations) | S_Spec (0.30), S_Immuno (0.20) | IS621 |
| `research_discovery` | Academic discovery; maturity less critical | S_Prog (0.30), S_Spec (0.20) | IS621 |
| `megabase_rearrangement` | > 100 kb inversions or translocations | S_Prog (0.30), S_Cargo (0.25), S_DSB (0.20) | ISCro4 |
| `therapeutic_excision_bcl11a` | BCL11A +58 enhancer excision for SCD/β-thal | S_DSB (0.25), S_Spec (0.25), S_Immuno (0.20) | IS621 |

Weight files: [`pen_score/data/use_case_profiles.yaml`](pen_score/data/use_case_profiles.yaml).

---

## Pre-Registered Validation

All five predictions were registered **before any score was computed** (tag `pre-registration-v1.0.2`, locked 2026-05-13T17:56:38 UTC, commit `af184cc`). This lock ensures results are not reverse-engineered.

| ID | Prediction | Result | Observed |
|----|-----------|--------|----------|
| P1 | evoCAST in top-5 DSB-free editors (AAV profile) | **PASS** | Rank 5/13 DSB-free (threshold: top-5) |
| P2 | IS621 in top-3 programmable DSB-free editors | **PASS** | Rank 1/7; bootstrap CI = [1, 1] |
| P3 | SpCas9 in bottom 30% overall | **PASS** | 31.0 % of editors score lower; rank 20/29 |
| P4 | enNlovFz2 S_Deliv > NlovFz2 S_Deliv | **NOT EVALUABLE** | Both are sentinels - sequences not yet deposited |
| P5 | SpuFz1_V4 S_Spec > SpuFz1 S_Spec | **PASS** | 1.0000 vs 0.9999 (specificity-bias factor applied) |

**4 / 4 evaluable predictions PASS.** P4 is a documented data-availability boundary (sentinels listed in SCORE_PROVENANCE.md section 4.5 before the pre-registration lock), not a model failure.

Bootstrap: 1,000 iterations, seed = 42, σ = 0.02.

---

## Key Findings

1. **The top 13 editors are all DSB-free.** Clean stratification by mechanism class is the headline quantitative result. This was not imposed by the weights - it emerged from data. Editors 1-13 have S_DSB >= 0.9; editors 14-29 are transposases and nucleases.

2. **IS621 ranks #1 across all seven use cases** with a bootstrap CI of [1, 1] - meaning it ranks first in all 1,000 bootstrap resamples. Its advantage is simultaneous strength on four axes: no DSB (S_DSB = 1.0), megabase cargo (S_Cargo = 1.0), single-AAV compatible (S_Deliv = 0.94), ATP-free (S_Energy = 1.0).

3. **SpCas9 is penalised primarily for mechanism, not maturity.** SpCas9 scores 0.89 on S_Mature (very well established) and 0.9999 on S_Spec. Its PenScore collapses to 0.4017 because S_DSB = 0.0 - the highest-weight axis - and its 997-aa size makes single-AAV delivery marginal.

4. **ISCro4 (formerly IS622) is IS621's nearest neighbour.** UniProt D2TGM5, demonstrated 0.93 Mb inversions in HEK293T (Pelea et al. 2026 *Science*). PenScore 0.9181; its only weakness is S_Mature = 0.0 (published 2026, no prior clinical citations). Use `exclude_axes=["S_Mature"]` for a biophysical comparison (PenScore ~ 0.94).

5. **CAST systems are penalised for energy dependence.** evoCAST, CAST_VK, and CAST_IF all have S_Energy = 0.0 because TnsC ATPase activity is required. This suppresses their composite score despite strong DSB, programmability, and specificity axes.

---

## Repository Layout

```
pen-score/
│
├── pen_score/                       ← Installable Python package
│   ├── api.py                       ← Scorer, EditorMetadata, get_editor_metadata()
│   │                                   score_editor() · select_editor() · get_scorecard()
│   │                                   exclude_axes parameter for biophysical comparisons (v0.1.2)
│   │                                   get_editor_metadata() for PEN-COMPARE v3.2 (v0.1.3)
│   ├── cli.py                       ← Click CLI: score-editor / select / score-all-editors
│   ├── axes/
│   │   ├── dsb.py                   ← S_DSB   (MECH-CLASS v0.5.4 Tier-A gate)
│   │   ├── spec.py                  ← S_Spec  (CRISPOR / pysam / GRCh38)
│   │   ├── cargo.py                 ← S_Cargo (literature log-sigmoid)
│   │   ├── deliv.py                 ← S_Deliv (UniProt length sigmoid, centre 900 aa)
│   │   ├── immuno.py                ← S_Immuno (netMHCpan-4.1 / netMHCIIpan-4.0)
│   │   ├── prog.py                  ← S_Prog  (MECH-CLASS Tier-B)
│   │   ├── mature.py                ← S_Mature (PubMed citation count, log-normalised)
│   │   ├── energy.py                ← S_Energy (Walker A/B motif scan, v0.1.1)
│   │   └── d7_homolumo.py           ← Supplementary HOMO-LUMO gap (xTB; Docker/VM only)
│   ├── scorer/
│   │   ├── composite.py             ← compute_pen_score(): weighted avg + exclude_axes + missing
│   │   ├── bootstrap.py             ← bootstrap_ranking_ci(), bootstrap_axis_ci() (seed=42)
│   │   └── ranker.py                ← rank_editors() with DSB-free and use-case filters
│   ├── data/
│   │   ├── editor_universe.yaml     ← 29 editors + 2 sentinels; v1.0.7 (ISCro4 canonical)
│   │   ├── axis_definitions.yaml    ← Mathematical specs for all 8 axes (v1.0.1)
│   │   ├── use_case_profiles.yaml   ← 7 clinical weight profiles, all sum to 1.00 (v1.0.2)
│   │   └── pre_registration.yaml    ← P1-P5 predictions + SHA-256 lock (v1.0.2)
│   └── utils/
│       ├── pdb.py                   ← PDB / AlphaFold structure fetching helpers
│       ├── pubmed.py                ← NCBI E-utilities client for S_Mature citation counts
│       └── uniprot.py               ← UniProt REST client for accession validation
│
├── scripts/                         ← Numbered pipeline scripts; run in order
│   ├── 01_curate_editor_universe.py ← Editor curation and UniProt validation
│   ├── 10-17_compute_S_*.py         ← One script per axis
│   ├── 18_score_perry2025_panel.py  ← Perry 2025 IS110 ortholog panel (v0.1.2)
│   ├── 20_assemble_scorecard.py     ← Merge 8 axis parquets → public_scorecard.parquet
│   ├── 21_inter_axis_correlation.py ← 8x8 Spearman ρ matrix
│   ├── 22_bootstrap_rankings.py     ← 1000-iter bootstrap (seed=42, σ=0.02)
│   ├── 30-35_test_pred_*.py         ← P1-P5 pre-registered prediction evaluation
│   ├── 40_generate_html_browser.py  ← Self-contained scorecard browser (24 KB)
│   └── 41_generate_outreach_materials.py
│
├── tests/
│   ├── unit/                        ← 267 tests pass · 93 % coverage · 0 failures
│   │   ├── test_api.py              ← Scorer composite, reasoning, select_editor, exclude_axes
│   │   ├── test_axes.py             ← S_Cargo, S_Deliv, S_Prog, S_Energy (no external deps)
│   │   ├── test_axes_mocked.py      ← S_Spec / S_Immuno / S_Mature / S_DSB (mocked)
│   │   ├── test_bootstrap.py        ← bootstrap_ranking_ci, bootstrap_axis_ci
│   │   ├── test_cli_extended.py     ← Click CliRunner: all commands + flags (v0.1.3)
│   │   ├── test_d7_homolumo.py      ← ImportError path; xtb block is pragma: no cover
│   │   ├── test_energy.py           ← Walker A/B regex + CAST overrides (v0.1.1)
│   │   ├── test_get_editor_metadata.py ← EditorMetadata, alias resolution (v0.1.3)
│   │   ├── test_loader.py           ← editor universe + use-case profiles YAML loading
│   │   ├── test_pdb_utils.py        ← fetch_pdb_structure, fetch_alphafold_structure (mocked)
│   │   ├── test_profiles.py         ← 7 profiles: weight sums, semantic constraints
│   │   ├── test_ranker.py           ← rank_editors(), _build_reasoning() (v0.1.3)
│   │   ├── test_scorer_extended.py  ← Scorer.score_editor from scorecard, get_scorecard (v0.1.3)
│   │   └── test_utils.py            ← UniProt + PubMed REST clients (mocked)
│   ├── regression/
│   │   └── test_scorecard_regression.py ← Golden values: IS621 rank 1, SpCas9 S_DSB=0
│   └── test_placeholder.py          ← Version smoke test; Scorer + CLI importability
│
├── docs/                            ← Sphinx documentation → GitHub Pages
│   ├── scorecards/index.html        ← Interactive browser: sort · filter · switch use cases
│   └── ...
│
├── containers/
│   ├── biophysics/Dockerfile        ← xTB + netMHCpan + pysam environment (VM only)
│   └── spec/Dockerfile              ← CRISPOR + BWA + GRCh38 environment
│
├── .github/workflows/
│   ├── ci.yml                       ← ruff → tests + coverage → Codecov upload
│   └── docs.yml                     ← Sphinx → gh-pages branch
│
├── MODEL_CARD.md                    ← Capabilities, per-axis limitations, update conditions
├── SCORE_PROVENANCE.md              ← Full audit trail: data sources, formulas, edge cases
├── VALIDATION.md                    ← Pre-registration criteria, P1-P5 operational definitions
├── CHANGELOG.md                     ← Version history (v0.0.1 → v0.1.3)
├── pyproject.toml                   ← Package config, ruff/mypy/pytest settings
└── LICENSE                          ← MIT
```

---

## Score Axis Details

### S_DSB - DSB Safety (MECH-CLASS v0.5.4)

MECH-CLASS classifies the catalytic mechanism of each protein. The v0.5.2 IS110 Tier-A gate uses the dual PF01548 + PF02371 Pfam signature to hard-assign `DSB_FREE` status to IS110 bridge recombinases, overriding LightGBM when OOD inputs make it unreliable.

### S_Energy - Energy Independence (v0.1.1)

Walker A (`G[A-Z]{4}GK[ST]`) and Walker B (`[LVIMF]{4}DE`) motif scans on the canonical UniProt FASTA. Presence of either motif → S_Energy = 0.0 (ATP-dependent). Neither motif → 1.0. Multi-subunit CAST systems receive `walker_motif_override: true` in `editor_universe.yaml` because TnsC ATPase is a separate subunit not captured by the Cas12k/TnsB primary accession.

### S_Immuno - Immunogenicity

netMHCpan-4.1 predicts MHC class I epitopes across HLA-A02:01, HLA-A01:01, HLA-B07:02, HLA-B44:02 (European reference panel). netMHCIIpan-4.0 predicts class II across DRB1_0101, DRB1_0301, DRB1_0401. Combined epitope load is normalised against the 95th-percentile ceiling. **SpCas9 S_Immuno = 0.0 is correct biology**: n_I = 49, n_II = 409 epitopes; combined load saturates the normalisation ceiling.

---

## Reproducibility

### Requirements

```
Python       >= 3.10
mech-class   >= 0.5.4, < 0.6.0   # ISCro4 holdout probe rename
genome-atlas >= 0.7.2, < 0.8.0   # ISCro4 canonical in foundational_systems
pen-score    == 0.1.3
```

### Reproduce the scorecard from axis parquets (~10 seconds)

```bash
git clone https://github.com/ahmedanees-m/pen-score.git
cd pen-score
pip install ".[dev]"

python scripts/20_assemble_scorecard.py           # merge 8 axis parquets
python scripts/35_summarise_predictions.py        # evaluate P1-P5
```

### Re-run the full scoring pipeline from scratch

```bash
# Requires: mech-class >= 0.5.4, genome-atlas >= 0.7.2, BWA+GRCh38, netMHCpan-4.1
python scripts/10_compute_S_DSB.py          # ~5 min  - MECH-CLASS inference
python scripts/11_compute_S_Spec.py         # ~2 h    - BWA-MEM off-target scan
python scripts/12_compute_S_Cargo.py        # ~1 min  - literature values
python scripts/13_compute_S_Deliv.py        # ~1 min  - UniProt length lookup
python scripts/14_compute_S_Immuno.py       # ~3 h    - netMHCpan / netMHCIIpan
python scripts/15_compute_S_Prog.py         # ~5 min  - MECH-CLASS Tier-B
python scripts/16_compute_S_Mature.py       # ~5 min  - NCBI E-utilities
python scripts/17_compute_S_Energy.py       # ~3 min  - Walker A/B regex scan
python scripts/20_assemble_scorecard.py     # ~10 s   - merge all axes
python scripts/22_bootstrap_rankings.py     # ~2 min  - 1000-iter bootstrap, seed=42
```

All bootstrap operations use **seed = 42** (fixed per PEN-STACK convention).

---

## Upstream Dependencies

| Tool | Version | What PEN-SCORE uses |
|------|---------|----------------------|
| [GENOME-ATLAS](https://github.com/ahmedanees-m/genome-atlas) | v0.7.2 | UniProt accession verification · ISCro4/D2TGM5 canonical · SIMILAR_TO secondary edges |
| [MECH-CLASS](https://github.com/ahmedanees-m/mech-class) | v0.5.4 | Tier-A mechanism class → S_DSB · Tier-B programmability → S_Prog · IS110 Tier-A hard gate |

---

## Development

```bash
git clone https://github.com/ahmedanees-m/pen-score.git
cd pen-score
pip install ".[dev]"

pytest tests/unit/ tests/test_placeholder.py -v --cov=pen_score
ruff check pen_score/ tests/
mypy pen_score/ --ignore-missing-imports
```

---

## Documentation

- **Score provenance:** [`SCORE_PROVENANCE.md`](SCORE_PROVENANCE.md) - data sources, formulas, 27 accession corrections, edge cases for every axis
- **Model card:** [`MODEL_CARD.md`](MODEL_CARD.md) - intended use, per-axis limitations, inherited MECH-CLASS limitations
- **Interactive browser:** [`docs/scorecards/index.html`](docs/scorecards/index.html) - sort by any axis, filter by mechanism class, switch use cases; no internet required
- **Sphinx docs:** [https://ahmedanees-m.github.io/pen-score/](https://ahmedanees-m.github.io/pen-score/)

---

## Citation

If you use pen-score, please cite:

> Mahaboob Ali, A. A. (2026). PEN-SCORE: Multi-axis scoring framework for programmable genome editors. https://github.com/ahmedanees-m/pen-score

```bibtex
@software{pen_score_2026,
  author  = {Mahaboob Ali, Anees Ahmed},
  title   = {{PEN-SCORE}: Multi-axis scoring framework for programmable genome editors},
  version = {0.1.3},
  year    = {2026},
  url     = {https://github.com/ahmedanees-m/pen-score},
  license = {MIT}
}
```

Please also cite the upstream tools your analysis depends on:

```bibtex
@software{genome_atlas_2026,
  author  = {Mahaboob Ali, Anees Ahmed},
  title   = {{GENOME-ATLAS}: Heterogeneous GNN knowledge graph for genome editors},
  version = {0.7.2},
  year    = {2026},
  url     = {https://github.com/ahmedanees-m/genome-atlas}
}

@software{mech_class_2026,
  author  = {Mahaboob Ali, Anees Ahmed},
  title   = {{MECH-CLASS}: Mechanism classifier for programmable genome editors},
  version = {0.5.4},
  year    = {2026},
  url     = {https://github.com/ahmedanees-m/mech-class}
}
```

---

## License

MIT - see [LICENSE](LICENSE).

---

*Part of [PEN-STACK](https://github.com/ahmedanees-m) - Programmable Enzyme Networks, Systematic Tool for Atlas and Knowledge. Unified computational infrastructure for non-destructive genome engineering.*
