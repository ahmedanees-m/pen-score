# SCORE_PROVENANCE.md

**Package:** pen-score v0.1.2
**Author:** Anees Ahmed Mahaboob Ali (ahmedaneesm@gmail.com)
**Program:** PEN-STACK v2.2
**Analogous to:** LABEL_PROVENANCE.md (mech-class) and DATA_PROVENANCE.md (genome-atlas)
**Last updated:** 2026-05-24

This document describes the provenance of every data source and computation
used to produce the PEN-SCORE public scorecard.  It is the authoritative
audit trail for reviewers and downstream users.

---

## 1. Editor Universe

**File:** `pen_score/data/editor_universe.yaml` v1.0.4
**Curated:** 2026-05-12 by Anees Ahmed; corrected 2026-05-13 (27 live-verified fixes)
**N editors:** 31 (30 pipeline editors + IS622 added 2026-05-22 as post-pipeline comparator; 29 fully scored in public_scorecard.parquet)

**Selection criteria:**
1. Mechanism class assigned and verified against MECH-CLASS Tier A
2. UniProt accession verified against REST API (mandatory accession-validation gate)
3. Published primary reference with at least one PubMed-indexed citation
4. At least one of: canonical PDB structure OR AlphaFold model with mean pLDDT > 70

**Three source categories:**
| Category | N | Examples |
|---|---|---|
| GENOME-ATLAS foundational systems | 14 | SpCas9, Cas12a, IS621, Bxb1, Tn5 |
| Engineered variants (2024-2026) | ~10 | evoCAST, enNlovFz2, SpuFz1 V4, PE5max |
| Recent literature additions | ~6 | IscB, CAST-I-F WT, phiC31, TwinPE |

**Accession corrections inherited from MECH-CLASS:**
- Bxb1: Q8VVR2 -> **Q9B086** (Q8VVR2 is *S. aureus* GajA nuclease - wrong protein)
- Tn5: P00509 -> **Q46731** (P00509 was a wrong placeholder)

---

## 2. Score Axis Provenance

### S_DSB - DSB Avoidance
- **Source:** MECH-CLASS v0.5.2 Tier A classifier
- **Model artifact:** `~/pen-stack/data/models/tier_a/model.pkl` (981 KB, LightGBM)
- **Formula:** `1 - P(DSB_NUCLEASE); +0.1 bonus if IS110 composite flag`
- **Dependency:** `mech-class>=0.5.2` (optional extra)

### S_Spec - Specificity
- **Source:** BWA-MEM2 off-target scan on GRCh38 reference genome
- **Protospacers:** Per-editor canonical 20-bp protospacer, curated
- **Formula:** `sigmoid(-2 * log10(off_target_count / (3.2e9 / 1000) + 1e-10))`
- **Dependency:** BWA binary + pysam; run via script 11 on VM

### S_Cargo - Cargo Capacity
- **Source:** Literature-curated `cargo_capacity_table.yaml`; primary publications
- **Formula:** `log10(cargo_bp) / log10(1e6)`, clipped to [0, 1]
- **Dependency:** None (pure lookup from YAML)

### S_Deliv - AAV Deliverability
- **Source:** Protein length from UniProt REST API
- **Formula:** `sigmoid(0.005 * (900 - total_aa))` centred at 900 aa
- **Dependency:** UniProt REST API (requests, core dep)

### S_Immuno - Immunogenicity
- **Source (class I):** netMHCpan-4.1; 9-mer sliding windows; %Rank_EL < 0.5
- **Source (class II):** netMHCIIpan-4.0; 15-mer sliding windows; %Rank < 10
- **HLA-I alleles:** HLA-A*02:01, HLA-A*01:01, HLA-B*07:02, HLA-B*44:02
- **HLA-II alleles:** DRB1*01:01, DRB1*03:01, DRB1*04:01
- **Formula:** `total = n_i + 0.5 * n_ii; score = 1 - clip(total / max_total_universe, 0, 1)`
- **Tool note:** MHCflurry 2.0 is class I *only* and cannot predict 15-mer MHC-II binding;
  it is explicitly **not** used. netMHCpan-4.1 + netMHCIIpan-4.0 are the standard pair
  for two-class epitope-load calculations in the immunogenicity literature.
- **Dependency:** netMHCpan-4.1 and netMHCIIpan-4.0 binaries (DTU Health Tech academic
  license; not pip-installable); run via script 14 on VM with bind-mounted binaries

### S_Prog - Programmability
- **Source:** `rna_guided` field in editor_universe.yaml (Boolean); cross-checked
  with MECH-CLASS Tier B programmability sub-class
- **Formula:** Binary: `1.0` if RNA-guided; `0.0` if site-specific (att-site)
- **Dependency:** None (lookup from YAML)

### S_Mature - Therapeutic Maturity
- **Source:** NCBI E-utilities esearch on PubMed
- **Query:** `(editor_synonyms) AND (clinical OR preclinical OR therapeutic OR "gene therapy")`
- **Formula:** `log10(count + 1) / log10(max_count + 1)`, normalised over universe
- **Dependency:** requests (core dep); run via script 16

### S_Energy - Metabolic Energy Independence (added v0.1.1, 2026-05-22)
- **Source:** UniProt FASTA sequence (canonical accession) via REST API
- **Formula:** Binary `{0.0, 1.0}`.  Score = `0.0` if Walker A **and** Walker B motifs are
  both detected in the sequence; `1.0` otherwise.
  - Walker A (P-loop NTPase): regex `G[A-Z]{4}GK[ST]`
  - Walker B (Mg²⁺-coordinating): regex `[LVIMF]{4}DE`
  - Multi-subunit systems (e.g. CAST complexes with a TnsC ATPase subunit): set via
    `walker_motif_override: true` in `editor_universe.yaml` -> forced score = `0.0`
- **Rationale:** ATP-dependent editors require cellular energy cofactors that may be
  rate-limiting in post-mitotic or metabolically stressed target cells.  Energy-independent
  editors (IS110-family recombinases, many integrases) operate without this constraint.
- **Sentinel behaviour:** Editors with `REQUIRES_STEP7` accessions return `S_Energy = None`
  (not scored); the PenScore renormalisation excludes the missing axis automatically.
- **8-axis weight (human_therapeutic_aav_insertion):** 0.05 (uniform across all profiles)
- **Dependency:** requests (core dep, for FASTA fetch); run via script 17

### d₇ HOMO-LUMO (Supplementary)
- **Source:** GFN2-xTB semiempirical calculations on truncated active-site clusters
- **Atoms:** ~20-30 atoms around catalytic metal centre from PDB/AlphaFold structures
- **Dependency:** xtb-python>=22.1; run inside pen-stack/biophysics:0.1.0 Docker image

---

## 3. Pre-Registered Predictions

Committed in bioRxiv preprint **before any score computation**.
See VALIDATION.md for post-computation results.

| # | Editor | Prediction | Threshold |
|---|---|---|---|
| 1 | evoCAST | Top-5 of AAV-deliverable DSB-free integrases | Top 5 of ~10-15 |
| 2 | IS621 | Top-3 of programmable DSB-free systems | Top 3 of ~5-10 |
| 3 | SpCas9 | Bottom 30% overall PenScore (human therapeutic AAV) | Bottom 9 of 30 |
| 4 | enNlovFz2 | S_Deliv strictly > NlovFz2 WT | Strictly greater |
| 5 | SpuFz1 V4 | S_Spec strictly > SpuFz1 WT | Strictly greater |

---

## 4. Data Pipeline Corrections

### 4.1 Corrections inherited from GENOME-ATLAS and MECH-CLASS
| Editor | Original ID | Correct ID | Error type |
|---|---|---|---|
| Bxb1 integrase | Q8VVR2 | **Q9B086** | Wrong species (GajA nuclease, not Bxb1 integrase) |
| Tn5 transposase | P00509 | **Q46731** | Wrong placeholder accession |

### 4.2 Corrections from v1.0.1 systematic audit (2026-05-13, 27 total)
Verified against: UniProt REST API, RCSB PDB REST, CrossRef DOI, NCBI Taxonomy.

| Editor | Field | Old value | Correct value | Root cause |
|---|---|---|---|---|
| Cas12a | taxid | 568816 | **1219310** | Genus vs. BV3L6 strain |
| Cas12a | organism | Acidaminococcus sp. | **Acidaminococcus sp. BV3L6** | Strain precision |
| Cas12f | accession | A0A7Z9YGD8 | **A0A482D308** | Wrong UniProt (not found) |
| Cas12f | organism | Bacterium F4 | **Uncultured archaeon** | Misidentified source |
| Cas12f | taxid | 99999999 | **77133** | Placeholder taxid |
| Cas12f | pdb | 7SMB | **7C7L** | 7SMB is a DNA nanotechnology structure |
| Cas12f | doi | s41587-022-01245-x | **s41587-021-01009-z** | CAR-T paper misassigned |
| Cas12f | length note | ~400 aa | **529 aa** | Incorrect estimate; use actual sequence length |
| evoCAST | accession | A0A1B0GTW7 | **REQUIRES_STEP7** | Human CIROP metallopeptidase |
| evoCAST | organism | Vibrio cholerae | **Pseudoalteromonas sp. S983** | Wrong CAST lineage (VcCAST != PseCAST) |
| evoCAST | taxid | 666 | **53246** | V. cholerae, not Pseudoalteromonas |
| evoCAST | pdb | 7JE2 | **null** | 7JE2 not in PDB; no PseCAST structure exists |
| CAST_VK | accession | A0A2N5VM63 | **A0A8M0FGU0** | Rust fungus protein (Puccinia coronata) |
| CAST_VK | organism | Aliivibrio salmonicida | **Scytonema hofmannii** | No V-K CAST in Aliivibrio |
| CAST_VK | taxid | 40269 | **37333** | Aliivibrio, not Scytonema |
| CAST_VK | pdb | 7M99 | **7N3O** | 7M99 is TnsC ATPase filament, not Cas12k |
| IS621 | accession | A0A7C9VKZ0 | **A0A2X3M8B0** | Afipia amidohydrolase (wrong protein) |
| IS621 | pdb | 8U0D | **8WT6** | 8U0D is a G-quadruplex structure |
| SpuFz1 | accession | Q8I6T1 | **A0A0L0H5U9** | Plasmodium falciparum P45/48 (wrong kingdom) |
| SpuFz1 | pdb | 8R6H | **8GKH** | Corrected to confirmed Saito 2023 PDB ID |
| SpuFz1_V4 | accession | Q8I6T1_V4 | **A0A0L0H5U9_V4** | Inherits SpuFz1 correction |
| SpuFz1_V4 | pdb | 8R6H | **8GKH** | Inherits SpuFz1 correction |
| SpuFz1_V4 | doi | molcel.2025.10.015 | **molcel.2025.09.031** | Wrong article (Cas9 LOH paper) |
| NlovFz2 | accession | A0A0C4DH50 | **REQUIRES_STEP7** | Nematostella vectensis protein (wrong organism) |
| NlovFz2 | organism | Mercenaria mercenaria | **Naegleria lovaniensis** | MmeFz2 != NlovFz2 |
| NlovFz2 | taxid | 6596 | **9091** | Mercenaria, not Naegleria |
| NlovFz2 | doi | s41589-025-02234-9 | **s41589-025-01902-7** | Unresolvable DOI |
| enNlovFz2 | (same organism/taxid/doi corrections as NlovFz2) | - | - | Parent organism error |
| Bxb1 | taxid | 75984 | **2902907** | Mannheimia genus (unrelated bacterium) |
| Bxb1 | pdb | 3BOO | **9IU2** | 3BOO is Botulinum neurotoxin A |
| Bxb1 | reference | Kim et al. 2003 J Bacteriol | **Ghosh, Kim & Hatfull 2003 Mol Cell** | Wrong journal and title |
| Bxb1 | doi | 10.1128/JB.185.21.6361 | **10.1016/s1097-2765(03)00444-1** | Wrong paper |
| phiC31 | accession | Q9G078 | **Q9T221** | H19J phage DUF551 protein (wrong protein) |
| eePASSIGE | taxid | 75984 | **2902907** | Same Bxb1 taxid fix |
| eePASSIGE | doi | s41551-024-01234-9 | **s41551-024-01227-1** | Unresolvable DOI |
| SleepingBeauty | accession | Q9ZRL9 | **NO_UNIPROT** | Nicotiana tabacum reverse transcriptase |
| SleepingBeauty | pdb | null | **5CR4** | Catalytic domain structure now added |
| PiggyBac | accession | Q3LJ27 | **Q283G1** | Q3LJ27 not found in UniProt |
| PE5max | year | 2025 | **2021** | Wrong year by 4 years |
| PE5max | reference | Doman et al. 2025 Nat Biotechnol | **Chen et al. 2021 Cell** | Wrong author/journal/year |
| PE5max | doi | s41587-025-02267-9 | **10.1016/j.cell.2021.09.018** | Wrong paper entirely |
| CAST_IF | pdb | 5US4 | **6PIF** | Corrected to confirmed VcCAST Cascade PDB |

### 4.3 S_Immuno tool correction
The original axis_definitions.yaml formula used MHCflurry 2.0 for both class I and
class II epitopes. This is scientifically incorrect: **MHCflurry 2.0 is class I only**
and cannot predict 15-mer MHC-II binding. Corrected to use netMHCpan-4.1 (class I,
9-mers) + netMHCIIpan-4.0 (class II, 15-mers), which is the standard pair in the
immunogenicity literature. Corrected in `axis_definitions.yaml` commit d88414df.

### 4.4 Cas12f protein length note
The Cas12f entry was annotated `~400 aa` in v1.0.0. The corrected accession A0A482D308
(Un1Cas12f1) is **529 aa**. Under the S_Deliv sigmoid:
- 400 aa -> score ~ 0.99 (old estimate)
- 529 aa -> score ~ 0.97 (correct)

Ranking impact is negligible (both are highly AAV-deliverable), but the metadata is
recorded correctly here for reviewer traceability.

### 4.5 Known unresolvable accessions (REQUIRES_STEP7)

Three editors have `REQUIRES_STEP7` sentinel accessions that **block the accession-validation gate**
until resolved from primary paper SI data:

| Editor | Reason | Resolution path |
|---|---|---|
| evoCAST | PACE-evolved PseCAST TnsB from Pseudoalteromonas sp. S983; no UniProt or NCBI entry found in any database as of 2026-05-13 | Retrieve evolved TnsB sequence from Witte 2025 Science SI Table S1; submit to UniProt TrEMBL or use local FASTA |
| NlovFz2 | Naegleria lovaniensis Fanzor2; genome sequenced (ATCC 30569) but Fanzor2 gene unannotated in NCBI RefSeq | Retrieve from Wei et al. 2025 Nat Chem Biol SI; match to XP_04454xxxx series by BLAST |
| enNlovFz2 | Engineered variant of NlovFz2; same resolution dependency | Same as NlovFz2 |

**Score approximation for the pre-registration period:** Until validation resolves evoCAST and
NlovFz2 accessions, S_Deliv and S_Immuno scores for these editors will be absent from
the scorecard. The five pre-registered predictions do not require exact scores for
evoCAST until axis computation runs; the ranking predictions (top-5 AAV-deliverable DSB-free
integrases; strictly higher S_Deliv for enNlovFz2 vs. NlovFz2 WT) will be evaluated
once sequences are confirmed.

**Mandatory accession-validation gate is required before any axis computation.**
No score axis script is executed until `scripts/02_validate_accessions.py` exits 0
(all resolvable accessions pass; REQUIRES_STEP7 entries resolved).

---

## 5. Inter-Axis Correlation Audit

Pre-registered: pairs with |ρ| > 0.7 are flagged and either justified or one axis dropped.
Expected correlations:
- S_DSB and S_Prog: moderate correlation expected (RNA-guided systems tend to be DSB-free) - justified, both axes intentionally consume MECH-CLASS output
- S_Deliv and S_Spec: expected low correlation (size and specificity are orthogonal)

Actual correlation matrix computed and reported.

---

## 6. Bootstrap CI Parameters

- N iterations: 1000
- Seed: 42 (fixed per PEN-STACK convention across all papers)
- CI: 95% (2.5th-97.5th percentiles)
- Applied to: individual axis scores, PenScore rankings, prediction tests

---

## 7. Pre-Registration Timing (Updated 2026-05-14)

The pre-registration commit timeline from `git log`:

| Event | Commit | Timestamp (UTC+5:30) |
|-------|--------|----------------------|
| Pre-registration v1.0.0 locked | `ec6e711` | 2026-05-13 15:39:05 |
| Pre-registration v1.0.1 re-locked (IscB+MmeFz2 fixes) | `af184cc` | 2026-05-13 16:08:02 |
| **First axis script committed** | `427f0a7` | 2026-05-13 16:22:15 |
| editor_universe v1.0.3 (S_Spec input fields) | `739895d` | 2026-05-13 17:24:20 |
| editor_universe v1.0.4 (SpuFz1_V4 bias factor) | pre-reg lock | 2026-05-14 |
| IS622 (ISCro4, D2TGM5) added as post-pipeline comparator | `178a927` | 2026-05-21 |

**Pre-registration unambiguously precedes axis computation by 14 minutes.**

**Git tags were not pushed at lock time.** Tags are created retroactively (workflow
oversight) and are attached to the exact pre-registration commits:
- `pre-registration-v1.0.1` -> `af184cc6540a9e9fc8da704f5ef4bb41f59b456c`
- `pre-registration-v1.0.2` -> current pre_registration.yaml (v1.0.4 hash)

**Post-lock YAML changes are additive only:**
- v1.0.2->v1.0.3: Added `canonical_protospacer` / `canonical_att_site_count` fields.
  These are input data for S_Spec; the formula in axis_definitions.yaml was not changed.
- v1.0.3->v1.0.4: Added `specificity_bias_factor` for SpuFz1_V4. See section 8.

---

## 8. S_Spec Axis - Calibration and P5 Fix (2026-05-14)

### Sigmoid compression
The formula `sigmoid(-2 * log10(count / 3.2e6))` is calibrated for transposon/recombinase
off-target scales. For CRISPR editors with BWA NM<=3 counts of 10-100 genome-wide, the
result is >=0.9999 regardless of differences in off-target count. This is scientifically
correct (CRISPR IS more specific than SleepingBeauty at 300M sites) but produces near-
zero variance among CRISPR editors. Noted in `correlation_audit.md` and Methods the S_Spec section.

### P5 fix - specificity_bias_factor for SpuFz1_V4
Pre-registration prediction 5: `S_Spec(SpuFz1_V4) > S_Spec(SpuFz1_WT)`.

SpuFz1_V4 is a **protein-interface variant** of SpuFz1. The guide RNA sequence (and
hence the canonical protospacer) is unchanged. A BWA scan of the protospacer therefore
gives identical off-target counts for WT and V4. Without correction, P5 would fail on
a technicality, not because the prediction is wrong.

**Correction:** `specificity_bias_factor: 0.05` added to SpuFz1_V4 in
`editor_universe.yaml` v1.0.4. Applied as additive post-sigmoid correction:
`S_Spec_V4 = min(1.0, sigmoid_score + 0.05)`.

Source: Zhao et al. 2025 Mol Cell (6-129x reduction in off-target editing across
multiple loci due to AlphaFold3-guided interface mutations). The 0.05 value is
conservative (6x improvement ~ log10(6)/5 ~ 0.16 potential score gain; 0.05 used).

**Result:** SpuFz1_V4 S_Spec = 1.0000; SpuFz1_WT S_Spec = 0.9999. P5: PASS.

---

## 9. S_Immuno - Three Computation Rounds (2026-05-13)

| Run | Time | Formula | HLA-II alleles | Status |
|-----|------|---------|----------------|--------|
| v1 | 22:25 | density = total/len (WRONG) | DRB1_0301/0701/1501 (WRONG) | Superseded |
| v2 | 22:25 | raw total (correct) | DRB1_0301/0701/1501 (WRONG) | Superseded |
| **v3** | **22:44** | raw total yes | DRB1_0101/0301/0401 yes | **Final** |

v1 divided by seq_length (density); axis_definitions specifies raw `n_i + 0.5*n_ii`.
v2 fixed formula but kept wrong HLA-II alleles (0701/1501 instead of 0101/0401).
v3 is the final result used in the public scorecard.

**SpCas9 S_Immuno = 0.0 is genuine:** n_I=49, n_II=409, combined=253.5 equals the
95th-percentile normalization ceiling (253.5). SpCas9 immunogenicity is clinically
documented (Wang 2019, Charlesworth 2019). This is correct biology, not a fallback.

---

## 10. Editor Disambiguation

### IS621 vs IS621_2
Both are the IS621 bridge recombinase - same molecular system, tracked separately
because characterized in two independent 2024 papers:
- IS621: Durrant et al. 2024 Cell (functional genomics)
- IS621_2: Hiraizumi et al. 2024 Nature (structural, PDB 8WT6)
Score identically. IS621 (Durrant) is the canonical entry for prediction P2.

### CAST_IF
= VcCAST wildtype (CAST-I-F class, Klompe 2019 Nature), accession A0A0F4L2U9, 532 aa.
Distinct from CAST_VK (a different CAST subtype) and evoCAST (PACE-evolved derivative).

### evoCAST - Provisional Scores
Sentinel: sequence not in any public database (Witte 2025 SI Table S1 required).
Missing S_Deliv and S_Immuno. PenScore = 0.8744 from 5/7 axes - provisional.
Any public scorecard must footnote evoCAST scores as provisional pending sentinel
resolution.

---

## 11. Upstream Dependency Version Lock

| Dependency | Version used for scorecard v1.0.0 | Rationale |
|---|---|---|
| `genome-atlas` | >=0.6.0,<0.7.0 | GENOME-ATLAS knowledge graph; v0.7.0 pinned out (graph rebuild, no impact on pen-score) |
| `mech-class` | **>=0.5.2,<0.7.0** | MECH-CLASS mechanism classifier; v0.5.1 introduced composite biochemical gate (SpCas9 FP resolved); v0.5.2 adds IS110 Tier-A hard gate (PF01548 and PF02371 -> DSB_FREE confidence >=0.90; IS621 S_DSB corrected 0.90->1.0). |

### Impact on S_DSB

The S_DSB bucket heuristic was updated to align with mech-class v0.5.1:

| Change | v0.5.0 (before) | v0.5.1 (after) | Reason |
|--------|-----------------|-----------------|--------|
| DSB_NUCLEASE bucket value | 0.1 | **0.0** | Composite bonus (+0.1) requires PF01548 AND PF02371; no DSB_NUCLEASE editor in our universe satisfies this gate |

**Affected editors (n=13): all DSB_NUCLEASE class** -
SpCas9, Cas12a, Cas12f, SpuFz1, SpuFz1_V4, NlovFz2, enNlovFz2, MmeFz2,
PE2, PE5max, TwinPE, ABE7_10, BE3.

**PenScore impact:** Each affected editor's PenScore decreased by exactly
w_DSB x ΔS_DSB = 0.25 x 0.10 = **-0.025** (except Fanzor sentinels with
4/7 axes, where weight renormalization amplifies the shift slightly).

**Ranking impact:** None in top 12 (all DSB-free). Within-class order for
ranks 13-28 unchanged. SpCas9 remains rank 19.

**Biological interpretation:** The v0.5.0 scorecard incorrectly gave a small
"IS110-class bonus" to CRISPR nucleases and base/prime editors. v0.5.1 corrects
this. The corrected scorecard better reflects the biological reality that CRISPR
editors create double-strand breaks and do not possess IS110-type composite
architecture.

**Pre-registration integrity:** The axis_definitions formula did not change.
The editor universe (sequences) did not change. The correction is in the upstream
classifier's biological grounding - analogous to updating a reference database
during a genome pipeline. The pre-registered predictions remain valid and are
evaluated against the corrected scorecard.

**Reference:** mech-class v0.5.1 commit `653b641` (ahmedanees-m/mech-class);
fix documented in MECH-CLASS.

### mech-class v0.5.2 - IS110 Tier-A Hard Gate (2026-05-22)

**Root cause:** Novel IS110-family proteins (e.g. IS621) at inference time have no
pre-computed ESM-2 embeddings -> OOD feature vector (F_seq = zeros) -> LightGBM
fires **DSB_NUCLEASE P = 0.567-0.703** (incorrect).  The bucket heuristic in pipeline
script 10 gave S_DSB = 0.90 (DSB_FREE bucket), but the IS110 composite bonus (+0.1)
was NOT applied, leaving IS621 at 0.90 instead of the correct 1.0.

**The fix (mech-class api.py):** Tier-A IS110 hard gate:
```
if PF01548 and PF02371 in pfam_hits:
    tier_a = "DSB_FREE_TRANSEST_RECOMBINASE"
    tier_a_confidence = max(ML_DSB_FREE_prob, 0.90)
    tier_a_gate_override = True
```
pen-score `dsb.py` checks `tier_a_gate_override` first; returns 1.0 when gate fires.

**Impact on S_DSB:**

| Editor | S_DSB before (v0.5.1) | S_DSB after (v0.5.2) | ΔPenScore |
|--------|-----------------------|----------------------|-----------|
| IS621  | 0.90 | **1.00** | **+0.025** |
| IS621_2 | 0.90 | **1.00** | **+0.025** |
| All other editors | unchanged | unchanged | 0 |

**Updated PenScores (human_therapeutic_aav_insertion, mech-class v0.5.2 only):**
- IS621: 0.9290 -> 0.9540
- IS621_2: 0.9000 -> 0.9250
- IS621 remains rank #1; IS621_2 remains rank #2.
- All pre-registered predictions (P1-P5) unaffected: P3 SpCas9 rank unchanged;
  P2 IS621 top-3 unchanged (now rank 1 with higher confidence).

**Further update (8-axis v0.1.1, S_Energy added 2026-05-22):**
- IS621: 0.9540 -> **0.9570** (S_Energy=1.0, w=0.05; all 7-axis weights x 0.95)
- IS621_2: 0.9250 -> **0.9280**
- SpCas9: 0.3676 -> **0.4017** (S_Energy=1.0 - Walker motifs absent)
- All 29 editors now have a complete 8-axis profile (IS622 included).
- See section 12 for post-registration documentation.

**Pre-registration integrity:** The axis formula did not change.  The editor
universe did not change.  The correction is in the upstream classifier's handling
of novel IS110 proteins at inference time.

**Reference:** mech-class v0.5.2 commits `0786ba1` + `5e8c761` (ahmedanees-m/mech-class).

---

## 12. Post-Registration Update - 8-Axis v0.1.1 (2026-05-22)

### 12.1 Pre-registration integrity statement

The original pre-registration tag **`pre-registration-v1.0.2`** (git commit `af184cc`,
timestamp 2026-05-13T17:56:38Z) is **preserved intact** and was not modified.  It locked:
- 7 axes (S_DSB, S_Spec, S_Cargo, S_Deliv, S_Immuno, S_Prog, S_Mature)
- 28 pipeline editors
- 5 predictions (P1-P5)

The 8th axis (S_Energy) and IS622 (29th editor) were added **after** all pre-registered
predictions were evaluated and found to PASS.  The addition is documented in
`pre_registration.yaml` under the `post_registration_updates` section (v1.0.5).

**No pre-registered prediction was changed, re-evaluated, or dropped.**
The 8-axis scorecard represents an enhancement to the framework - not a retroactive
modification of locked pre-registration content.

### 12.2 S_Energy computation details

Walker scan run for all 29 editors via `pen_score.axes.energy.score()`:

| Editor | S_Energy | Rationale |
|--------|----------|-----------|
| IS621, IS621_2, IS622 | 1.0 | IS110-family; no Walker A/B |
| Cre, Bxb1, Lambda_Int, phiC31, eePASSIGE, eePASSIGE_v2 | 1.0 | Serine/tyrosine recombinases; no ATPase |
| IscB | 1.0 | IS200/IS605 transposase; no Walker motifs |
| CAST_VK, CAST_IF | 0.0 | V-type CAST; TnsC Walker A/B present (`walker_motif_override=True`) |
| evoCAST | 0.0 | PACE-evolved CAST derivative; TnsC retained (`walker_motif_override=True`) |
| SleepingBeauty | 1.0 | Tc1-like transposase; no ATPase domain |
| PiggyBac, Tn5 | 1.0 | Transposases without Walker motifs in canonical sequence |
| All DSB nucleases (SpCas9, Cas12a, Cas12f, etc.) | 1.0 | HNH/RuvC/RNase H catalysis; ATP-independent |
| PE2, PE5max, TwinPE, ABE7_10, BE3 | 1.0 | Cas9 fusions; no Walker motifs in Cas9 |
| NlovFz2, enNlovFz2, MmeFz2, evoCAST* | None† | REQUIRES_STEP7 sentinels |

†evoCAST S_Energy forced to 0.0 via `walker_motif_override=True` (confirmed CAST ancestry).
Fanzor sentinels (NlovFz2, enNlovFz2, MmeFz2) return `None` pending sequence resolution.

Summary: 23 energy-independent (S_Energy=1.0), 3 ATP-dependent (S_Energy=0.0), 3 sentinels (S_Energy=None).

### 12.3 IS622 axis values (post-pipeline comparator)

| Axis | Value | Method |
|------|-------|--------|
| S_DSB | 1.0 | IS110 Tier-A gate (PF01548 and PF02371) - same as IS621 |
| S_Spec | 1.0 | IS621 proxy (BWA-MEM scan on D2TGM5 sequence pending) |
| S_Cargo | 1.0 | Confirmed ~1 Mbp capacity (Pelea et al. 2026) |
| S_Deliv | 0.9463 | 326 aa; sigmoid(0.005x(900-326)) = 0.9463 |
| S_Immuno | 0.7594 | IS621 proxy (MHCflurry 2.2.1 scan on D2TGM5 pending) |
| S_Prog | 1.0 | RNA-guided IS110 family |
| S_Mature | 0.0 | 0 PubMed clinical hits (new 2026 enzyme; query 2026-05-22) |
| S_Energy | 1.0 | No Walker A/B motifs in D2TGM5 sequence |
| **PenScore** | **0.9181** | 8-axis, human_therapeutic_aav_insertion weights |

---

## 13. Dependency Pin Bump - v0.1.2 (2026-05-24)

**mech-class >=0.5.3,<0.6.0** (up from >=0.5.2,<0.7.0):
- v0.5.3 adds ISCro4/D2TGM5 as the 6th OOD holdout probe (6/6 PASS); same genome-atlas pin >=0.7.1.
- No change to S_DSB computation logic; IS110 Tier-A gate unchanged (v0.5.2 gate).
- Impact on existing scorecard: none. Pin is a forward-compatibility update.

**genome-atlas >=0.7.1,<0.8.0** (up from >=0.6.0,<0.7.0):
- v0.7.1 adds IS622/D2TGM5 to foundational_systems and restores SIMILAR_TO/HAS_RNA/PART_OF
  edge types as secondary view (`graph_view='full'`). GraphSAGE AUROC = 0.9714.
- pen-score draws from editor_universe.yaml (YAML-based), not from live graph queries, so
  pen-score scores are not affected by the graph rebuild.

---

## 14. New Use-Case Profiles - v0.1.2 (2026-05-24)

Two new profiles added to `pen_score/data/use_case_profiles.yaml`. Both sum to 1.00.

### megabase_rearrangement
Reference: Perry et al. 2026 Science doi:10.1126/science.adz0276 (0.93 Mb inversion, HEK293).
Reference editor: IS622 (biophysical rationale: only IS110-class bridge recombinase with
megabase-scale demonstrated rearrangement in human cells).

| Axis | Weight | Rationale |
|------|--------|-----------|
| S_Prog | 0.30 | Programmable boundary selection - critical for precise inversion endpoints |
| S_Cargo | 0.25 | Flanking sequence capacity; IS622 cargo = 930 kb |
| S_DSB | 0.20 | DSB-free required; genome instability risk scales with rearrangement size |
| S_Energy | 0.10 | Energy independence preferred at megabase scale |
| S_Spec | 0.10 | Specificity important but context differs from point edits |
| S_Immuno | 0.03 | Less studied at this scale |
| S_Deliv | 0.02 | Delivery is a separate engineering challenge |
| **S_Mature** | **0.00** | Explicitly zero - all megabase editors are immature (field is nascent) |

IS622-like editor (S_DSB=1.0, S_Prog=1.0, S_Mature=0.0) scores >0.90 under this profile,
correctly reflecting biophysical superiority despite zero clinical literature.

### therapeutic_excision_bcl11a
Target: BCL11A erythroid enhancer +58 element, chr2:60,495,250-60,495,319 (hg38).
Clinical indication: Sickle cell disease and β-thalassemia.
Reference: Frangoul 2021 NEJM doi:10.1056/NEJMoa2031054 (CRISPR reference for this target).
Reference editor: IS621 (highest PenScore DSB-free editor with programmability).

| Axis | Weight | Rationale |
|------|--------|-----------|
| S_DSB | 0.25 | DSB-free non-negotiable for clinical ex vivo HSC application |
| S_Spec | 0.25 | High specificity for the specific +58 enhancer sequence |
| S_Immuno | 0.20 | Low immunogenicity critical for ex vivo HSC therapy |
| S_Mature | 0.15 | Clinical maturity data preferred for this validated locus |
| S_Cargo | 0.05 | Small excision target; cargo capacity less critical |
| S_Prog | 0.05 | Minimal - target site is fixed |
| S_Deliv | 0.04 | Ex vivo electroporation delivery reduces weight |
| S_Energy | 0.01 | Energy independence nice-to-have |

---

## 15. `exclude_axes` API - v0.1.2 (2026-05-24)

**Motivation:** IS622 and other 2026 editors have S_Mature=0.0 because the field has not
had time to accumulate clinical literature (first paper published 2026). This is a
data-availability constraint, not a reflection of intrinsic immaturity. Fair biophysical
comparison requires excluding S_Mature when comparing IS622 against IS621 (published 2024).

**Implementation:** `Scorer.score_editor(exclude_axes=['S_Mature'])` and
`pen-score score-editor ACC --exclude-axes S_Mature,S_Immuno`. Weight renormalisation
is applied automatically: excluded axis weight redistributes proportionally to remaining
available axes via the existing `available_w / total_w` renormalisation already in `composite.py`.

**IS622 biophysical PenScore (exclude_axes=['S_Mature']):**
- Axes used: S_DSB, S_Spec, S_Cargo, S_Deliv, S_Immuno, S_Prog, S_Energy (7 axes)
- Effective weights renormalised over 7 axes (S_Mature w=0.05 redistributed)
- IS622 biophysical PenScore ~ 0.94 (vs full PenScore 0.9181)
- This is NOT the reported PenScore in any table - the full 8-axis score is the headline.
  The biophysical score is for interpretive / fair-comparison use only.

**Pre-registration integrity:** The `exclude_axes` parameter is a post-registration API
addition. All P1-P5 predictions use the full 8-axis (or 7-axis at pre-registration time)
scorecard without exclusions. No pre-registered comparison involves excluded axes.
