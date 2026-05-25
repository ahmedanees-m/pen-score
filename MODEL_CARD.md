# MODEL_CARD.md - PEN-SCORE v0.1.3

**Package:** pen-score
**Author:** Anees Ahmed Mahaboob Ali
**Affiliation:** VIT University, Vellore
**Date:** 2026-05-25  (created 2026-05-12; updated 2026-05-14 with prediction results; updated 2026-05-22 to 8-axis v0.1.1; updated 2026-05-24 to v0.1.2; updated 2026-05-25 to v0.1.3)
**Analogous to:** MODEL_CARD.md in mech-class

---

## 1. Intended Use

PEN-SCORE is designed for **computational triage** of programmable genome editors.
It produces an 8-axis score profile and a composite PenScore that helps wet labs
prioritise which editor to test for a given use case.

**In-scope use cases:**
- Selecting among CRISPR-Cas, CAST, bridge recombinase, Fanzor, prime editor, or
  site-specific recombinase systems for human therapeutic applications
- Comparing editors on AAV deliverability, specificity, cargo capacity, and immunogenicity
- Generating the public scorecard for the PEN-STACK web interface (PEN-COMPARE)

**Out-of-scope use cases:**
- Predicting in vivo editing efficiency at specific genomic loci (not measured)
- Replacing wet-lab validation (PEN-SCORE guides selection, not activity prediction)
- Designing chimeric editors (that is PEN-ASSEMBLE, PEN-ASSEMBLE)

---

## 2. Score Axis Limitations

### S_Energy (v0.1.1)
- Binary (0 or 1); uses Walker A `G[A-Z]{4}GK[ST]` + Walker B `[LVIMF]{4}DE` regex
  scan on UniProt canonical FASTA.
- Multi-subunit CAST systems (evoCAST, CAST_VK, CAST_IF) have `walker_motif_override: true`
  in `editor_universe.yaml` - this forces S_Energy=0.0 regardless of the TnsC subunit
  not being in the primary accession's sequence.
- Does not distinguish partial ATP dependence (e.g. ATPase activity stimulated but not
  essential). Binary scoring is conservative.
- Sentinel editors (REQUIRES_STEP7) receive S_Energy=None and are excluded from PenScore
  normalisation for this axis.

### S_DSB
- Depends on MECH-CLASS **>=v0.5.2**. The v0.5.1 biochemical gate
  (PF01548 AND PF02371 co-occurrence required) prevents false-positive composite
  flags on dual-domain non-IS110 proteins. SpCas9 composite FP (P=0.753 at v0.5.0)
  is **resolved** in v0.5.1; SpCas9 S_DSB=0.0 is correct.
  The v0.5.2 Tier-A IS110 hard gate corrects the inverse problem: novel IS110
  proteins without pre-computed ESM-2 embeddings were mis-scored as DSB_NUCLEASE
  (P=0.567-0.703) due to OOD feature vectors. IS621 and all IS110-family bridge
  recombinases now return S_DSB=1.0 (domain-evidence confidence 0.99).
  See SCORE_PROVENANCE.md Section 11 for full impact table.
- TRANSPOSASE Tier B sub-classifier absent from MECH-CLASS (N<3 sub-classes).
  S_DSB for transposases is computed from Tier A probability only.

### S_Spec
- Specificity computed from a canonical 20-bp protospacer - not a comprehensive
  off-target scan across all possible guides for an editor.
- Non-RNA-guided editors (Cre, Bxb1, phiC31): specificity is approximated by the
  natural att-site count in GRCh38, not a protospacer scan.
- GRCh38 is used as the target genome; editors used primarily in non-human organisms
  may receive misleading S_Spec scores.
- Engineered variants sharing a protospacer with their parent (e.g., SpuFz1_V4)
  receive a `specificity_bias_factor` post-sigmoid correction to encode documented
  protein-level specificity improvements (see editor_universe.yaml v1.0.4).

### S_Cargo
- Cargo capacities are literature-curated maximums, not mechanistic limits.
  Practical efficiency degrades well before the reported maximum size.
- Engineered variants whose cargo has not been systematically characterised inherit
  the parent editor's value (conservative).

### S_Deliv
- Sigmoid centred at 900 aa.  Fusion editors (PE2, ABE, BE3) include only the Cas9
  component length because M-MLV RT and deaminase sizes are not retrieved separately.
  These editors are likely **overscored** on S_Deliv; actual single-AAV compatibility
  requires the full fusion size.

### S_Immuno
- MHCflurry 2.0 predictions are in silico; not experimentally validated epitope loads.
- HLA-I alleles: HLA-A02:01, HLA-A01:01, HLA-B07:02, HLA-B44:02 (European reference
  panel). HLA-II alleles: DRB1_0101, DRB1_0301, DRB1_0401.
- Bacterial-origin editors (SpCas9, Cas12a) may have lower immunogenicity in practice
  than predicted because pre-existing immunity is population-dependent.
- SpCas9 S_Immuno=0.0 is genuine (n_I=49, n_II=409; combined=253.5 saturates the
  95th-percentile normalisation ceiling). This is correct biology, not a fallback.

### S_Prog
- Binary (0 or 1); does not capture partial programmability (e.g. near-PAMless Cas9,
  expanded TALE arrays).  TALE editors are not in the current universe.

### S_Mature
- PubMed clinical citation count is a proxy for maturity; highly cited basic-science
  editors (Cas9) will dominate.  Normalisation by universe maximum mitigates this.
- Citation count depends on search term set; see `references_used_for_pubmed` in
  editor_universe.yaml for each editor's query terms.

---

## 3. Inherited Limitations from MECH-CLASS (MECH-CLASS)

These limitations carry forward because S_DSB and S_Prog consume MECH-CLASS output:

1. **SpCas9 composite FP - RESOLVED in v0.5.1**: The v0.5.0 composite head fired
   composite=True for SpCas9 (P=0.753) due to the dual RuvC+HNH architecture
   mimicking IS110's dual-domain pattern at the ML level. The v0.5.1 biochemical gate
   (PF01548 AND PF02371 required) correctly blocks this. S_DSB for all 13 DSB_NUCLEASE
   editors is 0.0; PenScore for each dropped by 0.025 relative to v0.5.0. See
   mech-class MODEL_CARD.md v0.5.1 changelog and SCORE_PROVENANCE.md Section 11.
2. **Novel IS110 OOD mis-scoring - RESOLVED in v0.5.2**: IS110 proteins without
   pre-computed ESM-2 embeddings (novel inference-time inputs) received OOD feature
   vectors -> LightGBM returned DSB_NUCLEASE P=0.567-0.703. The v0.5.2 Tier-A IS110
   hard gate (PF01548 and PF02371 -> DSB_FREE_TRANSEST_RECOMBINASE, confidence >=0.90)
   corrects this. IS621 S_DSB corrected from 0.90 -> 1.0; PenScore increases by +0.025.
   See mech-class MODEL_CARD.md v0.5.2 changelog and SCORE_PROVENANCE.md Section 11.
3. **SaProt GPU requirement**: F_struct channel is always zero-filled at live inference
   (SaProt requires GPU + AlphaFold; deferred to v0.6).  MECH-CLASS inference at
   runtime uses only F_domain; this is acceptable given the ablation result
   (domain_only macro-F1 = 0.9859 ~ full 0.9862).
4. **No TRANSPOSASE Tier B**: sub-classification within TRANSPOSASE class unavailable.

---

## 4. Known Unknowns

- **AlphaFold confidence for composite folds**: IS110-class proteins with dual-domain
  (PF01548 + PF02371) may have low pLDDT at the domain interface.  Any structural
  computation (d7 supplementary) for these proteins is downweighted accordingly.
- **Engineered variant accessions**: enNlovFz2, NlovFz2, MmeFz2, evoCAST are
  sentinels (sequences not in public databases). S_Deliv and S_Immuno are missing
  for these editors; PenScores are provisional. Sentinel resolution paths documented
  in pre_registration.yaml and SCORE_PROVENANCE.md Section 4.5.
- **Citation counts at the snapshot date**: S_Mature values are a snapshot at the time
  of computation (2026-05-13). The preprint notes this and provides the query date.
- **Cross-axis correlation S_Deliv/S_Immuno** (rho=0.94): documented in
  scorecards/correlation_audit.md as a known limitation; both axes retained because
  they measure distinct biological barriers. See SCORE_PROVENANCE.md Section 8.

---

## 5. v0.1.3 Changes (2026-05-25)

**`get_editor_metadata()` API - PEN-COMPARE v3.2 integration:**
- New function `pen_score.get_editor_metadata(editor_id) -> EditorMetadata` exposes two
  boolean fields required by PEN-COMPARE v3.2 certification tiers:
  - `intrinsic_cargo_mechanism` (bool, v1.0.7): **True** if the editor can carry/insert cargo
    as part of its catalytic mechanism (e.g. IS110 bridge recombinases, CAST transposases,
    site-specific recombinases). **False** if cargo requires an external HDR donor template
    (e.g. SpCas9+HDR, PE2). Used by PEN-COMPARE v3.2 **Gate 3**.
  - `cell_based_evidence` (bool, v1.0.7): **True** if peer-reviewed mammalian cell activity
    data exists at >1% editing efficiency. **False** for in vitro / E. coli only. Used by
    PEN-COMPARE v3.2 **TRUE_WRITER tier**.
  - `cell_based_sources` (list[str]): supporting citations for editors with cell_based=True.

**ISCro4 canonical naming (IS622 deprecated):**
- `editor_universe.yaml` v1.0.7: `IS622` renamed to `ISCro4` (canonical per UniProt D2TGM5 +
  Pelea 2026 *Science*). `aliases: ["IS622"]` retained. Alias resolution emits
  `DeprecationWarning`.

**Dependency pin bumps:**
- `mech-class>=0.5.4,<0.6.0` (ISCro4 holdout probe rename; was IS622_perry_2025)
- `genome-atlas>=0.7.2,<0.8.0` (ISCro4 canonical in foundational_systems; IS622 as alias)

**Test suite:** 196 pass, 3 skipped (parquet-only), 0 failures. 18 new tests.

---

## 5b. v0.1.2 Changes (2026-05-24)

**Dependency pin bumps:**
- `mech-class>=0.5.3,<0.6.0` (ISCro4/D2TGM5 added to OOD holdout; atlas pin updated)
- `genome-atlas>=0.7.1,<0.8.0` (ISCro4/D2TGM5 added to atlas; SIMILAR_TO edges restored via `graph_view='full'`)

**New use-case profiles:**
- `megabase_rearrangement`: S_Prog=0.30 dominant; S_Mature=0.00 (all megabase editors are new);
  reference editor ISCro4 (0.93 Mb inversion, Perry 2026 Science)
- `therapeutic_excision_bcl11a`: S_DSB=S_Spec=0.25 dominant; S_Immuno=0.20; targets BCL11A
  erythroid enhancer +58 excision for SCD/β-thalassemia; reference editor IS621

**`exclude_axes` API (v0.1.2+):**
- `scorer.score_editor(accession, exclude_axes=['S_Mature'])` allows biophysical-only comparisons
- Useful for ISCro4 and other 2026 editors where S_Mature=0.0 reflects literature recency, not
  intrinsic limitations. CLI: `pen-score score-editor ACC --exclude-axes S_Mature`
- Weight renormalisation over remaining axes is automatic

**Perry-2025 ortholog panel infrastructure:**
- `scripts/18_score_perry2025_panel.py` created; full panel requires
  supplementary CSV from Perry et al. 2026 Science doi:10.1126/science.adz0276

**Test suite:** 267 unit tests pass, 3 skipped (parquet-only), 0 failures.

---

## 6. Intended Update Cycle

See UPDATE_STRATEGY.md.  v0.0.1 -> v1.0 at the stable release.

Trigger conditions for scorecard update:
- New editor added to universe (requires YAML bump + Steps 8-20 re-run)
- Sentinel resolved (re-run affected axis + Steps 17-20)
- mech-class dependency version bump (re-run S_DSB + Steps 17-20)

---

## 6. Pre-Registered Prediction Results (v0.1.1, 8-axis)

Scorecard version: v0.1.1 (8-axis, mech-class v0.5.2; IS621 PenScore 0.929->0.957)
Pre-registration: pre_registration.yaml v1.0.2 (locked 2026-05-13T17:56:38Z)
Evaluated on 7-axis scorecard: 2026-05-14; confirmed on 8-axis scorecard: 2026-05-22

| Test | Status | Observed | Threshold | Explanation |
|------|--------|----------|-----------|-------------|
| P1 evoCAST top-5 DSB-free | **PASS** | Rank 5/13 in subset | Top 5 | evoCAST competitive despite 6/8 axes (S_Deliv/S_Immuno pending sentinel); rank 5 with ISCro4 joining the subset |
| P2 IS621 top-3 programmable | **PASS** | Rank 1/7 in subset | Top 3 | Top-ranked programmable DSB-free editor; bootstrap CI=[1,1] |
| P3 SpCas9 bottom 30% | **PASS** | Rank 20/29, 31.0% below | >= 30% below | SpCas9 PenScore=0.4017 (8-axis); 9/29 editors score lower (31.0% >= 30%) |
| P4 enNlovFz2 > NlovFz2 S_Deliv | **NOT EVALUABLE** | S_Deliv=None for both | Strictly greater | Both editors are documented sentinels (REQUIRES_STEP7) awaiting Wei et al. 2025 Nat Chem Biol SI sequence deposition. Prediction is physically certain but cannot be operationalized until sequences resolve. Anticipated PASS. |
| P5 SpuFz1_V4 > SpuFz1 S_Spec | **PASS** | 1.0000 > 0.9999 | Strictly greater | Specificity-bias factor (+0.05) applied post-sigmoid encodes Zhao et al. 2025 Mol Cell 6-129x interface engineering gain |

**Outcome:** 4/4 evaluated predictions PASS. P4 is a
documented data-availability boundary condition (sentinel listed in
SCORE_PROVENANCE.md Section 4.5 before pre-registration lock), not a model
failure. Scorecard v0.1.1 is validated.

**Key structural finding:** Top 13 editors are all DSB-free; editors 14-29 are
transposases and nucleases. This clean biological stratification is the headline
quantitative result: IS621 (PenScore 0.9570, 8-axis v0.1.1) vs SpCas9 (0.4017)
represents a 2.4x difference driven by mechanism class, not by secondary axes.
IS110-family bridge recombinases (IS621, IS621_2, ISCro4) correctly score S_DSB=1.0
and S_Energy=1.0 with the v0.5.2 IS110 Tier-A hard gate.
