"""
scripts/41_generate_outreach_materials.py
=========================================
Generate 5 wet-lab outreach Markdown summaries for potential collaborators.

Each summary includes:
  (a) Lab research focus and relevance to PEN-SCORE
  (b) Editors from the public scorecard relevant to the lab's work
  (c) Top-5 recommendations with axis-by-axis reasoning
  (d) PEN-ASSEMBLE collaboration invitation

Output: ~/pen-stack/data/pen-score/outreach/<slug>.md  (5 files)
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Load scorecard
# ---------------------------------------------------------------------------

SCORECARD_PATH = Path.home() / "pen-stack/data/pen-score/scorecards/public_scorecard.parquet"
OUTREACH_DIR   = Path.home() / "pen-stack/data/pen-score/outreach"
OUTREACH_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_parquet(SCORECARD_PATH)

# ---------------------------------------------------------------------------
# Axis helpers
# ---------------------------------------------------------------------------

AXES   = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature"]
LABELS = {
    "S_DSB":    "DSB safety",
    "S_Spec":   "Specificity",
    "S_Cargo":  "Cargo capacity",
    "S_Deliv":  "Deliverability",
    "S_Immuno": "Immunogenicity",
    "S_Prog":   "Programmability",
    "S_Mature": "Tech maturity",
}

USE_CASE_WEIGHTS: dict[str, dict[str, float]] = {
    "human_therapeutic_aav_insertion": {
        "S_DSB":0.25,"S_Spec":0.15,"S_Cargo":0.20,"S_Deliv":0.20,
        "S_Immuno":0.10,"S_Prog":0.05,"S_Mature":0.05
    },
    "human_therapeutic_electroporation": {
        "S_DSB":0.25,"S_Spec":0.25,"S_Cargo":0.20,"S_Deliv":0.05,
        "S_Immuno":0.10,"S_Prog":0.05,"S_Mature":0.10
    },
    "large_cargo_integration": {
        "S_DSB":0.20,"S_Spec":0.10,"S_Cargo":0.40,"S_Deliv":0.10,
        "S_Immuno":0.05,"S_Prog":0.10,"S_Mature":0.05
    },
    "base_editing_small_correction": {
        "S_DSB":0.15,"S_Spec":0.40,"S_Cargo":0.05,"S_Deliv":0.15,
        "S_Immuno":0.15,"S_Prog":0.05,"S_Mature":0.05
    },
    "research_discovery": {
        "S_DSB":0.15,"S_Spec":0.15,"S_Cargo":0.15,"S_Deliv":0.10,
        "S_Immuno":0.05,"S_Prog":0.30,"S_Mature":0.10
    },
}

SENTINEL_IDS = {"evoCAST", "NlovFz2", "enNlovFz2", "MmeFz2", "SleepingBeauty"}


def penscore_uc(row: pd.Series, uc: str) -> float:
    weights = USE_CASE_WEIGHTS[uc]
    num, den = 0.0, 0.0
    for ax, w in weights.items():
        v = row.get(ax)
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            num += w * v
            den += w
    return round(num / den, 4) if den > 0 else float("nan")


def fmt_score(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "-"
    return f"{v:.4f}"


def score_row(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "N/A"
    if v >= 0.75:
        return f"{v:.4f} +"
    if v >= 0.40:
        return f"{v:.4f} ~"
    return f"{v:.4f} -"


def reasoning_bullets(row: pd.Series, uc: str) -> list[str]:
    weights = USE_CASE_WEIGHTS[uc]
    bullets = []
    for ax in AXES:
        v = row.get(ax)
        is_null = v is None or (isinstance(v, float) and np.isnan(v))
        label = LABELS[ax]
        w = weights.get(ax, 0)
        if is_null:
            bullets.append(f"  - **{label}**: not computed (sentinel or not publicly sequenced)")
        elif v >= 0.75:
            bullets.append(f"  - **{label}**: {v:.4f} - strength (w={w:.2f})")
        elif v <= 0.35:
            bullets.append(f"  - **{label}**: {v:.4f} - weakness (w={w:.2f})")
    return bullets


def top_k_uc(df: pd.DataFrame, uc: str, k: int = 5, dsb_free: bool = False) -> pd.DataFrame:
    sub = df.copy()
    if dsb_free:
        sub = sub[sub["S_DSB"].notna() & (sub["S_DSB"] >= 0.85)]
    sub = sub.copy()
    sub["_ps"] = sub.apply(lambda r: penscore_uc(r, uc), axis=1)
    return sub.nlargest(k, "_ps").reset_index(drop=True)


def editor_table(rows: pd.DataFrame, uc: str) -> str:
    weights = USE_CASE_WEIGHTS[uc]
    header = "| Rank | Editor | PenScore | DSB | Spec | Cargo | Deliv | Immuno | Prog | Mature | Note |\n"
    header += "|------|--------|----------|-----|------|-------|-------|--------|------|--------|------|\n"
    lines = [header]
    for i, (_, row) in enumerate(rows.iterrows(), 1):
        ps = penscore_uc(row, uc)
        note = "provisional (sentinel)" if row["editor_id"] in SENTINEL_IDS else ""
        lines.append(
            f"| {i} | **{row['editor_id']}** | {fmt_score(ps)} "
            f"| {score_row(row.get('S_DSB'))} "
            f"| {score_row(row.get('S_Spec'))} "
            f"| {score_row(row.get('S_Cargo'))} "
            f"| {score_row(row.get('S_Deliv'))} "
            f"| {score_row(row.get('S_Immuno'))} "
            f"| {score_row(row.get('S_Prog'))} "
            f"| {score_row(row.get('S_Mature'))} "
            f"| {note} |\n"
        )
    return "".join(lines)


# ---------------------------------------------------------------------------
# Lab profiles
# ---------------------------------------------------------------------------

LABS = [
    {
        "slug":  "hsu_arc",
        "title": "Hsu Lab (Arc Institute) - Site-Specific Integration & AAV Delivery",
        "pi":    "Patrick Hsu, Arc Institute",
        "focus_summary": (
            "The Hsu Lab pioneers RNA-guided genome editing with a focus on "
            "site-specific genomic integration and AAV-based delivery. Their work on "
            "CAST transposases (Tn7-like RNA-guided systems) and compact nucleases "
            "makes them a natural audience for DSB-free integrase scoring."
        ),
        "primary_uc": "human_therapeutic_aav_insertion",
        "dsb_free": True,
        "lab_editors": ["CAST_IF", "CAST_VK", "evoCAST", "IS621", "IS621_2"],
        "collab_angle": (
            "The Hsu Lab has characterized CAST systems extensively. PEN-SCORE ranks "
            "IS621 as the top DSB-free editor for AAV therapeutic insertion, with "
            "CAST_IF and evoCAST in the top-5. We invite discussion on co-authorship "
            "in a PEN-ASSEMBLE benchmarking study that would compare intracellular "
            "integration efficiency to PEN-SCORE predictions."
        ),
        "collaboration_opportunity": (
            "We propose a prospective wet-lab validation: measure IS621 vs evoCAST "
            "integration efficiency at a safe-harbour locus (AAVS1 or Rosa26 ortholog) "
            "in HEK293T cells and correlate with PenScore differences. If IS621 "
            "outperforms evoCAST in functional assays as predicted (rank 1 vs rank 3 "
            "in our framework), this would constitute the first prospective experimental "
            "validation of the PEN-SCORE composite metric."
        ),
    },
    {
        "slug":  "liu_broad",
        "title": "Liu Lab (Broad Institute) - Base Editing & Prime Editing",
        "pi":    "David Liu, Broad Institute / Harvard",
        "focus_summary": (
            "The Liu Lab invented base editing (CBEs, ABEs) and prime editing (PE2, "
            "PE5max, TwinPE). PEN-SCORE includes five Liu-lab editors: BE3, ABE7.10, "
            "PE2, PE5max, and TwinPE. These tools score highly on specificity but face "
            "strong DSB-safety penalties (S_DSB=0.0 - all are SpCas9-based nucleases "
            "that create nicks or DSBs) and poor deliverability (S_Deliv=0.09 - >5 kb "
            "combined payload)."
        ),
        "primary_uc": "base_editing_small_correction",
        "dsb_free": False,
        "lab_editors": ["BE3", "ABE7_10", "PE2", "PE5max", "TwinPE"],
        "collab_angle": (
            "Our analysis identifies a therapeutic window where base editors (ABE7.10, "
            "BE3) are optimal for small precise corrections but are outranked by "
            "DSB-free integrases for any insertion task. We believe this quantitative "
            "trade-off analysis - with specific axis scores - could be useful for the "
            "Liu Lab's ongoing therapeutic platform development and potential "
            "regulatory submissions."
        ),
        "collaboration_opportunity": (
            "We propose a head-to-head comparison: prime editing (PE5max) vs IS621 "
            "for a correction task that both can technically perform (e.g., SERPINA1 "
            "single-nucleotide correction). PenScore predicts PE5max at 0.34 vs IS621 "
            "at 0.93 for the AAV therapeutic use case. An experimental test of "
            "editing efficiency, off-target rate, and immunogenic response would "
            "directly validate or challenge our composite weighting scheme."
        ),
    },
    {
        "slug":  "zhang_broad",
        "title": "Zhang Lab (Broad Institute) - CRISPR Discovery & Miniaturized Nucleases",
        "pi":    "Feng Zhang, Broad Institute / McGovern",
        "focus_summary": (
            "The Zhang Lab is responsible for founding discoveries in CRISPR-Cas9 "
            "application, Cas12, and IS110-family miniaturized nucleases (IscB, IsrB). "
            "Notably, IscB - a compact RNA-guided recombinase-nuclease - scores "
            "dramatically higher than SpCas9 in the PEN-SCORE framework (0.78 vs 0.37) "
            "primarily because of its superior deliverability (S_Deliv=0.92 vs 0.09). "
            "Cas12f scores similarly well on delivery (0.86). These findings may be "
            "relevant to the Zhang Lab's ongoing miniaturized-editor pipeline."
        ),
        "primary_uc": "research_discovery",
        "dsb_free": False,
        "lab_editors": ["SpCas9", "Cas12a", "Cas12f", "IscB", "SpuFz1", "SpuFz1_V4"],
        "collab_angle": (
            "IscB is classified in PEN-SCORE as DSB_FREE_TRANSEST_RECOMBINASE (IS110 "
            "family), not as a classical nuclease - giving it S_DSB=0.9 and placing "
            "it in the top-12 DSB-free editors. If our mech-class classification is "
            "correct, IscB's therapeutic potential has been underestimated. We seek "
            "biochemical input from the Zhang Lab on whether IscB's transposon-derived "
            "mechanism truly avoids DSBs in human cells."
        ),
        "collaboration_opportunity": (
            "Proposed collaboration: genome-wide DSB mapping (GUIDE-seq or CHANGE-seq) "
            "for IscB, Cas12f, and SpCas9 at matched target sites. PEN-SCORE predicts "
            "IscB >> Cas12f >> SpCas9 on combined DSB-safety + deliverability. "
            "Correlation with measured off-target rates and integration frequency "
            "would validate the S_DSB and S_Spec axis weighting."
        ),
    },
    {
        "slug":  "sternberg_columbia",
        "title": "Sternberg Lab (Columbia) - CRISPR Mechanisms & RNA-Guided Transposition",
        "pi":    "Samuel Sternberg, Columbia University",
        "focus_summary": (
            "The Sternberg Lab is the leading group on mechanistic understanding of "
            "CAST (CRISPR-associated transposase) systems, particularly Tn7-like "
            "RNA-guided transposases that insert large DNA cargos without creating "
            "double-strand breaks. PEN-SCORE's large-cargo use case is designed "
            "precisely for this class of editors. CAST_IF (0.80) and CAST_VK (0.78) "
            "both rank highly; evoCAST (0.87 in AAV context, P1 prediction: PASS) "
            "is the top CAST-family system overall."
        ),
        "primary_uc": "large_cargo_integration",
        "dsb_free": True,
        "lab_editors": ["CAST_IF", "CAST_VK", "evoCAST", "Tn5"],
        "collab_angle": (
            "The Sternberg Lab's detailed biochemical characterization of CAST "
            "mechanisms is directly relevant to PEN-SCORE's S_DSB and S_Prog axes. "
            "We would welcome expert review of our IS110-class classification criteria "
            "(mech-class v0.5.1) and the S_DSB bucket heuristic for Tn7-like systems. "
            "Inaccuracies in classification propagate to all downstream therapeutic "
            "rankings - external validation is essential."
        ),
        "collaboration_opportunity": (
            "We propose a targeted collaboration: measure large-cargo insertion "
            "efficiency (10 kb, 50 kb, 100 kb payloads) for CAST_IF vs evoCAST vs "
            "IS621 in human cells. PEN-SCORE predicts IS621 > evoCAST > CAST_IF for "
            "the large-cargo profile (S_Cargo: 1.0 vs 0.78 vs 0.67). A discrepancy "
            "between predicted and measured efficiency would identify axes requiring "
            "recalibration - particularly whether programmability (S_Prog) should "
            "be weighted more heavily for Tn7-class systems."
        ),
    },
    {
        "slug":  "sjnahs_vit",
        "title": "SJNAHS / VIT - Internal Summary & PEN-ASSEMBLE Roadmap",
        "pi":    "Anees Ahmed Mahaboob Ali, SJNAHS / VIT",
        "focus_summary": (
            "This document summarises the PEN-SCORE public scorecard (v0.0.1, PEN-SCORE) "
            "for internal use at SJNAHS/VIT. It serves as the foundation for PEN-ASSEMBLE "
            "(prospective experimental validation) and as a reference for grant "
            "applications, ethics submissions, and wet-lab planning."
        ),
        "primary_uc": "human_therapeutic_aav_insertion",
        "dsb_free": False,
        "lab_editors": None,   # all editors
        "collab_angle": (
            "PEN-ASSEMBLE will test the top PEN-SCORE predictions experimentally. The key "
            "experiments are: (1) IS621 vs SpCas9 integration efficiency at AAVS1 in "
            "HEK293T cells; (2) IS621 vs evoCAST cargo capacity at 5 kb / 20 kb; "
            "(3) MHC epitope validation for IS621 vs SpCas9 immunogenic peptides "
            "predicted by NetMHCpan in HLA-A*02:01 / HLA-B*07:02 donors."
        ),
        "collaboration_opportunity": (
            "Pre-registered experimental plan (to be submitted to bioRxiv with "
            "this paper as pre-registration): three-arm cell assay, n=3 biological "
            "replicates, Sanger + NGS validation, ELISA immunogenicity readout. "
            "Target: submit PEN-ASSEMBLE within 18 months of PEN-SCORE publication."
        ),
    },
]


# ---------------------------------------------------------------------------
# Generate Markdown
# ---------------------------------------------------------------------------

def generate_md(lab: dict, df: pd.DataFrame) -> str:
    uc = lab["primary_uc"]
    uc_label = uc.replace("_", " ").title()

    # Relevant editors
    if lab["lab_editors"] is not None:
        relevant = df[df["editor_id"].isin(lab["lab_editors"])].copy()
    else:
        relevant = df.copy()

    # Top-5 recommendations
    top5 = top_k_uc(df, uc, k=5, dsb_free=lab["dsb_free"])

    lines = []
    lines.append(f"# PEN-SCORE Outreach: {lab['title']}\n")
    lines.append(f"**PI / Contact:** {lab['pi']}  \n")
    lines.append(f"**Primary Use Case:** {uc_label}  \n")
    lines.append(f"**Generated:** 2026-05-14 by `scripts/41_generate_outreach_materials.py`  \n")
    lines.append(f"**PEN-SCORE version:** v0.0.1 (mech-class v0.5.1)  \n\n")
    lines.append("---\n\n")

    # (a) Lab focus
    lines.append("## A. Lab Focus & PEN-SCORE Relevance\n\n")
    lines.append(lab["focus_summary"] + "\n\n")

    # (b) Editors relevant to this lab
    if lab["lab_editors"] is not None and len(relevant) > 0:
        lines.append("## B. Editors From This Lab in the Public Scorecard\n\n")
        lines.append("The following editors developed or characterised by this group are included "
                     "in the PEN-SCORE v0.0.1 public scorecard (28 editors total):\n\n")
        lines.append(editor_table(relevant, uc))
        lines.append("\n*Score legend: + strength (>=0.75) | ~ acceptable (>=0.40) | - weakness (<0.40) | . not computed*\n\n")
    else:
        lines.append("## B. Full Scorecard Overview\n\n")
        all_sorted = top_k_uc(df, uc, k=28, dsb_free=False)
        lines.append(editor_table(all_sorted, uc))
        lines.append("\n*Score legend: + strength (>=0.75) | ~ acceptable (>=0.40) | - weakness (<0.40) | . not computed*\n\n")

    # (c) Top-5 recommendations with reasoning
    dsb_note = " (DSB-free only)" if lab["dsb_free"] else ""
    lines.append(f"## C. Top-5 Recommended Editors for *{uc_label}*{dsb_note}\n\n")

    for i, (_, row) in enumerate(top5.iterrows(), 1):
        ps = penscore_uc(row, uc)
        sentinel_note = " *(provisional - no public sequence)*" if row["editor_id"] in SENTINEL_IDS else ""
        lines.append(f"### {i}. {row['editor_id']} - PenScore {fmt_score(ps)}{sentinel_note}\n\n")
        bullets = reasoning_bullets(row, uc)
        if bullets:
            lines.append("**Axis-by-axis reasoning:**\n")
            lines.extend([b + "\n" for b in bullets])
        else:
            lines.append("All axes within acceptable range.\n")
        lines.append("\n")

    # (d) Collaboration invitation
    lines.append("## D. Collaboration Angle\n\n")
    lines.append(lab["collab_angle"] + "\n\n")

    lines.append("## E. PEN-ASSEMBLE Collaboration Opportunity\n\n")
    lines.append(lab["collaboration_opportunity"] + "\n\n")

    lines.append("---\n\n")
    lines.append("## Contact\n\n")
    lines.append("**Anees Ahmed Mahaboob Ali**  \n")
    lines.append("SJNAHS / VIT  \n")
    lines.append("ahmedaneesm@gmail.com  \n")
    lines.append("GitHub: [ahmedanees-m/pen-score](https://github.com/ahmedanees-m/pen-score)  \n\n")
    lines.append("*This document was generated automatically from the PEN-SCORE v0.0.1 public "
                 "scorecard. All scores are pre-registered (pre-registration-v1.0.2, "
                 "locked 2026-05-13T17:56:38Z) and have not been manually tuned.*\n")

    return "".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

generated = []
for lab in LABS:
    md = generate_md(lab, df)
    out_path = OUTREACH_DIR / f"{lab['slug']}.md"
    out_path.write_text(md, encoding="utf-8")
    generated.append((lab["slug"], len(md)))
    print(f"Written: {out_path}  ({len(md):,} chars)")

print(f"\nTotal: {len(generated)} files written to {OUTREACH_DIR}")
print("Files:", [g[0] for g in generated])
