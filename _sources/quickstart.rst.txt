Quick Start
===========

Installation
------------

.. code-block:: bash

    pip install pen-score

With machine-learning extras (LightGBM, scikit-learn):

.. code-block:: bash

    pip install "pen-score[ml]"

Python API
----------

Score the full editor universe for a clinical use case:

.. code-block:: python

    from pen_score.api import Scorer

    scorer = Scorer.load()

    # Top-5 editors for AAV-based human therapeutics
    ranked = scorer.select_editor(
        use_case="human_therapeutic_aav_insertion",
        top_k=5,
    )
    print(ranked[["editor_id", "PenScore"]].to_string(index=False))
    # editor_id  PenScore
    # IS621        0.9290
    # enNlovFz2    0.8891
    # NlovFz2      0.8759
    # evoCAST      0.8620
    # MmeFz2       0.8404

DSB-free filtering:

.. code-block:: python

    dsb_free = scorer.select_editor(
        use_case="human_therapeutic_aav_insertion",
        top_k=5,
        require_dsb_free=True,   # S_DSB >= 0.85
    )

Biophysical-only score (exclude S_Mature for fair comparison of new editors):

.. code-block:: python

    # IS622 has S_Mature=0 because it is brand-new (2026).
    # Excluding S_Mature gives a fair biophysical comparison.
    result = scorer.score_editor(
        accession="D2TGM5",  # IS622 / ISCro4
        use_case="human_therapeutic_aav_insertion",
        exclude_axes=["S_Mature"],
    )
    print(result.pen_score)       # >0.93 biophysical PenScore
    print(result.axes.S_Mature)   # 0.0 - still in axes output, just not in composite

CLI
---

.. code-block:: bash

    # Top-5 for AAV insertion
    pen-score select --use-case human_therapeutic_aav_insertion --top-k 5

    # DSB-free only, large cargo
    pen-score select --use-case large_cargo_integration --top-k 5 --require-dsb-free

    # Score a single editor
    pen-score score-editor Q99ZW2 --use-case human_therapeutic_aav_insertion

    # Score IS622 excluding S_Mature (biophysical comparison)
    pen-score score-editor D2TGM5 --exclude-axes S_Mature

Use cases
---------

+---------------------------------------------+-----------------------------------------------+
| Key                                         | Description                                   |
+=============================================+===============================================+
| ``human_therapeutic_aav_insertion``         | AAV-compatible human gene therapy             |
+---------------------------------------------+-----------------------------------------------+
| ``human_therapeutic_electroporation``       | Ex vivo electroporation; high specificity     |
+---------------------------------------------+-----------------------------------------------+
| ``large_cargo_integration``                 | Payload > 4.5 kb (lentiviral / LNP)           |
+---------------------------------------------+-----------------------------------------------+
| ``base_editing_small_correction``           | Small point-mutation correction               |
+---------------------------------------------+-----------------------------------------------+
| ``research_discovery``                      | Academic discovery; maturity less important   |
+---------------------------------------------+-----------------------------------------------+
| ``megabase_rearrangement``                  | >100 kb inversions / translocations (v0.1.2)  |
+---------------------------------------------+-----------------------------------------------+
| ``therapeutic_excision_bcl11a``             | BCL11A +58 enhancer excision, SCD / β-thal   |
+---------------------------------------------+-----------------------------------------------+

Interactive scorecard
---------------------

Browse all 31 curated editors with live use-case switching in the
`interactive scorecard browser <scorecards/index.html>`_.
