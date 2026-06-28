# Functional-category analyses

This directory contains the curated, Gene Ontology (GO), and Reactome
category analyses used by manuscript Results 2.3 and Methods 4.13.

## Primary workflow

1. `2_curated_gene_effects.R` summarizes fixed curated sets from the
   cross-dataset gene-effect table.
2. `3_prepare_go_gene_sets.py` constructs GO sets from fixed roots and
   `is_a`/`part_of` descendants.
3. `4_go_category_effects.R` summarizes GO category effects.
4. `5_prepare_reactome_gene_sets.py` constructs Reactome sets from fixed
   pathway identifiers.
5. `6_reactome_category_effects.R` summarizes Reactome category effects.
6. `7_union_leave_one_out.R` performs leave-one-gene-out checks for the
   prespecified metabolic-support unions.

`1_curated_score_probe.R` is an earlier z-score-based exploratory probe. It
is retained for provenance but is not part of `RUN_ORDER.tsv` and is not the
source of the manuscript's primary category estimates.

## Effect definition

The primary input `ortholog_gene_effects.csv` contains donor/sample-level
aged-minus-young log2(CPM + 0.5) differences. Category effects are medians of
those gene effects, first within each dataset-boundary block, then within each
dataset across its shared CMB/ESB boundaries. They are not standardized
effects.

GO and Reactome are two annotation systems applied to the same three dataset
blocks. They are not independent biological replications.

## Required inputs

- `ORTHOLOG_GENE_EFFECTS_FILE`
- `ORTHOLOG_ONE_TO_ONE_FILE`
- `CURATED_GENE_SET_MANIFEST`
- `FUNCTIONAL_CATEGORY_ANNOTATION_ROOT` containing the four files listed in
  `ANNOTATION_INPUTS.tsv`

All outputs are written below `FUNCTIONAL_CATEGORY_ROOT`.
