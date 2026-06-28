# Inputs and provenance

## Public datasets

The workflow uses GEO accessions GSE287203, GSE270931, GSE303193, GSE254315,
and GSE182786. Raw data and large processed objects are not stored in Git.

## Frozen analysis inputs

`run_all.sh` expects `MANUSCRIPT_INPUT_ROOT` to contain `raw_data/`, `rds/`,
and `tables/`. Locked processed objects are linked read-only into the isolated
run tree. This preserves the exact downstream inputs used for the manuscript.

## Functional annotations

The GO and Reactome files used on 2026-06-19 are identified by URL and SHA-256
in `Analysis Script/12_functional_categories/ANNOTATION_INPUTS.tsv`. Because
the upstream `current` files can change, long-term publication reproducibility
requires depositing the four exact files in a versioned archive such as
Zenodo and citing its DOI.

## GSE254315 limitation

The locked GSE254315 object and all downstream donor-level tables are retained,
but the manuscript-aligned object-construction script was not found after a
full local and server search. The retained legacy script uses incompatible age
cutoffs and is not represented as final code. Until the object is independently
rebuilt and validated, the repository supports exact downstream reproduction
from the locked object, not end-to-end reconstruction from raw GSE254315 files.
