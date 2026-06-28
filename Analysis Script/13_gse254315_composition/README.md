# GSE254315 composition sensitivity

These scripts support manuscript Results 2.2 and Supplementary Table S21.

`1_composition_sensitivity.R` uses the frozen CMB (clusters 1, 6, 18) and ESB
(clusters 1, 6, 18, 7) to calculate donor-level cluster proportions,
within-cluster pseudobulk readouts, and equal-cluster donor scores. The primary
equal-cluster threshold is 25 cells per donor-cluster; thresholds 1, 50, and
100 are descriptive sensitivities.

`2_composition_adjusted_regression.R` fits donor-level age-only and
age-plus-composition linear models. Cluster 1 is the reference component and
all other fixed-boundary cluster proportions are included as covariates. The
regression uses all 23 donors.

Persistence after adjustment is an exclusionary sensitivity result. It does
not separate confounding from age-mediated composition change and does not
establish a cell-intrinsic mechanism.

## Provenance limitation

The analysis requires `GSE254315_LOCKED_RDS`, which contains 12 Young donors
(age <=44) and 11 Aged donors (age >=52). The original script that created
this March 2026 clustered object could not be recovered. A later script with
different age cutoffs (<=45 and >=55) is intentionally excluded because it
does not reproduce the manuscript cohort.
