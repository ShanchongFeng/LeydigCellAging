# Results, Methods, and code audit

## Frozen snapshot

- Manuscript: `Manuscript_v21.docx`, SHA-256
  `c5e2c054afd76990f16b31ed24f1cd47e4ce3a8cd70c10ab7599c2d051734cea`
- Supplement: `Tables.xlsx`, SHA-256
  `70636fdd5948dd1fe801148cd34ff7bc879cde40b844fc41cbe57860d9fddc71`
- GitHub starting commit: `93824bba31e9133a967f841beaf94ee52b1d78b4`
- Figure appearance was outside this audit.

## Overall conclusion

The main donor/sample-aware framework is largely traceable, but v21 is not
yet methodologically self-consistent. Two missing analysis families have now
been added to the repository. One upstream provenance gap remains: the exact
GSE254315 raw-to-clustered-object construction script is unavailable.

## Findings requiring manuscript revision

### 1. The HMGCS2 effect is not standardized

The abstract and Methods 4.13 call the value -1.134 a "signed standardized
effect". The code uses raw aged-minus-young log2(CPM + 0.5) differences. For
HMGCS2, the two boundaries are first summarized within each dataset and the
three dataset values are then summarized by their median.

Use: `dataset-blocked median aged-minus-young logCPM difference = -1.134`.
Methods 4.13 should define the boundary, dataset, and cross-dataset medians.

### 2. The GSE254315 object cannot be rebuilt from the repository

Methods 4.6 describes donor-level 10x import, QC, normalization, reduction,
and clustering. The repository starts the formal analysis from a locked RDS.
A full server search found only a superseded script using age <=45 and >=55,
not the manuscript grouping of <=44 and >=52. That script must remain
excluded.

Preferred resolution: independently reconstruct the object from raw data,
freeze parameters, and verify cluster identities and all downstream values.
Interim resolution: disclose that downstream analyses use a retained processed
object and state this provenance limitation in Code Availability.

### 3. The GSE254315 composition conclusion is overgeneralized

Results 2.2 says steroidogenic execution weakens after composition adjustment
at both boundaries. Table S21 does not support a general statement:

- HMGCS2 becomes more negative in CMB (-0.980 to -1.368) and ESB (-1.289 to
  -1.766), but neither survives the regression-family BH correction (q=0.390).
- GO steroidogenic execution weakens only in ESB (+0.321, P=0.047 to +0.209,
  P=0.270).
- In CMB it becomes stronger (+0.210 to +0.304), and curated steroidogenesis
  does not show consistent attenuation.

Recommended wording: report the ESB-specific attenuation and state that the
steroidogenic composition sensitivity is readout- and boundary-dependent.

Methods 4.15 also needs a dedicated GSE254315 paragraph: equal-cluster primary
threshold 25 cells per donor-cluster (only 8 complete donors), descriptive
thresholds 1/50/100, and all-23-donor age-plus-nonreference-cluster-proportion
linear models with cluster 1 as reference. The current generic 50/100-cell
description does not cover Table S21.

### 4. One module/program Results sentence has the wrong direction

Results 2.4 states that CMB metabolic-support and steroidogenic-related
readouts are aged-lower. In Table S11, CMB ketogenesis/FAO is -0.153, but CMB
steroidogenesis is +0.142. Table S10 also has greenyellow at +0.020 in CMB.

Replace the blanket statement with the actual split directions. The
normal-BMI ESB results (-0.459 steroidogenesis and -0.349 ketogenesis/FAO) can
remain as a stratum-specific finding.

### 5. scFEA multiple-testing wording is numerically wrong and test-sensitive

The export contains 168 modules, but only 164 finite Welch P values and 161
finite Spearman P values enter BH correction. Methods 4.17 and Table S13 must
state those denominators instead of 168 and 168.

An exact two-sided label-permutation audit over the 12 sample means gives
M_169 P=0.00202, but BH across the 164 finite module tests is approximately
0.331. The Welch FDR=0.021 therefore depends on the parametric test choice.
M_169 should remain explicitly exploratory and test-sensitive.

### 6. P values and supporting tables need clearer labels

- Results 2.2 should call the GSE182786 CMB P=0.012 and ESB P=0.002
  two-sided exact label-permutation P values.
- The HMGCS2 -1.134 and HMGCS2-axis leave-out values come from curated
  category outputs but are absent from S18/S19. Add those rows or a dedicated
  supplementary block.
- Supplementary mappings are stale: S6-S9, S10-S13, S14-S17, and S18-S20
  still point to the pre-reordering Results section numbers.

### 7. Annotation and code availability need publication-grade provenance

GO and Reactome are two annotation systems applied to the same datasets, not
independent biological evidence. Preserve that distinction in the abstract
and Results. Deposit the exact four annotation files with fixed hashes in a
versioned archive. Add a Code Availability statement with the repository URL,
release or commit, environment requirements, and the GSE254315 limitation.

## Repository remediation in this audit

- Added the functional-category scripts and fixed their hard-coded server
  paths.
- Added GSE254315 composition and adjusted-regression scripts.
- Added fixed annotation URLs/hashes, section READMEs, method-code mapping,
  provenance documentation, and steps 42-49 to the runner.
- Retained the exploratory curated score probe but excluded it from the
  manuscript run order.
- Did not add the incompatible legacy GSE254315 reconstruction script.

## Verification

- All R files parsed, both new Python files compiled, and both shell runners
  passed `bash -n`.
- An isolated server rerun reproduced all 18 functional-category tables
  byte-for-byte.
- The GSE254315 composition rerun reproduced all tables byte-for-byte except
  the regression CSV; its maximum numeric difference was `1.998401e-14`
  (`age_t`), consistent with floating-point linear algebra variation.
