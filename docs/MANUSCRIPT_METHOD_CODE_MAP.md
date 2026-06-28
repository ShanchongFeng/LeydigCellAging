# Manuscript-method-code map

Audit target: `Manuscript_v21.docx` and `Tables.xlsx`, frozen 2026-06-28.

| Results | Methods | Primary code | Audit status |
|---|---|---|---|
| 2.1 bulk perturbation | 4.2 | `1_bulk/` | Present and aligned |
| 2.1 mouse state localization | 4.3-4.4 | `2_mouse/` | Present and aligned; group-level only |
| 2.1 topology diagnostic | 4.3 | `2_mouse/4_topology.R` | Present and correctly descriptive |
| 2.2 GSE303193 | 4.5, 4.11-4.12 | `4_mouse_external/`, `7_robustness/` | Present and aligned |
| 2.2 GSE254315 primary | 4.6, 4.11-4.12 | `5_human_targeted/`, `7_robustness/` | Downstream code present; raw-to-clustered object construction unavailable |
| 2.2 GSE254315 composition | 4.6, 4.15 | `13_gse254315_composition/` | Code added; Methods need the 25-cell and regression details |
| 2.2 GSE182786 target test | 4.8, 4.10-4.12 | `6_human/`, `7_robustness/` | Present; Results should name exact permutation P values |
| 2.3 functional categories | 4.13 | `10_null/1_orthologs.R`, `12_functional_categories/` | Code added; effect is logCPM difference, not standardized |
| 2.4 hdWGCNA/modules | 4.7, 4.9-4.10 | `3_hdwgcna/`, `6_human/module_checks/` | Present; one Results sentence does not match S10/S11 directions |
| 2.4 portable signature/null | 4.16 | `9_signature/`, `10_null/` | Present and aligned |
| 2.5 scFEA | 4.17 | `11_scfea/` | Present; 164/161 finite tests, not 168/168 |
| 2.6 robustness | 4.10-4.15 | `6_human/`, `7_robustness/`, `8_state_composition/` | Present; GSE254315-specific composition procedure needed in Methods |

The main runner executes 49 primary steps. The early exploratory curated
z-score probe is retained but excluded from the primary run order.
