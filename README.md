# Leydig Cell Aging

This repository contains analysis scripts for a multi-cohort transcriptomic study of Leydig cell aging. The scripts cover bulk RNA-seq, mouse and human single-cell RNA-seq, hdWGCNA, cross-dataset robustness checks, functional-category analyses, state and cluster-composition sensitivity, portable signature tests, a size-matched random gene-set null analysis, and a scFEA-based Leydig CMB flux supplement.

## Scope

The `Analysis Script/` directory contains:

1. analysis scripts supporting the methods and results described in the manuscript;
2. the ordered workflow wrapper and run manifest;
3. file checksums for the analysis scripts.

Input data are not included in this repository.

## Reproduction Command

`Analysis Script/RUN_ORDER.tsv` records the 49 executable primary analysis steps, their environments, dependencies, and primary output locations. `Analysis Script/run_all.sh` creates an isolated run tree, links the required inputs read-only, and writes all new outputs and logs beneath:

`<SUBMISSION_RUN_ROOT>` or, if unset, `<MANUSCRIPT_INPUT_ROOT>/submission_reproduction_20260605`

Set the site-specific path placeholders before running:

| Placeholder | Meaning |
|---|---|
| `MANUSCRIPT_INPUT_ROOT` | Read-only input root containing `raw_data/`, `rds/`, and `tables/` |
| `LOCKED_INPUT_ROOT` | Linked copy of locked input files used by downstream scripts |
| `LEGACY_INPUT_ROOT` | Root containing cached GEO archives such as `GSE303193_RAW.tar` and `GSE254315_RAW.tar`; defaults to `<MANUSCRIPT_INPUT_ROOT>/legacy_inputs` |
| `GSE303193_RAW_TAR` | Optional direct path to `GSE303193_RAW.tar`, overriding `LEGACY_INPUT_ROOT` |
| `GSE254315_RAW_TAR` | Optional direct path to `GSE254315_RAW.tar`, overriding `LEGACY_INPUT_ROOT` |
| `SUBMISSION_RUN_ROOT` | Optional isolated output/log root |
| `SUBMISSION_WORK_ROOT` | Optional isolated working tree root; defaults to `<SUBMISSION_RUN_ROOT>/work` |
| `BULK_GSE287203_ROOT` | Bulk GSE287203 DESeq2 and enrichment output root |
| `GSE270931_SCRNA_ROOT` | Mouse GSE270931 scRNA processing, state, and Harmony-control root |
| `HDWGCNA_ROOT` | Mouse Leydig hdWGCNA reconstruction and module-enrichment root |
| `GSE303193_EXTERNAL_MOUSE_ROOT` | External mouse GSE303193 reconstruction and boundary-sensitivity root |
| `GSE303193_PSEUDOBULK_ROOT` | GSE303193 sample-level pseudobulk output root |
| `GSE254315_HUMAN_TARGETED_ROOT` | Human GSE254315 acquisition and metadata-inspection root |
| `GSE254315_PSEUDOBULK_ROOT` | GSE254315 donor-level pseudobulk output root |
| `GSE182786_HUMAN_AGING_ROOT` | Human GSE182786 acquisition and Seurat reconstruction root |
| `GSE182786_BOUNDARY_STRESS_ROOT` | GSE182786 boundary and BMI-sensitivity root |
| `GSE182786_SCHEME3_ROOT` | GSE182786 module-summary, LOO, threshold, and sparsity root |
| `DONOR_SAMPLE_DOWNSAMPLING_ROOT` | Donor/sample-balanced downsampling root |
| `CLUSTER_BOUNDARY_AUDIT_ROOT` | Cluster-boundary and cross-dataset check root |
| `STATE_STRUCTURE_ROOT` | State-composition, topology, and signature-state output root |
| `PORTABLE_SIGNATURE_PSEUDOBULK_ROOT` | Cross-dataset signature pseudobulk root |
| `ORTHOLOG_CONSERVATION_ROOT` | HomoloGene ortholog-conservation root |
| `ENSEMBL_ORTHOLOG_ROOT` | Ensembl release-115 ortholog retrieval root |
| `PORTABLE_SIGNATURE_AUDIT_ROOT` | HomoloGene signature check root |
| `ENSEMBL_SIGNATURE_AUDIT_ROOT` | Ensembl signature check root |
| `RANDOM_SIGNATURE_NULL_ROOT` | Size-matched random signature null root |
| `SCFEA_REPO_ROOT` | Local scFEA repository root containing `src/scFEA.py` and `data/` |
| `SCFEA_OUTPUT_ROOT` | scFEA input, raw flux, summary table, figure, and log output root |
| `SCFEA_PYTHON_BIN` | Python executable with scFEA dependencies; defaults to `PYTHON_BIN` |
| `SCFEA_MAX_CELLS_PER_SAMPLE` | Maximum Leydig CMB cells retained per sample before scFEA input export; defaults to `400` |
| `SCFEA_EPOCHS` | scFEA training epochs; defaults to `100` |
| `FUNCTIONAL_CATEGORY_ROOT` | Curated, GO, and Reactome category output root |
| `FUNCTIONAL_CATEGORY_ANNOTATION_ROOT` | Directory containing the four fixed GO/Reactome annotation inputs |
| `CURATED_GENE_SET_MANIFEST` | Frozen curated gene-set manifest |
| `ORTHOLOG_GENE_EFFECTS_FILE` | Cross-dataset donor/sample-level gene-effect table |
| `GSE254315_COMPOSITION_ROOT` | GSE254315 composition-sensitivity output root |
| `GSE254315_LOCKED_RDS` | Frozen manuscript-aligned GSE254315 Seurat object |
| `CONDA_EXE` | Conda executable; defaults to `conda` on `PATH` |
| `PYTHON_BIN` | Python executable for the Python summary scripts; defaults to `python` on `PATH` |

Run with:

```bash
MANUSCRIPT_INPUT_ROOT="<MANUSCRIPT_INPUT_ROOT>" \
LEGACY_INPUT_ROOT="<LEGACY_INPUT_ROOT>" \
SCFEA_REPO_ROOT="<SCFEA_REPO_ROOT>" \
CONDA_EXE="<CONDA_EXE>" \
PYTHON_BIN="<PYTHON_BIN>" \
SCFEA_PYTHON_BIN="<SCFEA_PYTHON_BIN>" \
bash "Analysis Script/run_all.sh"
```

To resume after a repaired failed step while retaining earlier outputs:

```bash
MANUSCRIPT_INPUT_ROOT="<MANUSCRIPT_INPUT_ROOT>" \
LEGACY_INPUT_ROOT="<LEGACY_INPUT_ROOT>" \
START_AT=18 \
bash "Analysis Script/run_all.sh"
```

The execution status table is written to `<SUBMISSION_RUN_ROOT>/run_status.tsv`; individual logs are written to `<SUBMISSION_RUN_ROOT>/logs/`.

## Directory Map

| Directory | Manuscript analysis covered |
|---|---|
| `Analysis Script/1_bulk` | Bulk RNA-seq differential expression and pathway enrichment |
| `Analysis Script/2_mouse` | Primary mouse scRNA-seq processing, Leydig-state analysis, Harmony control, and topology checks |
| `Analysis Script/3_hdwgcna` | Mouse Leydig hdWGCNA network construction and module enrichment |
| `Analysis Script/4_mouse_external` | External mouse dataset reconstruction, boundary sensitivity, and sample-level pseudobulk |
| `Analysis Script/5_human_targeted` | Human targeted dataset acquisition, metadata inspection, and donor-level pseudobulk |
| `Analysis Script/6_human` | Human testis dataset processing, boundary checks, BMI downsampling, and module robustness analyses |
| `Analysis Script/7_robustness` | Donor/sample downsampling, cluster-boundary checks, and cross-dataset robustness summaries |
| `Analysis Script/8_state_composition` | Age-associated Leydig-state composition analysis |
| `Analysis Script/9_signature` | Pseudobulk export, portable signature scoring, and ortholog-based signature checks |
| `Analysis Script/10_null` | Ortholog universe construction and size-matched random signature null analysis |
| `Analysis Script/11_scfea` | GSE182786 Leydig CMB scFEA input, flux inference, and sample-level module summary |
| `Analysis Script/12_functional_categories` | Curated, GO, Reactome, dataset-blocked, and leave-one-gene-out category analyses |
| `Analysis Script/13_gse254315_composition` | Donor-cluster composition, equal-cluster scores, and composition-adjusted regressions |

The Results-Methods-code map and the 2026-06-28 audit are in `docs/`.

## Known Provenance Limitation

The exact script that created the March 2026 manuscript-aligned GSE254315 clustered object was not recovered after local and server searches. The locked object contains 12 Young donors (age <=44) and 11 Aged donors (age >=52), and the repository reproduces all downstream analyses from that object. A later script using incompatible age cutoffs was excluded. See `docs/INPUTS_AND_PROVENANCE.md`.

## Data Availability

The public datasets analyzed in this study are available from the NCBI Gene Expression Omnibus under accession numbers GSE287203, GSE270931, GSE303193, GSE254315, and GSE182786. Summary tables supporting the reported analyses are provided with the manuscript as Supplementary Tables. Analysis scripts are provided in this repository under `Analysis Script/`.

## Checksums

`Analysis Script/SHA256SUMS.txt` records SHA-256 hashes using paths relative to `Analysis Script/`. To verify the files after download:

```bash
cd "Analysis Script"
python3 - <<'PY'
from pathlib import Path
import hashlib

expected = {}
for line in Path("SHA256SUMS.txt").read_text().splitlines():
    if not line.strip():
        continue
    digest, name = line.split(maxsplit=1)
    expected[name] = digest

bad = [
    name
    for name, digest in expected.items()
    if not Path(name).is_file() or hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest
]
if bad:
    raise SystemExit("Checksum mismatch: " + ", ".join(bad))
print("Checksums verified")
PY
```

## License

This repository is released under the MIT License.
