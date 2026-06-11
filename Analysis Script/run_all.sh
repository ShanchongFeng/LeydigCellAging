#!/usr/bin/env bash
set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANUSCRIPT_INPUT_ROOT="${MANUSCRIPT_INPUT_ROOT:-}"
if [[ -z "$MANUSCRIPT_INPUT_ROOT" ]]; then
  echo "Set MANUSCRIPT_INPUT_ROOT to the read-only input root containing raw_data/, rds/, and tables/." >&2
  exit 2
fi
MANUSCRIPT_INPUT_ROOT="$(cd "$MANUSCRIPT_INPUT_ROOT" && pwd)"
RUN_ROOT="${SUBMISSION_RUN_ROOT:-$MANUSCRIPT_INPUT_ROOT/submission_reproduction_20260605}"
START_AT="${START_AT:-}"
SUBMISSION_WORK_ROOT="${SUBMISSION_WORK_ROOT:-$RUN_ROOT/work}"
LEGACY_INPUT_ROOT="${LEGACY_INPUT_ROOT:-$MANUSCRIPT_INPUT_ROOT/legacy_inputs}"
GSE303193_RAW_TAR="${GSE303193_RAW_TAR:-$LEGACY_INPUT_ROOT/04_mouse_external_gse303193/raw_data/GSE303193_RAW.tar}"
GSE254315_RAW_TAR="${GSE254315_RAW_TAR:-$LEGACY_INPUT_ROOT/05_human_gse254315/raw_data/GSE254315_RAW.tar}"

export HOME="$SUBMISSION_WORK_ROOT/home"
export LOCKED_INPUT_ROOT="${LOCKED_INPUT_ROOT:-$SUBMISSION_WORK_ROOT/LOCKED_INPUTS}"
export BULK_GSE287203_ROOT="${BULK_GSE287203_ROOT:-$SUBMISSION_WORK_ROOT/BULK_GSE287203}"
export GSE270931_SCRNA_ROOT="${GSE270931_SCRNA_ROOT:-$SUBMISSION_WORK_ROOT/GSE270931_SCRNA}"
export GSE270931_RAW_ROOT="${GSE270931_RAW_ROOT:-$LOCKED_INPUT_ROOT/raw_data/GSE270931_mouse_scRNA}"
export GSE270931_SCRIPT_ROOT="${GSE270931_SCRIPT_ROOT:-$SUBMISSION_WORK_ROOT/GSE270931_SCRIPT_SOURCE}"
export HDWGCNA_ROOT="${HDWGCNA_ROOT:-$SUBMISSION_WORK_ROOT/HDWGCNA_MOUSE_LEYDIG}"
export GSE270931_HDWGCNA_INPUT_RDS="${GSE270931_HDWGCNA_INPUT_RDS:-$GSE270931_SCRNA_ROOT/rds/GSE270931_filtered_clustered_seurat.rds}"
export HDWGCNA_PRENETWORK_INPUT_RDS="${HDWGCNA_PRENETWORK_INPUT_RDS:-$LOCKED_INPUT_ROOT/rds/primary_locked/hdWGCNA/leydig_scRNA_phase2a__phase2b__hdwgcna__rds__00_input_obj.rds}"
export HDWGCNA_MODULE_ASSIGNMENT_CSV="${HDWGCNA_MODULE_ASSIGNMENT_CSV:-$HDWGCNA_ROOT/tables/hdwgcna_module_assignment.csv}"
export GSE303193_EXTERNAL_MOUSE_ROOT="${GSE303193_EXTERNAL_MOUSE_ROOT:-$SUBMISSION_WORK_ROOT/GSE303193_EXTERNAL_MOUSE}"
export GSE303193_PSEUDOBULK_ROOT="${GSE303193_PSEUDOBULK_ROOT:-$SUBMISSION_WORK_ROOT/GSE303193_PSEUDOBULK}"
export GSE254315_HUMAN_TARGETED_ROOT="${GSE254315_HUMAN_TARGETED_ROOT:-$SUBMISSION_WORK_ROOT/GSE254315_HUMAN_TARGETED}"
export GSE254315_PSEUDOBULK_ROOT="${GSE254315_PSEUDOBULK_ROOT:-$SUBMISSION_WORK_ROOT/GSE254315_PSEUDOBULK}"
export GSE182786_HUMAN_AGING_ROOT="${GSE182786_HUMAN_AGING_ROOT:-$SUBMISSION_WORK_ROOT/GSE182786_HUMAN_AGING}"
export GSE182786_BOUNDARY_STRESS_ROOT="${GSE182786_BOUNDARY_STRESS_ROOT:-$SUBMISSION_WORK_ROOT/GSE182786_BOUNDARY_STRESS}"
export GSE182786_SCHEME3_ROOT="${GSE182786_SCHEME3_ROOT:-$SUBMISSION_WORK_ROOT/GSE182786_SCHEME3}"
export DONOR_SAMPLE_DOWNSAMPLING_ROOT="${DONOR_SAMPLE_DOWNSAMPLING_ROOT:-$SUBMISSION_WORK_ROOT/DONOR_SAMPLE_DOWNSAMPLING}"
export CLUSTER_BOUNDARY_AUDIT_ROOT="${CLUSTER_BOUNDARY_AUDIT_ROOT:-$SUBMISSION_WORK_ROOT/CLUSTER_BOUNDARY_AUDIT}"
export STATE_STRUCTURE_ROOT="${STATE_STRUCTURE_ROOT:-$SUBMISSION_WORK_ROOT/STATE_STRUCTURE}"
export PORTABLE_SIGNATURE_PSEUDOBULK_ROOT="${PORTABLE_SIGNATURE_PSEUDOBULK_ROOT:-$SUBMISSION_WORK_ROOT/PORTABLE_SIGNATURE_PSEUDOBULK}"
export ORTHOLOG_CONSERVATION_ROOT="${ORTHOLOG_CONSERVATION_ROOT:-$SUBMISSION_WORK_ROOT/ORTHOLOG_CONSERVATION}"
export ENSEMBL_ORTHOLOG_ROOT="${ENSEMBL_ORTHOLOG_ROOT:-$STATE_STRUCTURE_ROOT/reference/ENSEMBL_RELEASE115_ORTHOLOGS}"
export PORTABLE_SIGNATURE_AUDIT_ROOT="${PORTABLE_SIGNATURE_AUDIT_ROOT:-$STATE_STRUCTURE_ROOT/tables/PORTABLE_SIGNATURE_AUDIT}"
export ENSEMBL_SIGNATURE_AUDIT_ROOT="${ENSEMBL_SIGNATURE_AUDIT_ROOT:-$STATE_STRUCTURE_ROOT/tables/ENSEMBL_SIGNATURE_AUDIT}"
export RANDOM_SIGNATURE_NULL_ROOT="${RANDOM_SIGNATURE_NULL_ROOT:-$STATE_STRUCTURE_ROOT/tables/RANDOM_SIGNATURE_NULL}"
export SCFEA_OUTPUT_ROOT="${SCFEA_OUTPUT_ROOT:-$SUBMISSION_WORK_ROOT/SCFEA_LEYDIG_CMB}"
export SCFEA_REPO_ROOT="${SCFEA_REPO_ROOT:-$MANUSCRIPT_INPUT_ROOT/external_tools/scFEA}"
export SCFEA_MAX_CELLS_PER_SAMPLE="${SCFEA_MAX_CELLS_PER_SAMPLE:-400}"
export SCFEA_EPOCHS="${SCFEA_EPOCHS:-100}"
export MPLCONFIGDIR="$RUN_ROOT/matplotlib"

LOCKED_ROOT="$LOCKED_INPUT_ROOT"
LOG_ROOT="$RUN_ROOT/logs"
STATUS_FILE="$RUN_ROOT/run_status.tsv"
PREVIOUS_STATUS_FILE="$RUN_ROOT/run_status.previous.tsv"

CONDA_EXE="${CONDA_EXE:-conda}"
PYTHON_BIN="${PYTHON_BIN:-python}"
SCFEA_PYTHON_BIN="${SCFEA_PYTHON_BIN:-$PYTHON_BIN}"
export SCFEA_PYTHON_BIN
BIO_R=("$CONDA_EXE" run --no-capture-output -n bio_r Rscript)
HD_R=("$CONDA_EXE" run --no-capture-output -n hdWGCNA Rscript)

mkdir -p "$HOME" "$LOCKED_ROOT" "$LOG_ROOT" "$MPLCONFIGDIR"

link_input() {
  local source_path="$1"
  local target_path="$2"
  if [[ ! -e "$source_path" ]]; then
    echo "Missing required input: $source_path" >&2
    exit 2
  fi
  mkdir -p "$(dirname "$target_path")"
  if [[ -e "$target_path" && ! -L "$target_path" ]]; then
    echo "Refusing to replace non-symlink input target: $target_path" >&2
    exit 2
  fi
  ln -sfn "$source_path" "$target_path"
}

link_future_output() {
  local source_path="$1"
  local target_path="$2"
  mkdir -p "$(dirname "$target_path")"
  if [[ -e "$target_path" && ! -L "$target_path" ]]; then
    echo "Refusing to replace non-symlink output target: $target_path" >&2
    exit 2
  fi
  ln -sfn "$source_path" "$target_path"
}

setup_inputs() {
  mkdir -p "$LOCKED_ROOT/rds" "$LOCKED_ROOT/tables" "$GSE270931_SCRIPT_ROOT"
  link_input "$MANUSCRIPT_INPUT_ROOT/rds/primary_locked" "$LOCKED_ROOT/rds/primary_locked"
  link_input "$MANUSCRIPT_INPUT_ROOT/tables/primary_locked" "$LOCKED_ROOT/tables/primary_locked"
  link_input "$MANUSCRIPT_INPUT_ROOT/tables/source_copies" "$LOCKED_ROOT/tables/source_copies"
  link_input "$MANUSCRIPT_INPUT_ROOT/raw_data" "$LOCKED_ROOT/raw_data"
  link_input \
    "$SCRIPT_ROOT/2_mouse/1_seurat.R" \
    "$GSE270931_SCRIPT_ROOT/1_seurat.R"

  mkdir -p "$BULK_GSE287203_ROOT/data"
  link_input \
    "$MANUSCRIPT_INPUT_ROOT/raw_data/GSE287203_bulk/GSE287203_raw_counts.txt.gz" \
    "$BULK_GSE287203_ROOT/data/GSE287203_raw_counts.txt.gz"
  link_input \
    "$MANUSCRIPT_INPUT_ROOT/raw_data/GSE287203_bulk/metadata_clean.csv" \
    "$BULK_GSE287203_ROOT/data/metadata_clean.csv"

  local g303_raw="$GSE303193_EXTERNAL_MOUSE_ROOT/raw_data"
  mkdir -p "$g303_raw"
  link_input \
    "$GSE303193_RAW_TAR" \
    "$g303_raw/GSE303193_RAW.tar"

  local g254_raw="$GSE254315_HUMAN_TARGETED_ROOT/raw_data"
  mkdir -p "$g254_raw"
  link_input \
    "$GSE254315_RAW_TAR" \
    "$g254_raw/GSE254315_RAW.tar"
  while IFS= read -r -d '' input_file; do
    link_input "$input_file" "$g254_raw/$(basename "$input_file")"
  done < <(find "$MANUSCRIPT_INPUT_ROOT/raw_data/GSE254315_human_targeted" -maxdepth 1 -type f -print0)

  local g182_raw="$GSE182786_HUMAN_AGING_ROOT/raw_data"
  mkdir -p "$g182_raw"
  while IFS= read -r -d '' input_file; do
    link_input "$input_file" "$g182_raw/$(basename "$input_file")"
  done < <(find "$MANUSCRIPT_INPUT_ROOT/raw_data/GSE182786_human_main" -maxdepth 1 -type f -print0)

  local homologene_ref="$ORTHOLOG_CONSERVATION_ROOT/reference"
  mkdir -p "$homologene_ref"
  link_input \
    "$MANUSCRIPT_INPUT_ROOT/robustness_regulatory_context_20260529/03_orthologs/reference/homologene_build68.data" \
    "$homologene_ref/homologene_build68.data"
}

if [[ -f "$STATUS_FILE" ]]; then
  cp "$STATUS_FILE" "$PREVIOUS_STATUS_FILE"
fi
printf "step_id\tdescription\tstatus\tstarted_utc\tfinished_utc\texit_code\tlog\n" > "$STATUS_FILE"
started=0

run_step() {
  local step_id="$1"
  local description="$2"
  shift 2

  if [[ -n "$START_AT" && "$started" -eq 0 ]]; then
    if [[ "$step_id" == "$START_AT" ]]; then
      started=1
    else
      local retained_log="$LOG_ROOT/${step_id}.log"
      if [[ ! -s "$retained_log" ]]; then
        echo "Cannot retain step $step_id because its prior log is missing or empty: $retained_log" >&2
        exit 2
      fi
      printf "%s\t%s\tpassed_retained\t\t\t0\t%s\n" \
        "$step_id" "$description" "$retained_log" >> "$STATUS_FILE"
      return 0
    fi
  fi

  local log_file="$LOG_ROOT/${step_id}.log"
  local started_utc finished_utc exit_code
  started_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "[$started_utc] START $step_id $description"

  set +e
  "$@" > "$log_file" 2>&1
  exit_code=$?
  set -e

  finished_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [[ "$exit_code" -eq 0 ]]; then
    printf "%s\t%s\tpassed\t%s\t%s\t0\t%s\n" \
      "$step_id" "$description" "$started_utc" "$finished_utc" "$log_file" >> "$STATUS_FILE"
    echo "[$finished_utc] PASS  $step_id $description"
  else
    printf "%s\t%s\tfailed\t%s\t%s\t%s\t%s\n" \
      "$step_id" "$description" "$started_utc" "$finished_utc" "$exit_code" "$log_file" >> "$STATUS_FILE"
    echo "[$finished_utc] FAIL  $step_id $description (exit $exit_code)" >&2
    tail -n 80 "$log_file" >&2
    exit "$exit_code"
  fi
}

setup_inputs

run_step 1 "GSE287203 clean DESeq2" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/1_bulk/1_deseq2.R"
run_step 2 "GSE287203 enrichment" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/1_bulk/2_enrichment.R"

run_step 3 "GSE270931 primary scRNA-seq" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/2_mouse/1_seurat.R"
run_step 4 "GSE270931 state redistribution" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/2_mouse/2_states.R"
run_step 5 "GSE270931 Harmony control" \
  "${HD_R[@]}" "$SCRIPT_ROOT/2_mouse/3_harmony.R"
run_step 6 "Cross-cohort continuum topology" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/2_mouse/4_topology.R"

run_step 7 "Harmony Leydig subset for hdWGCNA" \
  "${HD_R[@]}" "$SCRIPT_ROOT/3_hdwgcna/1_subset.R"
run_step 8 "hdWGCNA network construction" \
  "${HD_R[@]}" "$SCRIPT_ROOT/3_hdwgcna/2_network.R"
run_step 9 "hdWGCNA post-network resume" \
  "${HD_R[@]}" "$SCRIPT_ROOT/3_hdwgcna/3_post_network.R"
run_step 10 "hdWGCNA module enrichment" \
  "${HD_R[@]}" "$SCRIPT_ROOT/3_hdwgcna/4_enrichment.R"

run_step 11 "GSE303193 archive inspection" \
  bash "$SCRIPT_ROOT/4_mouse_external/1_gse303193_fetch.sh"
run_step 12 "GSE303193 Seurat reconstruction" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/4_mouse_external/2_gse303193_seurat.R"
run_step 13 "GSE303193 boundary sensitivity" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/4_mouse_external/3_gse303193_boundary.R"
run_step 14 "GSE303193 pseudobulk" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/4_mouse_external/4_pseudobulk.R"

run_step 15 "GSE254315 archive inspection" \
  bash "$SCRIPT_ROOT/5_human_targeted/1_gse254315_fetch.sh"
run_step 16 "GSE254315 metadata inspection" \
  bash "$SCRIPT_ROOT/5_human_targeted/2_metadata.sh"
run_step 17 "GSE254315 pseudobulk" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/5_human_targeted/3_pseudobulk.R"

run_step 18 "GSE182786 archive inspection" \
  bash "$SCRIPT_ROOT/6_human/1_gse182786_fetch.sh"
run_step 19 "GSE182786 Seurat reconstruction" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/2_gse182786_seurat.R"
run_step 20 "GSE182786 boundary audit" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/3_gse182786_boundary.R"
run_step 21 "GSE182786 BMI downsampling" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/4_bmi.R"
run_step 22 "GSE182786 master module summary" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/module_checks/1_modules.R"
run_step 23 "GSE182786 leave-one-out stability" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/module_checks/2_loo.R"
run_step 24 "GSE182786 core versus extended effects" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/module_checks/3_core_ext.R"
run_step 25 "GSE182786 threshold and sparsity" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/module_checks/4_thresholds.R"
run_step 26 "GSE182786 gene-module decoupling" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/6_human/module_checks/5_decoupling.R"

run_step 27 "Donor/sample-balanced downsampling" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/7_robustness/1_downsampling.R"
run_step 28 "Downsampling summary" \
  "$PYTHON_BIN" "$SCRIPT_ROOT/7_robustness/2_downsampling_summary.py"
run_step 29 "Cluster-boundary audit" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/7_robustness/3_cluster_boundary.R"
run_step 30 "Cross-dataset checks" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/7_robustness/4_cross_dataset.R"

run_step 31 "State-composition analysis" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/8_state_composition/1_composition.R"
run_step 32 "State-composition summary" \
  "$PYTHON_BIN" "$SCRIPT_ROOT/8_state_composition/2_summary.py"

run_step 33 "Portable-signature pseudobulk export" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/9_signature/1_pseudobulk.R"
run_step 34 "HomoloGene ortholog conservation" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/10_null/1_orthologs.R"
run_step 35 "Ensembl release ortholog mapping" \
  "$PYTHON_BIN" "$SCRIPT_ROOT/9_signature/2_orthologs.py"
run_step 36 "HomoloGene portable-signature audit" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/9_signature/3_homologene.R"
run_step 37 "Ensembl portable-signature audit" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/9_signature/4_ensembl.R"
run_step 38 "Size-matched random-signature null" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/10_null/2_random_null.R"
run_step 39 "scFEA input" \
  "${BIO_R[@]}" "$SCRIPT_ROOT/11_scfea/1_input.R"
run_step 40 "scFEA flux" \
  bash "$SCRIPT_ROOT/11_scfea/2_flux.sh"
run_step 41 "scFEA summary" \
  "$SCFEA_PYTHON_BIN" "$SCRIPT_ROOT/11_scfea/3_summary.py"

echo "All submission analysis steps passed."
echo "Run root: $RUN_ROOT"
echo "Status:   $STATUS_FILE"
