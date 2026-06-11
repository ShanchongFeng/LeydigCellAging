#!/usr/bin/env bash
set -euo pipefail

SCFEA_REPO_ROOT="${SCFEA_REPO_ROOT:-}"
SCFEA_OUTPUT_ROOT="${SCFEA_OUTPUT_ROOT:-}"
SCFEA_PYTHON_BIN="${SCFEA_PYTHON_BIN:-python}"
SCFEA_EPOCHS="${SCFEA_EPOCHS:-100}"

if [[ -z "$SCFEA_REPO_ROOT" ]]; then
  echo "Set SCFEA_REPO_ROOT to the cloned scFEA repository root." >&2
  exit 2
fi
if [[ -z "$SCFEA_OUTPUT_ROOT" ]]; then
  echo "Set SCFEA_OUTPUT_ROOT to the isolated scFEA output root." >&2
  exit 2
fi

SCFEA_SCRIPT="$SCFEA_REPO_ROOT/src/scFEA.py"
if [[ ! -f "$SCFEA_SCRIPT" ]]; then
  echo "Missing scFEA.py: $SCFEA_SCRIPT" >&2
  exit 2
fi

mkdir -p "$SCFEA_OUTPUT_ROOT/raw" "$SCFEA_OUTPUT_ROOT/output"
cd "$SCFEA_OUTPUT_ROOT"

PYTHONWARNINGS="${PYTHONWARNINGS:-ignore::FutureWarning}" "$SCFEA_PYTHON_BIN" "$SCFEA_SCRIPT" \
  --data_dir "$SCFEA_REPO_ROOT/data" \
  --input_dir input \
  --res_dir raw \
  --test_file GSE182786_Leydig_CMB_0_17_scFEA_input_lognorm.csv \
  --moduleGene_file module_gene_m168.csv \
  --stoichiometry_matrix cmMat_c70_m168.csv \
  --cName_file cName_c70_m168.csv \
  --sc_imputation False \
  --output_flux_file raw/GSE182786_Leydig_CMB_0_17_flux_epoch100.csv \
  --output_balance_file raw/GSE182786_Leydig_CMB_0_17_balance_epoch100.csv \
  --train_epoch "$SCFEA_EPOCHS"

test -s raw/GSE182786_Leydig_CMB_0_17_flux_epoch100.csv
test -s raw/GSE182786_Leydig_CMB_0_17_balance_epoch100.csv
