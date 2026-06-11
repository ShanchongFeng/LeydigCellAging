#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${GSE254315_HUMAN_TARGETED_ROOT:-$HOME/GSE254315_HUMAN_TARGETED_ROOT}"
RAW_DIR="$BASE_DIR/raw_data"
LOG_DIR="$BASE_DIR/logs"
TABLE_DIR="$BASE_DIR/tables"

mkdir -p "$RAW_DIR" "$LOG_DIR" "$TABLE_DIR"

cd "$RAW_DIR"

RAW_TAR="$RAW_DIR/GSE254315_RAW.tar"
RAW_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE254315&format=file"

if [ ! -f "$RAW_TAR" ]; then
  wget -q -O "$RAW_TAR" "$RAW_URL"
fi


tar -tf "$RAW_TAR" > "$TABLE_DIR/0_raw_tar_all_entries.txt"

head -n 300 "$TABLE_DIR/0_raw_tar_all_entries.txt" > "$TABLE_DIR/0_raw_tar_first300.txt"

tar -xf "$RAW_TAR" -C "$RAW_DIR"

find "$RAW_DIR" -maxdepth 4 -type f | sed "s|$RAW_DIR/||" | sort > "$TABLE_DIR/1_extracted_files_maxdepth4.txt"

find "$RAW_DIR" -maxdepth 6 -type f | grep -Ei "matrix|barcodes|features|genes|mtx|h5|csv|tsv|txt|rds|h5ad|metadata|meta|annot|anno|cell" | sed "s|$RAW_DIR/||" | sort > "$TABLE_DIR/2_candidate_matrix_metadata_files.txt"

grep -Ei "donor|sample|age|young|old|aged|adult|testis|leydig|matrix|barcodes|features|mtx|h5|metadata|cell" "$TABLE_DIR/0_raw_tar_all_entries.txt" > "$TABLE_DIR/3_sample_donor_age_related_entries.txt" || true

