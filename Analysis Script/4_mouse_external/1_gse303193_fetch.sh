#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${GSE303193_EXTERNAL_MOUSE_ROOT:-$HOME/GSE303193_EXTERNAL_MOUSE_ROOT}"
RAW_DIR="$BASE_DIR/raw_data"
LOG_DIR="$BASE_DIR/logs"
TABLE_DIR="$BASE_DIR/tables"

mkdir -p "$RAW_DIR" "$LOG_DIR" "$TABLE_DIR"

cd "$RAW_DIR"

RAW_TAR="$RAW_DIR/GSE303193_RAW.tar"
RAW_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE303193&format=file"

if [ ! -f "$RAW_TAR" ]; then
  wget -q -O "$RAW_TAR" "$RAW_URL"
fi


tar -tf "$RAW_TAR" > "$TABLE_DIR/0_raw_tar_all_entries.txt"

head -n 200 "$TABLE_DIR/0_raw_tar_all_entries.txt" > "$TABLE_DIR/0_raw_tar_first200.txt"

tar -xf "$RAW_TAR" -C "$RAW_DIR"

find "$RAW_DIR" -maxdepth 3 -type f | sed "s|$RAW_DIR/||" | sort > "$TABLE_DIR/1_extracted_files_maxdepth3.txt"

find "$RAW_DIR" -maxdepth 3 -type f -exec ls -lh {} \; | sort -k 9 > "$TABLE_DIR/2_extracted_file_sizes.txt"

find "$RAW_DIR" -maxdepth 5 -type f | grep -Ei "matrix|barcodes|features|genes|mtx|h5|csv|tsv|txt|gz" | sed "s|$RAW_DIR/||" | sort > "$TABLE_DIR/3_candidate_matrix_files.txt"

grep -Ei "M5|M20|young|aged|old|sample|matrix|barcodes|features|mtx|h5" "$TABLE_DIR/0_raw_tar_all_entries.txt" > "$TABLE_DIR/4_sample_name_related_entries.txt" || true

