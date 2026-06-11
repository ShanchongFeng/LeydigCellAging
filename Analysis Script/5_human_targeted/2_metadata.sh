#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${GSE254315_HUMAN_TARGETED_ROOT:-$HOME/GSE254315_HUMAN_TARGETED_ROOT}"
RAW_DIR="$BASE_DIR/raw_data"
TABLE_DIR="$BASE_DIR/tables"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$RAW_DIR" "$TABLE_DIR" "$LOG_DIR"

cd "$RAW_DIR"

HTML="GSE254315_series.html"
SOFT="GSE254315_family.soft.gz"
MATRIX="GSE254315_series_matrix.txt.gz"

if [ ! -f "$HTML" ]; then
  wget -q -O "$HTML" "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254315"
fi

if [ ! -f "$SOFT" ]; then
  wget -q -O "$SOFT" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/soft/GSE254315_family.soft.gz"
fi

if [ ! -f "$MATRIX" ]; then
  wget -q -O "$MATRIX" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/matrix/GSE254315_series_matrix.txt.gz"
fi

zgrep -c '^\^SAMPLE' "$SOFT" > "$TABLE_DIR/4_soft_sample_count.txt"

zgrep -n '^\^SAMPLE\|!Sample_title' "$SOFT" > "$TABLE_DIR/5_soft_sample_titles.txt"

zgrep -n '^\^SAMPLE\|!Sample_title\|!Sample_geo_accession\|!Sample_source_name_ch1\|!Sample_characteristics_ch1\|!Sample_description\|!Sample_supplementary_file' "$SOFT" \
  > "$TABLE_DIR/6_soft_sample_detailed_metadata.txt"

zgrep -niE 'age|aged|young|old|donor|patient|sample|testis|testicular|adult|puberty|Leydig|cell type|celltype|organ|tissue|sex|male' "$SOFT" \
  > "$TABLE_DIR/7_soft_age_donor_related_lines.txt" || true

zgrep -niE 'age|aged|young|old|donor|patient|sample|testis|testicular|adult|puberty|Leydig|cell type|celltype|organ|tissue|sex|male' "$MATRIX" \
  > "$TABLE_DIR/8_series_matrix_age_donor_related_lines.txt" || true

