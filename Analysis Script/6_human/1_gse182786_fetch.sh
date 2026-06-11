#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${GSE182786_HUMAN_AGING_ROOT:-$HOME/GSE182786_HUMAN_AGING_ROOT}"
RAW_DIR="$BASE_DIR/raw_data"
LOG_DIR="$BASE_DIR/logs"
TABLE_DIR="$BASE_DIR/tables"

mkdir -p "$RAW_DIR" "$LOG_DIR" "$TABLE_DIR"

cd "$RAW_DIR"
if [ ! -f "GSE182786_series.html" ]; then
  wget -q -O GSE182786_series.html "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182786"
fi

if [ ! -f "GSE182786_family.soft.gz" ]; then
  wget -q -O GSE182786_family.soft.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182786/soft/GSE182786_family.soft.gz"
fi

if [ ! -f "GSE182786_series_matrix.txt.gz" ]; then
  wget -q -O GSE182786_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182786/matrix/GSE182786_series_matrix.txt.gz"
fi
if [ ! -f "GSE182786_Young_matrix.txt.gz" ]; then
  wget -q -O GSE182786_Young_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182786/suppl/GSE182786_Young_matrix.txt.gz"
fi

if [ ! -f "GSE182786_Older_Group1_matrix.txt.gz" ]; then
  wget -q -O GSE182786_Older_Group1_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182786/suppl/GSE182786_Older_Group1_matrix.txt.gz"
fi

if [ ! -f "GSE182786_Older_Group2_matrix.txt.gz" ]; then
  wget -q -O GSE182786_Older_Group2_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182786/suppl/GSE182786_Older_Group2_matrix.txt.gz"
fi

zgrep -c '^\^SAMPLE' GSE182786_family.soft.gz > "$TABLE_DIR/0_soft_sample_count.txt"

zgrep -n '^\^SAMPLE\|!Sample_title' GSE182786_family.soft.gz > "$TABLE_DIR/1_soft_sample_titles.txt"

zgrep -n '^\^SAMPLE\|!Sample_title\|!Sample_source_name_ch1\|!Sample_characteristics_ch1\|!Sample_supplementary_file' GSE182786_family.soft.gz \
  > "$TABLE_DIR/2_soft_sample_detailed_metadata.txt"

zgrep -niE 'age|aged|old|older|young|group|donor|sample|testis|testicular|Leydig|cell|matrix|metadata' GSE182786_family.soft.gz \
  > "$TABLE_DIR/3_soft_age_group_related_lines.txt" || true

zgrep -niE 'age|aged|old|older|young|group|donor|sample|testis|testicular|Leydig|cell|matrix|metadata' GSE182786_series_matrix.txt.gz \
  > "$TABLE_DIR/4_series_matrix_age_group_related_lines.txt" || true

gzip -t GSE182786_Young_matrix.txt.gz GSE182786_Older_Group1_matrix.txt.gz GSE182786_Older_Group2_matrix.txt.gz
set +o pipefail
for f in GSE182786_Young_matrix.txt.gz GSE182786_Older_Group1_matrix.txt.gz GSE182786_Older_Group2_matrix.txt.gz; do
  zcat "$f" | head -n 8

  zcat "$f" | head -n 3 | sed 's/\t/[TAB]/g'

  zcat "$f" | wc -l
  zcat "$f" | head -n 1 | awk -F'\t' '{print NF}'
done > "$TABLE_DIR/5_matrix_preview.txt"
set -o pipefail

md5sum GSE182786_*matrix.txt.gz > "$TABLE_DIR/6_matrix_md5sum.txt"

