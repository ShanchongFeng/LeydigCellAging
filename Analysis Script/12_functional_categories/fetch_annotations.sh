#!/usr/bin/env bash
set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${FUNCTIONAL_CATEGORY_ANNOTATION_ROOT:-$SCRIPT_ROOT/external_annotations}"
MANIFEST="$SCRIPT_ROOT/ANNOTATION_INPUTS.tsv"
mkdir -p "$OUT_DIR"

hash_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

while IFS=$'\t' read -r name url expected frozen_on; do
  [[ "$name" == "file" ]] && continue
  target="$OUT_DIR/$name"
  if [[ ! -f "$target" ]]; then
    curl --fail --location --retry 3 --output "$target" "$url"
  fi
  observed="$(hash_file "$target")"
  if [[ "$observed" != "$expected" ]]; then
    echo "Checksum mismatch for $name (frozen $frozen_on)." >&2
    echo "Expected: $expected" >&2
    echo "Observed: $observed" >&2
    echo "Use the archived frozen input; current upstream annotations may have changed." >&2
    exit 1
  fi
done < "$MANIFEST"

echo "Functional-category annotation inputs verified in $OUT_DIR"
