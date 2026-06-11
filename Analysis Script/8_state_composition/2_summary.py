#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(
    os.environ.get(
        "STATE_STRUCTURE_ROOT",
        str(Path.home() / "STATE_STRUCTURE_ROOT"),
    )
)
TABLE_DIR = ROOT / "tables" / "1_state_composition"
DOC_DIR = ROOT / "docs"


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[idx] * n / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


effects = read_csv(TABLE_DIR / "exact_permutation_bootstrap_effects.csv")
groups = defaultdict(list)
for row in effects:
    key = (
        row["dataset"],
        row["boundary"],
        row["denominator"],
        row["comparison"],
    )
    groups[key].append(row)

for rows in groups.values():
    qvals = bh_adjust([float(row["permutation_p"]) for row in rows])
    for row, qval in zip(rows, qvals):
        row["permutation_fdr_within_test_family"] = f"{qval:.10g}"

fieldnames = list(effects[0].keys())
write_csv(TABLE_DIR / "composition_effects_with_fdr.csv", effects, fieldnames)
