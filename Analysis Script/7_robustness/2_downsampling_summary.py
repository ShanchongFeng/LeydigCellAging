import os
from pathlib import Path

import pandas as pd


BASE = Path(
    os.environ.get(
        "DONOR_SAMPLE_DOWNSAMPLING_ROOT",
        str(Path.home() / "DONOR_SAMPLE_DOWNSAMPLING_ROOT"),
    )
)

TERMS = [
    "Hmgcs2",
    "HMGCS2",
    "Cyp11a1",
    "CYP11A1",
    "Insl3",
    "INSL3",
    "Star",
    "STAR",
    "hmgcs2_axis",
    "steroidogenesis",
    "tan_module1",
    "greenyellow_module1",
    "green_module1",
]


def contains_priority(readout: str) -> bool:
    return any(term in readout for term in TERMS)


def main() -> None:
    design = pd.read_csv(BASE / "downsampling_design.csv")
    summary = pd.read_csv(BASE / "downsampling_priority_readouts_summary.csv")
    priority = summary[summary["readout"].map(contains_priority)].copy()

    cols = [
        "dataset",
        "boundary",
        "readout",
        "median_effect",
        "ci95_low",
        "ci95_high",
        "pct_aged_lower",
        "pct_aged_higher",
        "direction_consistency",
        "direction",
        "median_p_value",
    ]
    cols = [c for c in cols if c in priority.columns]

    _ = design
    _ = priority[cols]


if __name__ == "__main__":
    main()
