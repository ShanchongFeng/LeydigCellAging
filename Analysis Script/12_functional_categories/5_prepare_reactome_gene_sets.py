#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
from pathlib import Path


WORK_DIR = Path(os.environ.get("FUNCTIONAL_CATEGORY_ROOT", Path.home() / "FUNCTIONAL_CATEGORY_ROOT")).expanduser()
EXT_DIR = Path(os.environ.get("FUNCTIONAL_CATEGORY_ANNOTATION_ROOT", WORK_DIR / "external_annotations")).expanduser()
TABLE_DIR = WORK_DIR / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

REACTOME = EXT_DIR / "NCBI2Reactome_All_Levels.txt"

NUMERIC_IDS = {
    "reactome_fatty_acid_beta_oxidation": [
        "389887",
        "390247",
        "77286",
        "77288",
        "77289",
    ],
    "reactome_ketone_body_metabolism": [
        "74182",
        "77108",
        "77111",
    ],
    "reactome_mitochondrial_energy": [
        "1428517",
        "611105",
        "71403",
    ],
    "reactome_cholesterol_handling": [
        "191273",
        "1655829",
        "6807047",
        "9029569",
        "9031525",
        "9969901",
    ],
    "reactome_oxidative_stress_redox": [
        "3299685",
        "156590",
        "174403",
    ],
    "reactome_steroidogenic_execution": [
        "193048",
        "196071",
    ],
    "reactome_steroid_metabolism_broad": [
        "8957322",
    ],
}

METABOLIC_COMPONENTS = [
    "reactome_fatty_acid_beta_oxidation",
    "reactome_ketone_body_metabolism",
    "reactome_mitochondrial_energy",
    "reactome_cholesterol_handling",
    "reactome_oxidative_stress_redox",
]

SPECIES_PREFIX = {
    "Homo sapiens": ("human", "R-HSA-"),
    "Mus musculus": ("mouse", "R-MMU-"),
}


def main():
    if not REACTOME.is_file():
        raise FileNotFoundError(f"Missing Reactome annotation input: {REACTOME}")

    category_pathway_ids = defaultdict(set)
    for category, ids in NUMERIC_IDS.items():
        for species_name, (_, prefix) in SPECIES_PREFIX.items():
            for numeric in ids:
                category_pathway_ids[category].add(prefix + numeric)

    pathway_meta = {}
    sets = defaultdict(lambda: defaultdict(set))
    rows_by_pathway = []

    with REACTOME.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            entrez, pathway_id, url, pathway_name, evidence, species_name = parts[:6]
            if species_name not in SPECIES_PREFIX:
                continue
            species, _ = SPECIES_PREFIX[species_name]
            pathway_meta[pathway_id] = (pathway_name, species_name, url)
            for category, pathway_ids in category_pathway_ids.items():
                if pathway_id in pathway_ids:
                    sets[(species, category)][entrez].add((pathway_id, pathway_name))

    metabolic_union = defaultdict(set)
    for species in ("human", "mouse"):
        for component in METABOLIC_COMPONENTS:
            metabolic_union[species].update(sets[(species, component)].keys())
        for entrez in metabolic_union[species]:
            sets[(species, "reactome_metabolic_support_union_prespecified")][entrez].add(
                ("UNION", "Prespecified Reactome metabolic-support union")
            )

    with (TABLE_DIR / "reactome_external_category_pathways_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["category", "reactome_id", "pathway_name", "species", "url"])
        for category, pathway_ids in sorted(category_pathway_ids.items()):
            for pid in sorted(pathway_ids):
                name, species_name, url = pathway_meta.get(pid, ("NOT_FOUND", "", ""))
                writer.writerow([category, pid, name, species_name, url])

    gene_rows = []
    for (species, category), entrez_map in sorted(sets.items()):
        for entrez, memberships in sorted(entrez_map.items(), key=lambda x: int(x[0]) if x[0].isdigit() else x[0]):
            pids = sorted({m[0] for m in memberships})
            names = sorted({m[1] for m in memberships})
            gene_rows.append([species, category, entrez, ";".join(pids), ";".join(names)])

    with (TABLE_DIR / "reactome_external_gene_sets_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", "category", "entrez_id", "reactome_ids", "reactome_names"])
        writer.writerows(gene_rows)

    counts = defaultdict(int)
    for species, category, entrez, _, _ in gene_rows:
        counts[(species, category)] += 1
    with (TABLE_DIR / "reactome_external_gene_set_sizes_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", "category", "n_entrez"])
        for (species, category), n in sorted(counts.items()):
            writer.writerow([species, category, n])

    print("Wrote Reactome external gene sets")
    print(TABLE_DIR / "reactome_external_gene_sets_20260619.csv")


if __name__ == "__main__":
    main()
