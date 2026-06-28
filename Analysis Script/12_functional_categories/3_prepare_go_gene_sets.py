#!/usr/bin/env python3
import csv
import gzip
import os
from collections import defaultdict, deque
from pathlib import Path


WORK_DIR = Path(os.environ.get("FUNCTIONAL_CATEGORY_ROOT", Path.home() / "FUNCTIONAL_CATEGORY_ROOT")).expanduser()
EXT_DIR = Path(os.environ.get("FUNCTIONAL_CATEGORY_ANNOTATION_ROOT", WORK_DIR / "external_annotations")).expanduser()
TABLE_DIR = WORK_DIR / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

OBO = EXT_DIR / "go-basic.obo"
GAF_FILES = {
    "human": EXT_DIR / "goa_human.gaf.gz",
    "mouse": EXT_DIR / "mgi.gaf.gz",
}


CATEGORY_ROOTS = {
    "go_fatty_acid_beta_oxidation": [
        "GO:0006635",
    ],
    "go_ketone_body_metabolism": [
        "GO:0046951",
        "GO:0046952",
    ],
    "go_mitochondrial_energy": [
        "GO:0006119",
        "GO:0022904",
        "GO:0006099",
    ],
    "go_cholesterol_handling": [
        "GO:0008203",
        "GO:0006695",
        "GO:0030301",
    ],
    "go_oxidative_stress_redox": [
        "GO:0006979",
        "GO:0072593",
        "GO:0006749",
    ],
    "go_steroidogenic_execution": [
        "GO:0120178",
        "GO:0006700",
        "GO:0006702",
    ],
    "go_steroid_biosynthesis_broad": [
        "GO:0006694",
    ],
}

METABOLIC_SUPPORT_COMPONENTS = [
    "go_fatty_acid_beta_oxidation",
    "go_ketone_body_metabolism",
    "go_mitochondrial_energy",
    "go_cholesterol_handling",
    "go_oxidative_stress_redox",
]


def parse_obo(path):
    terms = {}
    parents = defaultdict(set)
    current = None

    def flush(term):
        if term and not term.get("obsolete"):
            terms[term["id"]] = {
                "name": term.get("name", ""),
                "namespace": term.get("namespace", ""),
            }

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                flush(current)
                current = {}
                continue
            if current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:]
            elif line.startswith("is_obsolete: true"):
                current["obsolete"] = True
            elif line.startswith("is_a: "):
                parent = line.split()[1]
                if "id" in current:
                    parents[current["id"]].add(parent)
            elif line.startswith("relationship: part_of "):
                parent = line.split()[2]
                if "id" in current:
                    parents[current["id"]].add(parent)
        flush(current)

    children = defaultdict(set)
    for child, parent_set in parents.items():
        for parent in parent_set:
            children[parent].add(child)
    return terms, children


def descendants(root, children):
    seen = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        for child in children.get(node, ()):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return seen


def parse_gaf(path):
    go_to_genes = defaultdict(set)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            symbol = parts[2].strip()
            qualifier = parts[3].strip()
            go_id = parts[4].strip()
            aspect = parts[8].strip()
            if not symbol or not go_id:
                continue
            if "NOT" in qualifier.split("|"):
                continue
            if aspect != "P":
                continue
            go_to_genes[go_id].add(symbol)
    return go_to_genes


def main():
    required = [OBO, *GAF_FILES.values()]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing GO annotation input(s): " + ", ".join(missing))

    terms, children = parse_obo(OBO)

    category_terms = {}
    for category, roots in CATEGORY_ROOTS.items():
        included = set()
        for root in roots:
            if root not in terms:
                raise RuntimeError(f"GO root not found: {root}")
            included.update(descendants(root, children))
        category_terms[category] = included

    with (TABLE_DIR / "go_external_category_terms_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["category", "root_go_ids", "included_go_id", "included_go_name", "namespace"])
        for category, included in sorted(category_terms.items()):
            roots = ";".join(CATEGORY_ROOTS[category])
            for go_id in sorted(included):
                info = terms.get(go_id, {})
                writer.writerow([category, roots, go_id, info.get("name", ""), info.get("namespace", "")])

    gene_rows = []
    for species, gaf_path in GAF_FILES.items():
        go_to_genes = parse_gaf(gaf_path)
        category_to_genes = {}
        for category, go_ids in category_terms.items():
            genes = set()
            for go_id in go_ids:
                genes.update(go_to_genes.get(go_id, set()))
            category_to_genes[category] = genes

        union = set()
        for component in METABOLIC_SUPPORT_COMPONENTS:
            union.update(category_to_genes.get(component, set()))
        category_to_genes["go_metabolic_support_union_prespecified"] = union
        category_to_genes["go_metabolic_support_union_no_hmgcs2"] = {
            g for g in union if g.upper() != "HMGCS2"
        }

        for category, genes in sorted(category_to_genes.items()):
            for gene in sorted(genes):
                gene_rows.append([species, category, gene])

    with (TABLE_DIR / "go_external_gene_sets_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", "category", "gene"])
        writer.writerows(gene_rows)

    counts = defaultdict(int)
    for species, category, gene in gene_rows:
        counts[(species, category)] += 1
    with (TABLE_DIR / "go_external_gene_set_sizes_20260619.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", "category", "n_genes"])
        for (species, category), n in sorted(counts.items()):
            writer.writerow([species, category, n])

    print("Wrote GO external gene sets")
    print(TABLE_DIR / "go_external_gene_sets_20260619.csv")


if __name__ == "__main__":
    main()
