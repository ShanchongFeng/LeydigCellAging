#!/usr/bin/env python3
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

LOCKED_ROOT = Path(
    os.environ.get(
        "LOCKED_INPUT_ROOT",
        str(Path.home() / "LOCKED_INPUT_ROOT"),
    )
)
OUT = Path(
    os.environ.get(
        "ENSEMBL_ORTHOLOG_ROOT",
        str(Path.home() / "ENSEMBL_ORTHOLOG_ROOT"),
    )
)
TAN_FILE = (
    LOCKED_ROOT
    / "tables"
    / "primary_locked"
    / "hdWGCNA_enrichment"
    / "leydig_scRNA_phase2a__phase2b__hdwgcna__enrich__csv__module_tan_genes.csv"
)
REST = "https://rest.ensembl.org"
HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/json",
    "User-Agent": "leydig-aging-ortholog-audit/1.0",
}


def get_json(url, attempts=4, retry_http_400=False):
    last = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(request, timeout=45) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return {"_http_error": error.code}
            if error.code == 400 and not retry_http_400:
                return {"_http_error": error.code}
            last = error
        except Exception as error:
            last = error
        time.sleep(1.0 + attempt)
    if isinstance(last, urllib.error.HTTPError):
        return {"_http_error": last.code}
    raise RuntimeError(f"Request failed after {attempts} attempts: {url}: {last!r}")


def read_tan_genes():
    with TAN_FILE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    genes = []
    for row in rows:
        gene = row.get("gene_name") or row.get("gene") or next(iter(row.values()))
        if gene and gene not in genes:
            genes.append(gene)
    return genes


def resolve_human_symbol(ensembl_gene_id):
    payload = get_json(f"{REST}/lookup/id/{urllib.parse.quote(ensembl_gene_id)}")
    if "_http_error" in payload:
        return ""
    return payload.get("display_name", "")


def lookup_gene(ensembl_gene_id):
    payload = get_json(f"{REST}/lookup/id/{urllib.parse.quote(ensembl_gene_id)}")
    if "_http_error" in payload:
        return {}
    return payload


def resolve_mouse_gene(mouse_symbol):
    # Prefer the unambiguous current symbol lookup. The homology-by-symbol
    # endpoint can resolve short aliases such as C2 or Ostc to unrelated genes.
    symbol_url = (
        f"{REST}/lookup/symbol/mus_musculus/{urllib.parse.quote(mouse_symbol)}"
        "?expand=0"
    )
    direct = get_json(symbol_url)
    if "_http_error" not in direct and direct.get("id"):
        return direct.get("id", ""), direct.get("display_name", ""), "lookup_symbol"

    # For historical aliases such as Aes, use MGI xrefs and inspect the
    # resolved source records. A unique xref is accepted; otherwise prefer an
    # exact display-name match and fail closed if ambiguity remains.
    xref_url = (
        f"{REST}/xrefs/symbol/mus_musculus/{urllib.parse.quote(mouse_symbol)}"
        "?external_db=MGI"
    )
    xrefs = get_json(xref_url)
    if isinstance(xrefs, dict) and "_http_error" in xrefs:
        return "", "", "no_source_match"
    candidates = []
    for xref in xrefs:
        gene_id = xref.get("id", "")
        if not gene_id:
            continue
        info = lookup_gene(gene_id)
        if info.get("object_type") == "Gene":
            candidates.append((gene_id, info.get("display_name", "")))
    exact = [
        candidate for candidate in candidates
        if candidate[1].casefold() == mouse_symbol.casefold()
    ]
    if len(exact) == 1:
        return exact[0][0], exact[0][1], "mgi_xref_exact_display_name"
    if len(candidates) == 1:
        return candidates[0][0], candidates[0][1], "mgi_xref_unique_alias"
    return "", "", "ambiguous_or_missing_source"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    run_time = datetime.now(timezone.utc).isoformat()
    release_payload = get_json(f"{REST}/info/data")
    releases = release_payload.get("releases", [])
    release = releases[0] if releases else "unknown"
    tan_genes = read_tan_genes()

    all_rows = []
    for index, mouse_symbol in enumerate(tan_genes, start=1):
        source_id, current_mouse_symbol, source_resolution = resolve_mouse_gene(mouse_symbol)
        if not source_id:
            all_rows.append(
                {
                    "module": "hdwgcna_tan",
                    "mouse_symbol": mouse_symbol,
                    "mouse_current_symbol": current_mouse_symbol,
                    "mouse_ensembl_gene_id": "",
                    "human_symbol": "",
                    "human_ensembl_gene_id": "",
                    "orthology_type": "",
                    "method_link_type": "",
                    "taxonomy_level": "",
                    "source_resolution": source_resolution,
                    "query_status": "no_unambiguous_source_gene",
                    "ensembl_release": release,
                    "query_timestamp_utc": run_time,
                }
            )
            continue
        url = (
            f"{REST}/homology/id/mus_musculus/{urllib.parse.quote(source_id)}"
            "?target_species=homo_sapiens;type=orthologues;format=condensed"
        )
        # The Ensembl homology endpoint can occasionally return a transient
        # HTTP 400 for valid stable IDs, so retry before recording a failure.
        payload = get_json(url, retry_http_400=True)
        if "_http_error" in payload:
            all_rows.append(
                {
                    "module": "hdwgcna_tan",
                    "mouse_symbol": mouse_symbol,
                    "mouse_current_symbol": current_mouse_symbol,
                    "mouse_ensembl_gene_id": source_id,
                    "human_symbol": "",
                    "human_ensembl_gene_id": "",
                    "orthology_type": "",
                    "method_link_type": "",
                    "taxonomy_level": "",
                    "source_resolution": source_resolution,
                    "query_status": f"http_{payload['_http_error']}",
                    "ensembl_release": release,
                    "query_timestamp_utc": run_time,
                }
            )
            continue
        data = payload.get("data", [])
        if not data:
            all_rows.append(
                {
                    "module": "hdwgcna_tan",
                    "mouse_symbol": mouse_symbol,
                    "mouse_current_symbol": current_mouse_symbol,
                    "mouse_ensembl_gene_id": source_id,
                    "human_symbol": "",
                    "human_ensembl_gene_id": "",
                    "orthology_type": "",
                    "method_link_type": "",
                    "taxonomy_level": "",
                    "source_resolution": source_resolution,
                    "query_status": "no_source_match",
                    "ensembl_release": release,
                    "query_timestamp_utc": run_time,
                }
            )
            continue
        source_id = data[0].get("id", source_id)
        homologies = data[0].get("homologies", [])
        if not homologies:
            all_rows.append(
                {
                    "module": "hdwgcna_tan",
                    "mouse_symbol": mouse_symbol,
                    "mouse_current_symbol": current_mouse_symbol,
                    "mouse_ensembl_gene_id": source_id,
                    "human_symbol": "",
                    "human_ensembl_gene_id": "",
                    "orthology_type": "",
                    "method_link_type": "",
                    "taxonomy_level": "",
                    "source_resolution": source_resolution,
                    "query_status": "no_human_ortholog",
                    "ensembl_release": release,
                    "query_timestamp_utc": run_time,
                }
            )
            continue
        for homology in homologies:
            human_id = homology.get("id", "")
            all_rows.append(
                {
                    "module": "hdwgcna_tan",
                    "mouse_symbol": mouse_symbol,
                    "mouse_current_symbol": current_mouse_symbol,
                    "mouse_ensembl_gene_id": source_id,
                    "human_symbol": resolve_human_symbol(human_id) if human_id else "",
                    "human_ensembl_gene_id": human_id,
                    "orthology_type": homology.get("type", ""),
                    "method_link_type": homology.get("method_link_type", ""),
                    "taxonomy_level": homology.get("taxonomy_level", ""),
                    "source_resolution": source_resolution,
                    "query_status": "mapped",
                    "ensembl_release": release,
                    "query_timestamp_utc": run_time,
                }
            )
            time.sleep(0.08)
        time.sleep(0.08)

    all_file = OUT / "ensembl_tan_mouse_human_homology_all.csv"
    fieldnames = list(all_rows[0])
    with all_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    per_mouse_counts = {}
    for row in all_rows:
        if row["orthology_type"] == "ortholog_one2one" and row["human_symbol"]:
            per_mouse_counts.setdefault(row["mouse_symbol"], []).append(row)

    one_to_one_rows = []
    for mouse_symbol in tan_genes:
        hits = per_mouse_counts.get(mouse_symbol, [])
        if len(hits) == 1:
            row = dict(hits[0])
            row["mapped_one_to_one"] = True
        else:
            row = {
                "module": "hdwgcna_tan",
                "mouse_symbol": mouse_symbol,
                "mouse_current_symbol": "",
                "mouse_ensembl_gene_id": "",
                "human_symbol": "",
                "human_ensembl_gene_id": "",
                "orthology_type": "",
                "method_link_type": "",
                "taxonomy_level": "",
                "source_resolution": "",
                "query_status": "no_unique_one_to_one",
                "ensembl_release": release,
                "query_timestamp_utc": run_time,
                "mapped_one_to_one": False,
            }
        one_to_one_rows.append(row)

    map_file = OUT / "ensembl_tan_mouse_human_one_to_one_map.csv"
    map_fields = list(one_to_one_rows[0])
    with map_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=map_fields)
        writer.writeheader()
        writer.writerows(one_to_one_rows)

    manifest = {
        "query_timestamp_utc": run_time,
        "ensembl_release": release,
        "rest_base_url": REST,
        "homology_endpoint": "/homology/id/mus_musculus/{ensembl_gene_id}",
        "lookup_endpoint": "/lookup/id/{ensembl_gene_id}",
        "target_species": "homo_sapiens",
        "orthology_filter": "ortholog_one2one",
        "n_tan_mouse_genes": len(tan_genes),
        "n_unique_one_to_one": sum(bool(row["mapped_one_to_one"]) for row in one_to_one_rows),
    }
    (OUT / "ensembl_mapping_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
