#!/usr/bin/env Rscript

STATE_ROOT <- path.expand(Sys.getenv(
  "STATE_STRUCTURE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "STATE_STRUCTURE_ROOT")
))

Sys.setenv(
  WP4_ORTHOLOG_FILE = file.path(
    path.expand(Sys.getenv("ENSEMBL_ORTHOLOG_ROOT", unset = file.path(Sys.getenv("HOME"), "ENSEMBL_ORTHOLOG_ROOT"))),
    "ensembl_tan_mouse_human_one_to_one_map.csv"
  ),
  WP4_TABLE_DIR = path.expand(Sys.getenv("ENSEMBL_SIGNATURE_AUDIT_ROOT", unset = file.path(Sys.getenv("HOME"), "ENSEMBL_SIGNATURE_AUDIT_ROOT"))),
  WP4_FIG_PREFIX = "WP4_ensembl_release115_signature_gene_z_effect_forest",
  WP4_DEFINITION_FILE = "WP4_ensembl_release115_analysis_definition.md",
  WP4_SESSION_FILE = "4_ensembl_sessionInfo.txt",
  WP4_MAPPING_LABEL = "Ensembl REST homology release 115 ortholog_one2one"
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Could not resolve wrapper script path")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "3_homologene.R"))
