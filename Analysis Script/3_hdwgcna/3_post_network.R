#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

base_dir <- path.expand(Sys.getenv(
  "HDWGCNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT")
))
table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
rds_dir <- file.path(base_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(rds_dir, "16_after_ConstructNetwork_softpower7.rds")


if (!file.exists(input_rds)) {
  stop("Missing input RDS: ", input_rds)
}

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"



obj <- ModuleEigengenes(
  obj,
  group.by.vars = NULL,
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "17_after_ModuleEigengenes_raw_fallback.rds"))


obj <- ModuleConnectivity(
  obj,
  group.by = "network_group",
  group_name = "all",
  harmonized = FALSE,
  assay = "RNA",
  layer = "data",
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "18_after_ModuleConnectivity_raw_fallback.rds"))


mods <- GetModules(obj, wgcna_name = "leydig_hd")

write.csv(
  mods,
  file.path(table_dir, "hdwgcna_module_assignment.csv"),
  row.names = FALSE
)


module_size <- as.data.frame(table(mods$module))
colnames(module_size) <- c("module", "n_genes")
module_size <- module_size %>% arrange(desc(n_genes))

write.csv(
  module_size,
  file.path(table_dir, "hdwgcna_module_size_summary.csv"),
  row.names = FALSE
)



hub_df <- GetHubGenes(
  obj,
  n_hubs = 30,
  wgcna_name = "leydig_hd"
)

write.csv(
  hub_df,
  file.path(table_dir, "hdwgcna_hub_genes.csv"),
  row.names = FALSE
)


target_genes <- c(
  "Hmgcs2", "Cyp11a1", "Lhcgr", "Foxo3", "Nr5a1",
  "Star", "Hsd3b1", "Cyp17a1", "Insl3"
)

target_mods <- mods %>%
  filter(gene_name %in% target_genes)

write.csv(
  target_mods,
  file.path(table_dir, "hdwgcna_target_gene_modules.csv"),
  row.names = FALSE
)


mods_keep <- intersect(c("tan", "greenyellow", "pink", "green"), unique(mods$module))


for (m in mods_keep) {
  df_m <- mods %>% filter(module == m)
  write.csv(
    df_m,
    file.path(table_dir, paste0("hdwgcna_module_", m, "_genes.csv")),
    row.names = FALSE
  )
}

saveRDS(obj, file.path(rds_dir, "19_hdwgcna_final_object_raw_fallback.rds"))

