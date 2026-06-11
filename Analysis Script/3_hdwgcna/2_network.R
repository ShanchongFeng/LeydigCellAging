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
input_rds <- path.expand(Sys.getenv(
  "HDWGCNA_PRENETWORK_INPUT_RDS",
  unset = file.path(base_dir, "rds", "2_leydig_harmony_subset.rds")
))

table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
rds_dir <- file.path(base_dir, "rds")
tmp_dir <- file.path(base_dir, "tmp")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)


if (!file.exists(input_rds)) {
  stop("Input Harmony Leydig subset not found: ", input_rds)
}

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "group"


# Required metadata
if (!("group" %in% colnames(obj@meta.data))) {
  stop("Missing metadata column: group")
}

if (!("harmony" %in% names(obj@reductions))) {
  stop("Harmony reduction not found in object. Existing reductions: ",
       paste(names(obj@reductions), collapse = ", "))
}

saveRDS(obj, file.path(rds_dir, "10_hdwgcna_input_obj.rds"))


enableWGCNAThreads(nThreads = 8)

obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "11_after_SetupForWGCNA.rds"))

#    Use network_group = all

obj$network_group <- "all"



obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = "network_group",
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "network_group",
  min_cells = 50,
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "12_after_metacells_total_network.rds"))


obj <- NormalizeMetacells(
  obj,
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "13_after_NormalizeMetacells.rds"))


obj <- SetDatExpr(
  obj,
  group_name = "all",
  group.by = "network_group",
  assay = "RNA",
  layer = "data",
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "14_after_SetDatExpr.rds"))



obj <- TestSoftPowers(
  obj,
  networkType = "signed",
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "15_after_TestSoftPowers.rds"))

# Try exporting soft power table if available
soft_fields <- names(obj@misc$hdWGCNA$leydig_hd)

possible_power_tables <- c("power_table", "soft_power_table", "softPower", "sft")
for (nm in possible_power_tables) {
  if (nm %in% soft_fields) {
    x <- obj@misc$hdWGCNA$leydig_hd[[nm]]
    tryCatch({
      write.csv(as.data.frame(x), file.path(table_dir, paste0("hdwgcna_", nm, ".csv")), row.names = FALSE)
    }, error = function(e) {
    })
  }
}


soft_power <- 7

obj <- ConstructNetwork(
  obj,
  soft_power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.20,
  deepSplit = 3,
  overwrite_tom = TRUE,
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "16_after_ConstructNetwork_softpower7.rds"))


me_harmonized <- TRUE
obj <- tryCatch(
  ModuleEigengenes(
    obj,
    group.by.vars = "group",
    wgcna_name = "leydig_hd"
  ),
  error = function(e) {
    me_harmonized <<- FALSE
    ModuleEigengenes(
      obj,
      group.by.vars = NULL,
      wgcna_name = "leydig_hd"
    )
  }
)

saveRDS(obj, file.path(rds_dir, "17_after_ModuleEigengenes.rds"))


obj <- ModuleConnectivity(
  obj,
  group.by = "network_group",
  group_name = "all",
  harmonized = me_harmonized,
  assay = "RNA",
  layer = "data",
  wgcna_name = "leydig_hd"
)

saveRDS(obj, file.path(rds_dir, "18_after_ModuleConnectivity.rds"))


mods <- GetModules(
  obj,
  wgcna_name = "leydig_hd"
)

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

target_genes <- c("Hmgcs2", "Cyp11a1", "Lhcgr", "Foxo3", "Nr5a1", "Star", "Hsd3b1", "Cyp17a1", "Insl3")

target_mods <- mods %>%
  filter(gene_name %in% target_genes)

write.csv(
  target_mods,
  file.path(table_dir, "hdwgcna_target_gene_modules.csv"),
  row.names = FALSE
)


# Export key module genes if present
mods_keep <- intersect(c("tan", "greenyellow", "pink", "green"), unique(mods$module))


for (m in mods_keep) {
  df_m <- mods %>% filter(module == m)
  write.csv(
    df_m,
    file.path(table_dir, paste0("hdwgcna_module_", m, "_genes.csv")),
    row.names = FALSE
  )
}


tryCatch({
  mes <- GetMEs(obj, harmonized = TRUE, wgcna_name = "leydig_hd")
  write.csv(
    mes,
    file.path(table_dir, "hdwgcna_module_eigengenes_harmonized.csv"),
    row.names = TRUE
  )
}, error = function(e) {
})

tryCatch({
  mes_raw <- GetMEs(obj, harmonized = FALSE, wgcna_name = "leydig_hd")
  write.csv(
    mes_raw,
    file.path(table_dir, "hdwgcna_module_eigengenes_raw.csv"),
    row.names = TRUE
  )
}, error = function(e) {
})

saveRDS(obj, file.path(rds_dir, "19_hdwgcna_final_object.rds"))

