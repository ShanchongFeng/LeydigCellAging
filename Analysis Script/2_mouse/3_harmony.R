#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

ROOT <- path.expand(Sys.getenv(
  "GSE270931_SCRNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE270931_SCRNA_ROOT")
))
DIRS <- list(
  rds = file.path(ROOT, "rds"),
  tables = file.path(ROOT, "tables"),
  figures = file.path(ROOT, "figures"),
  logs = file.path(ROOT, "logs"),
  metadata = file.path(ROOT, "metadata"),
  session = file.path(ROOT, "session"),
  docs = file.path(ROOT, "docs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

session_file <- file.path(DIRS$session, "3_harmony_sessionInfo.txt")
on.exit({
  writeLines(capture.output(sessionInfo()), session_file)
}, add = TRUE)

input_rds <- file.path(DIRS$rds, "GSE270931_filtered_clustered_seurat.rds")
if (!file.exists(input_rds)) {
  stop("Missing input from main GSE270931 rerun: ", input_rds)
}


write.csv(
  data.frame(
    role = "filtered_clustered_seurat_from_locked_main",
    path = input_rds,
    exists = file.exists(input_rds)
  ),
  file.path(DIRS$metadata, "GSE270931_harmony_input_paths.csv"),
  row.names = FALSE
)

params <- data.frame(
  parameter = c(
    "seed", "integration_variable", "reduction.use", "dims.use",
    "FindNeighbors.reduction", "FindNeighbors.dims",
    "FindClusters.resolution", "RunUMAP.reduction", "RunUMAP.dims",
    "historical_harmony_leydig_boundary", "boundary_fallback",
    "FindMarkers.ident1", "FindMarkers.ident2", "FindMarkers.logfc.threshold", "FindMarkers.min.pct"
  ),
  value = c(
    "12345", "group", "pca", "1:20",
    "harmony", "1:20",
    "0.4", "harmony", "1:20",
    "clusters 0+15", "top two clusters by Leydig marker average expression if 0+15 missing",
    "Aged", "Young", "0", "0.1"
  )
)
write.csv(params, file.path(DIRS$metadata, "GSE270931_harmony_analysis_parameters.csv"), row.names = FALSE)

sce <- readRDS(input_rds)
DefaultAssay(sce) <- "RNA"


sce_harmony <- RunHarmony(
  object = sce,
  group.by.vars = "group",
  reduction.use = "pca",
  dims.use = 1:20
)

sce_harmony <- FindNeighbors(sce_harmony, reduction = "harmony", dims = 1:20)
sce_harmony <- FindClusters(sce_harmony, resolution = 0.4)
sce_harmony <- RunUMAP(
  sce_harmony,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap.harmony",
  reduction.key = "hUMAP_"
)

saveRDS(sce_harmony, file.path(DIRS$rds, "GSE270931_harmony_object.rds"))

cluster_table <- as.data.frame(table(sce_harmony$seurat_clusters))
colnames(cluster_table) <- c("cluster", "n_cells")
write.csv(cluster_table, file.path(DIRS$tables, "GSE270931_harmony_cluster_cell_counts.csv"), row.names = FALSE)

group_cluster <- as.data.frame(table(sce_harmony$group, sce_harmony$seurat_clusters))
colnames(group_cluster) <- c("group", "cluster", "n_cells")
group_cluster <- group_cluster %>%
  group_by(group) %>%
  mutate(proportion_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(proportion_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(group_cluster, file.path(DIRS$tables, "GSE270931_harmony_group_cluster_composition.csv"), row.names = FALSE)

p_harmony <- DimPlot(sce_harmony, reduction = "umap.harmony", label = TRUE) +
  ggtitle("GSE270931 Harmony UMAP")
ggsave(file.path(DIRS$figures, "GSE270931_harmony_umap.pdf"), p_harmony, width = 7, height = 6)

leydig_markers <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Hmgcs2")
target_genes <- c("Hmgcs2", "Cyp11a1", "Lhcgr", "Star", "Hsd3b1", "Hsd3b", "Foxo3", "Nr5a1", "Hdac1")
genes_marker_use <- leydig_markers[leydig_markers %in% rownames(sce_harmony)]
genes_target_use <- target_genes[target_genes %in% rownames(sce_harmony)]

avg_marker <- AverageExpression(
  sce_harmony,
  features = genes_marker_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA
avg_marker_df <- as.data.frame(avg_marker)
avg_marker_df$gene <- rownames(avg_marker_df)
avg_marker_df <- avg_marker_df %>% relocate(gene)
write.csv(avg_marker_df, file.path(DIRS$tables, "GSE270931_harmony_marker_average_expression_by_cluster.csv"), row.names = FALSE)

marker_by_cluster <- as.data.frame(t(as.matrix(avg_marker)))
marker_by_cluster$cluster <- rownames(marker_by_cluster)
marker_by_cluster$leydig_marker_score <- rowMeans(marker_by_cluster[, intersect(genes_marker_use, colnames(marker_by_cluster)), drop = FALSE], na.rm = TRUE)

Idents(sce_harmony) <- "seurat_clusters"
clusters_existing <- levels(Idents(sce_harmony))
historical <- c("0", "15")
if (all(historical %in% clusters_existing)) {
  leydig_clusters <- historical
  boundary_method <- "historical_0_15"
} else {
  leydig_clusters <- marker_by_cluster %>%
    arrange(desc(leydig_marker_score)) %>%
    slice_head(n = 2) %>%
    pull(cluster)
  boundary_method <- "fallback_top2_leydig_marker_score"
  warning("Historical Harmony clusters 0+15 were not both present; using top two marker-score clusters: ",
          paste(leydig_clusters, collapse = "+"))
}

write.csv(
  data.frame(
    boundary = "harmony_leydig",
    method = boundary_method,
    clusters = paste(leydig_clusters, collapse = "+"),
    stringsAsFactors = FALSE
  ),
  file.path(DIRS$metadata, "GSE270931_harmony_leydig_boundary.csv"),
  row.names = FALSE
)

leydig_harmony <- subset(sce_harmony, idents = leydig_clusters)
saveRDS(leydig_harmony, file.path(DIRS$rds, "GSE270931_harmony_leydig_subset.rds"))

leydig_counts <- as.data.frame(table(leydig_harmony$group, leydig_harmony$seurat_clusters))
colnames(leydig_counts) <- c("group", "cluster", "n_cells")
leydig_counts <- leydig_counts %>%
  group_by(group) %>%
  mutate(proportion_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(proportion_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(leydig_counts, file.path(DIRS$tables, "GSE270931_harmony_leydig_group_composition.csv"), row.names = FALSE)

Idents(leydig_harmony) <- "group"
deg <- FindMarkers(
  leydig_harmony,
  ident.1 = "Aged",
  ident.2 = "Young",
  logfc.threshold = 0,
  min.pct = 0.1
)
deg$gene <- rownames(deg)
deg <- deg %>% relocate(gene)
target_deg <- deg %>%
  filter(gene %in% genes_target_use) %>%
  mutate(interpretation_note = "Harmony cell-level technical robustness only; not an independent cohort")
write.csv(target_deg, file.path(DIRS$tables, "GSE270931_harmony_target_gene_DEG.csv"), row.names = FALSE)
write.csv(deg, file.path(DIRS$tables, "GSE270931_harmony_all_DEG.csv"), row.names = FALSE)

set.seed(12345)
module_sets <- list(
  SteroidScore = c("Star", "Cyp11a1", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1"),
  LipidScore = c("Hmgcs2", "Fdft1", "Sqle", "Hmgcr", "Acat1", "Cyp51"),
  OxphosScore = c("Ndufa1", "Ndufb5", "Cox5a", "Cox6a1", "Atp5f1a", "Uqcr10"),
  FoxoCampScore = c("Foxo3", "Prkaca", "Creb1", "Adcy3", "Pde4d", "Nr4a1")
)

score_params <- data.frame(
  dataset = "GSE270931",
  analysis = "Harmony_control",
  score_name = names(module_sets),
  genes = vapply(module_sets, paste, character(1), collapse = ";"),
  assay = "RNA",
  slot = "data",
  ctrl = 100,
  nbin = 24,
  seed = 12345,
  stringsAsFactors = FALSE
)
write.csv(score_params, file.path(DIRS$metadata, "module_score_parameters.csv"), row.names = FALSE)

for (nm in names(module_sets)) {
  genes_use <- module_sets[[nm]][module_sets[[nm]] %in% rownames(leydig_harmony)]
  if (length(genes_use) == 0) {
    warning("No genes available for module score: ", nm)
    next
  }
  set.seed(12345)
  leydig_harmony <- AddModuleScore(
    leydig_harmony,
    features = list(genes_use),
    name = nm,
    assay = "RNA",
    slot = "data",
    ctrl = 100,
    nbin = 24,
    seed = 12345
  )
}
saveRDS(leydig_harmony, file.path(DIRS$rds, "GSE270931_harmony_leydig_subset_with_module_scores.rds"))

score_cols <- intersect(paste0(names(module_sets), "1"), colnames(leydig_harmony@meta.data))
score_long <- leydig_harmony@meta.data %>%
  as.data.frame() %>%
  select(group, seurat_clusters, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "module_score", values_to = "score")

score_summary <- score_long %>%
  group_by(module_score, group) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

pvals <- score_long %>%
  group_by(module_score) %>%
  summarise(
    p_value = suppressWarnings(wilcox.test(score ~ group)$p.value),
    interpretation_note = "cell-level Wilcoxon; technical score trend only",
    .groups = "drop"
  )

write.csv(
  left_join(score_summary, pvals, by = "module_score"),
  file.path(DIRS$tables, "GSE270931_harmony_module_score_summary.csv"),
  row.names = FALSE
)


writeLines(capture.output(sessionInfo()), session_file)
