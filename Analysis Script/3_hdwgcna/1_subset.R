#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

base_dir <- path.expand(Sys.getenv(
  "HDWGCNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT")
))
input_rds <- path.expand(Sys.getenv(
  "GSE270931_HDWGCNA_INPUT_RDS",
  unset = file.path(
    path.expand(Sys.getenv(
      "GSE270931_SCRNA_ROOT",
      unset = file.path(Sys.getenv("HOME"), "GSE270931_SCRNA_ROOT")
    )),
    "rds", "GSE270931_filtered_clustered_seurat.rds"
  )
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
  stop("Input Seurat RDS not found: ", input_rds)
}

sce <- readRDS(input_rds)
DefaultAssay(sce) <- "RNA"


# Save input copy
saveRDS(sce, file.path(rds_dir, "0_input_filtered_clustered_seurat.rds"))


sce_harmony <- sce

sce_harmony <- RunHarmony(
  object = sce_harmony,
  group.by.vars = "group",
  reduction.use = "pca",
  dims.use = 1:20
)

sce_harmony <- FindNeighbors(
  sce_harmony,
  reduction = "harmony",
  dims = 1:20
)

sce_harmony <- FindClusters(
  sce_harmony,
  resolution = 0.4
)

sce_harmony <- RunUMAP(
  sce_harmony,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap.harmony",
  reduction.key = "hUMAP_"
)

saveRDS(sce_harmony, file.path(rds_dir, "1_harmony_clustered_seurat.rds"))


cluster_table <- as.data.frame(table(sce_harmony$seurat_clusters))
colnames(cluster_table) <- c("cluster", "n_cells")
write.csv(cluster_table, file.path(table_dir, "1_harmony_cluster_cell_counts.csv"), row.names = FALSE)

group_cluster_table <- as.data.frame(table(sce_harmony$group, sce_harmony$seurat_clusters))
colnames(group_cluster_table) <- c("group", "cluster", "n_cells")
write.csv(group_cluster_table, file.path(table_dir, "2_harmony_group_by_cluster_cell_counts.csv"), row.names = FALSE)

p_harmony_cluster <- DimPlot(
  sce_harmony,
  reduction = "umap.harmony",
  label = TRUE
) + ggtitle("Harmony UMAP by cluster")

p_harmony_group <- DimPlot(
  sce_harmony,
  reduction = "umap.harmony",
  group.by = "group"
) + ggtitle("Harmony UMAP by group")

ggsave(file.path(fig_dir, "1_harmony_umap_by_cluster.pdf"), p_harmony_cluster, width = 7, height = 6)
ggsave(file.path(fig_dir, "1_harmony_umap_by_cluster.png"), p_harmony_cluster, width = 7, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "2_harmony_umap_by_group.pdf"), p_harmony_group, width = 7, height = 6)
ggsave(file.path(fig_dir, "2_harmony_umap_by_group.png"), p_harmony_group, width = 7, height = 6, dpi = 300)

leydig_markers <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1")
target_genes <- c("Hmgcs2", "Lhcgr", "Star", "Cyp11a1", "Hsd3b1", "Nr5a1", "Foxo3")

genes_marker_use <- leydig_markers[leydig_markers %in% rownames(sce_harmony)]
genes_target_use <- target_genes[target_genes %in% rownames(sce_harmony)]

marker_presence <- data.frame(
  gene = unique(c(leydig_markers, target_genes)),
  present = unique(c(leydig_markers, target_genes)) %in% rownames(sce_harmony)
)

write.csv(marker_presence, file.path(table_dir, "3_harmony_marker_presence.csv"), row.names = FALSE)


p_marker_feature <- FeaturePlot(
  sce_harmony,
  features = genes_marker_use,
  reduction = "umap.harmony",
  ncol = 3
)

p_marker_dot <- DotPlot(
  sce_harmony,
  features = genes_marker_use
) + RotatedAxis() + ggtitle("Leydig markers after Harmony")

ggsave(file.path(fig_dir, "3_harmony_leydig_marker_featureplot.pdf"), p_marker_feature, width = 12, height = 8)
ggsave(file.path(fig_dir, "3_harmony_leydig_marker_featureplot.png"), p_marker_feature, width = 12, height = 8, dpi = 300)

ggsave(file.path(fig_dir, "4_harmony_leydig_marker_dotplot.pdf"), p_marker_dot, width = 9, height = 4)
ggsave(file.path(fig_dir, "4_harmony_leydig_marker_dotplot.png"), p_marker_dot, width = 9, height = 4, dpi = 300)

avg_marker <- AverageExpression(
  sce_harmony,
  features = genes_marker_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)

write.csv(
  avg_marker$RNA,
  file.path(table_dir, "4_harmony_leydig_marker_average_expression_by_cluster.csv")
)

#    Old log used cluster 0 + 15.
#    Here we use 0 + 15 if present.
Idents(sce_harmony) <- "seurat_clusters"

cluster_ids <- levels(Idents(sce_harmony))
harmony_leydig_clusters <- c("0", "15")

if (!all(harmony_leydig_clusters %in% cluster_ids)) {
  stop(
    "Expected Harmony Leydig clusters 0 and/or 15 not found. Existing clusters: ",
    paste(cluster_ids, collapse = ", "),
    "\nCheck figures/4_harmony_leydig_marker_dotplot.pdf before modifying boundary."
  )
}

leydig_harmony <- subset(sce_harmony, idents = harmony_leydig_clusters)

saveRDS(leydig_harmony, file.path(rds_dir, "2_leydig_harmony_subset_0_15.rds"))

# Compatibility name for downstream hdWGCNA scripts
saveRDS(leydig_harmony, file.path(rds_dir, "2_leydig_harmony_subset.rds"))


leydig_counts <- as.data.frame(table(leydig_harmony$group, leydig_harmony$seurat_clusters))
colnames(leydig_counts) <- c("group", "cluster", "n_cells")
write.csv(leydig_counts, file.path(table_dir, "5_harmony_leydig_subset_0_15_cell_counts.csv"), row.names = FALSE)

p_leydig_group <- DimPlot(
  leydig_harmony,
  reduction = "umap.harmony",
  group.by = "group"
) + ggtitle("Harmony Leydig subset 0+15 by group")

p_leydig_cluster <- DimPlot(
  leydig_harmony,
  reduction = "umap.harmony",
  label = TRUE
) + ggtitle("Harmony Leydig subset 0+15 by cluster")

ggsave(file.path(fig_dir, "5_harmony_leydig_subset_0_15_umap_by_group.pdf"), p_leydig_group, width = 6, height = 5)
ggsave(file.path(fig_dir, "5_harmony_leydig_subset_0_15_umap_by_group.png"), p_leydig_group, width = 6, height = 5, dpi = 300)

ggsave(file.path(fig_dir, "6_harmony_leydig_subset_0_15_umap_by_cluster.pdf"), p_leydig_cluster, width = 6, height = 5)
ggsave(file.path(fig_dir, "6_harmony_leydig_subset_0_15_umap_by_cluster.png"), p_leydig_cluster, width = 6, height = 5, dpi = 300)

genes_use <- genes_target_use[genes_target_use %in% rownames(leydig_harmony)]

p_target_dot <- DotPlot(
  leydig_harmony,
  features = genes_use,
  group.by = "group"
) + RotatedAxis() + ggtitle("Harmony Leydig target genes by group")

ggsave(file.path(fig_dir, "7_harmony_leydig_target_gene_dotplot_by_group.pdf"), p_target_dot, width = 9, height = 4)
ggsave(file.path(fig_dir, "7_harmony_leydig_target_gene_dotplot_by_group.png"), p_target_dot, width = 9, height = 4, dpi = 300)

p_target_feature <- FeaturePlot(
  leydig_harmony,
  features = genes_use,
  reduction = "umap.harmony",
  ncol = 3
)

ggsave(file.path(fig_dir, "8_harmony_leydig_target_gene_featureplot.pdf"), p_target_feature, width = 12, height = 8)
ggsave(file.path(fig_dir, "8_harmony_leydig_target_gene_featureplot.png"), p_target_feature, width = 12, height = 8, dpi = 300)

avg_target <- AverageExpression(
  leydig_harmony,
  features = genes_use,
  group.by = "group",
  assays = "RNA",
  layer = "data"
)

write.csv(
  avg_target$RNA,
  file.path(table_dir, "6_harmony_leydig_target_average_expression_by_group.csv")
)

Idents(leydig_harmony) <- "group"

deg_harmony <- FindMarkers(
  leydig_harmony,
  ident.1 = "Aged",
  ident.2 = "Young",
  logfc.threshold = 0,
  min.pct = 0.1
)

deg_harmony$Gene <- rownames(deg_harmony)

write.csv(
  deg_harmony,
  file.path(table_dir, "7_harmony_leydig_Aged_vs_Young_DEG_all.csv"),
  row.names = FALSE
)

target_deg <- deg_harmony %>%
  filter(Gene %in% genes_use)

write.csv(
  target_deg,
  file.path(table_dir, "8_harmony_leydig_target_genes_in_DEG.csv"),
  row.names = FALSE
)


pb_group <- AggregateExpression(
  leydig_harmony,
  features = genes_use,
  group.by = "group",
  assays = "RNA",
  return.seurat = FALSE
)

write.csv(
  pb_group$RNA,
  file.path(table_dir, "9_harmony_leydig_group_level_aggregated_expression.csv")
)

