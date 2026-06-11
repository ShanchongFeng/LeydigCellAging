#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(DESeq2)
  library(readr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

base_dir <- path.expand(Sys.getenv(
  "GSE303193_EXTERNAL_MOUSE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE303193_EXTERNAL_MOUSE_ROOT")
))

raw_dir <- file.path(base_dir, "raw_data")
table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
rds_dir <- file.path(base_dir, "rds")
tmp_dir <- file.path(base_dir, "tmp")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)


sample_manifest <- data.frame(
  gsm = c(
    "GSM9120206", "GSM9120208", "GSM9120210",
    "GSM9120212", "GSM9120214", "GSM9120216"
  ),
  sample_id = c("M5_1", "M5_2", "M5_3", "M20_1", "M20_2", "M20_3"),
  group = c("Young", "Young", "Young", "Aged", "Aged", "Aged"),
  prefix = c(
    "GSM9120206_M5_1",
    "GSM9120208_M5_2",
    "GSM9120210_M5_3",
    "GSM9120212_M20_1",
    "GSM9120214_M20_2",
    "GSM9120216_M20_3"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  sample_manifest,
  file.path(table_dir, "0_sample_manifest.csv"),
  row.names = FALSE
)


read_one_sample <- function(prefix, sample_id, group) {
  mtx_file <- file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz"))
  barcode_file <- file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_file <- file.path(raw_dir, paste0(prefix, "_features.tsv.gz"))

  if (!file.exists(mtx_file)) stop("Missing mtx: ", mtx_file)
  if (!file.exists(barcode_file)) stop("Missing barcodes: ", barcode_file)
  if (!file.exists(feature_file)) stop("Missing features: ", feature_file)


  mat <- ReadMtx(
    mtx = mtx_file,
    cells = barcode_file,
    features = feature_file,
    feature.column = 2,
    unique.features = TRUE
  )

  colnames(mat) <- paste0(sample_id, "_", colnames(mat))

  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj$sample_id <- sample_id
  obj$group <- group
  obj$age_month <- ifelse(group == "Young", 5, 20)


  obj
}

obj_list <- vector("list", nrow(sample_manifest))

for (i in seq_len(nrow(sample_manifest))) {
  obj_list[[i]] <- read_one_sample(
    prefix = sample_manifest$prefix[i],
    sample_id = sample_manifest$sample_id[i],
    group = sample_manifest$group[i]
  )
}

names(obj_list) <- sample_manifest$sample_id


obj <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  project = "GSE303193_mouse_testis"
)

DefaultAssay(obj) <- "RNA"

# Preserve per-sample Seurat v5 count layers through normalization and
# variable-feature selection, matching the retained manuscript object.


saveRDS(obj, file.path(rds_dir, "1_raw_merged_seurat.rds"))

qc_table <- obj@meta.data %>%
  group_by(sample_id, group) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(qc_table, file.path(table_dir, "1_raw_qc_by_sample.csv"), row.names = FALSE)


p_qc_raw <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "1_raw_qc_violin_by_sample.pdf"), p_qc_raw, width = 12, height = 5)
ggsave(file.path(fig_dir, "1_raw_qc_violin_by_sample.png"), p_qc_raw, width = 12, height = 5, dpi = 300)

# Main QC threshold
obj <- subset(
  obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 7000 &
           percent.mt < 10
)


qc_table_filtered <- obj@meta.data %>%
  group_by(sample_id, group) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(qc_table_filtered, file.path(table_dir, "2_filtered_qc_by_sample.csv"), row.names = FALSE)

saveRDS(obj, file.path(rds_dir, "2_filtered_preprocess_input.rds"))


set.seed(12345)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# Scale variable features only to reduce memory pressure
obj <- ScaleData(obj, features = VariableFeatures(obj))

obj <- RunPCA(obj)

p_elbow <- ElbowPlot(obj, ndims = 30)
ggsave(file.path(fig_dir, "2_elbow.pdf"), p_elbow, width = 6, height = 4)
ggsave(file.path(fig_dir, "2_elbow.png"), p_elbow, width = 6, height = 4, dpi = 300)

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)
obj <- RunUMAP(obj, dims = 1:20)

saveRDS(obj, file.path(rds_dir, "3_clustered_seurat_res0.4.rds"))


# Preserve the layered clustered object above for exact reconstruction checks,
# then join layers for downstream summaries and scoring.
obj <- JoinLayers(obj)

cluster_counts <- as.data.frame(table(obj$seurat_clusters))
colnames(cluster_counts) <- c("cluster", "n_cells")
write.csv(cluster_counts, file.path(table_dir, "3_cluster_cell_counts.csv"), row.names = FALSE)

group_cluster_counts <- as.data.frame(table(obj$group, obj$seurat_clusters))
colnames(group_cluster_counts) <- c("group", "cluster", "n_cells")
write.csv(group_cluster_counts, file.path(table_dir, "4_group_by_cluster_cell_counts.csv"), row.names = FALSE)

p_umap_cluster <- DimPlot(obj, reduction = "umap", label = TRUE) +
  ggtitle("GSE303193 UMAP by cluster")

p_umap_group <- DimPlot(obj, reduction = "umap", group.by = "group") +
  ggtitle("GSE303193 UMAP by group")

p_umap_sample <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
  ggtitle("GSE303193 UMAP by sample")

ggsave(file.path(fig_dir, "3_umap_by_cluster.pdf"), p_umap_cluster, width = 7, height = 6)
ggsave(file.path(fig_dir, "3_umap_by_cluster.png"), p_umap_cluster, width = 7, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "4_umap_by_group.pdf"), p_umap_group, width = 7, height = 6)
ggsave(file.path(fig_dir, "4_umap_by_group.png"), p_umap_group, width = 7, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "5_umap_by_sample.pdf"), p_umap_sample, width = 8, height = 6)
ggsave(file.path(fig_dir, "5_umap_by_sample.png"), p_umap_sample, width = 8, height = 6, dpi = 300)

leydig_markers <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1")
target_genes <- c("Hmgcs2", "Cyp11a1", "Cyp17a1", "Star", "Hsd3b1", "Lhcgr", "Nr5a1", "Foxo3", "Insl3")

marker_presence <- data.frame(
  gene = unique(c(leydig_markers, target_genes)),
  present = unique(c(leydig_markers, target_genes)) %in% rownames(obj)
)

write.csv(marker_presence, file.path(table_dir, "5_marker_presence.csv"), row.names = FALSE)


markers_use <- leydig_markers[leydig_markers %in% rownames(obj)]
targets_use <- target_genes[target_genes %in% rownames(obj)]

p_marker_feature <- FeaturePlot(obj, features = markers_use, ncol = 3)
ggsave(file.path(fig_dir, "6_leydig_marker_featureplot.pdf"), p_marker_feature, width = 12, height = 8)
ggsave(file.path(fig_dir, "6_leydig_marker_featureplot.png"), p_marker_feature, width = 12, height = 8, dpi = 300)

p_marker_dot <- DotPlot(obj, features = markers_use) + RotatedAxis()
ggsave(file.path(fig_dir, "7_leydig_marker_dotplot.pdf"), p_marker_dot, width = 10, height = 5)
ggsave(file.path(fig_dir, "7_leydig_marker_dotplot.png"), p_marker_dot, width = 10, height = 5, dpi = 300)

avg_marker <- AverageExpression(
  obj,
  features = markers_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)

avg_marker_mat <- avg_marker$RNA

write.csv(
  avg_marker_mat,
  file.path(table_dir, "6_leydig_marker_average_expression_by_cluster.csv")
)

# Leydig marker score by cluster
marker_score_df <- data.frame(
  cluster = colnames(avg_marker_mat),
  leydig_marker_mean = colMeans(avg_marker_mat[markers_use, , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(leydig_marker_mean))

write.csv(
  marker_score_df,
  file.path(table_dir, "7_leydig_marker_score_by_cluster.csv"),
  row.names = FALSE
)


Idents(obj) <- "seurat_clusters"

cluster_ids <- levels(Idents(obj))
old_leydig_clusters <- c("16", "17", "19")

if (all(old_leydig_clusters %in% cluster_ids)) {
  leydig_clusters_use <- old_leydig_clusters
  boundary_mode <- "old_boundary_16_17_19"
} else {
  leydig_clusters_use <- marker_score_df$cluster[1:min(3, nrow(marker_score_df))]
  boundary_mode <- paste0("auto_top_marker_clusters_", paste(leydig_clusters_use, collapse = "_"))
  warning(
    "Old GSE303193 clusters 16/17/19 not all found. Using auto top marker clusters: ",
    paste(leydig_clusters_use, collapse = ", ")
  )
}


leydig_obj <- subset(obj, idents = leydig_clusters_use)

saveRDS(leydig_obj, file.path(rds_dir, "4_leydig_subset_used.rds"))

# compatibility name if old boundary is reproduced
if (boundary_mode == "old_boundary_16_17_19") {
  saveRDS(leydig_obj, file.path(rds_dir, "4_leydig_cluster16_17_19.rds"))
}


leydig_counts <- as.data.frame(table(leydig_obj$sample_id, leydig_obj$group, leydig_obj$seurat_clusters))
colnames(leydig_counts) <- c("sample_id", "group", "cluster", "n_cells")
write.csv(leydig_counts, file.path(table_dir, "8_leydig_subset_cell_counts.csv"), row.names = FALSE)

p_leydig_group <- DimPlot(leydig_obj, reduction = "umap", group.by = "group") +
  ggtitle(paste0("GSE303193 Leydig subset: ", paste(leydig_clusters_use, collapse = "+")))

p_leydig_cluster <- DimPlot(leydig_obj, reduction = "umap", label = TRUE) +
  ggtitle("GSE303193 Leydig subset by cluster")

ggsave(file.path(fig_dir, "8_leydig_subset_umap_by_group.pdf"), p_leydig_group, width = 6, height = 5)
ggsave(file.path(fig_dir, "8_leydig_subset_umap_by_group.png"), p_leydig_group, width = 6, height = 5, dpi = 300)

ggsave(file.path(fig_dir, "9_leydig_subset_umap_by_cluster.pdf"), p_leydig_cluster, width = 6, height = 5)
ggsave(file.path(fig_dir, "9_leydig_subset_umap_by_cluster.png"), p_leydig_cluster, width = 6, height = 5, dpi = 300)

p_target_dot <- DotPlot(leydig_obj, features = targets_use, group.by = "group") +
  RotatedAxis() +
  ggtitle("GSE303193 Leydig target genes by group")

ggsave(file.path(fig_dir, "10_leydig_target_gene_dotplot_by_group.pdf"), p_target_dot, width = 10, height = 4)
ggsave(file.path(fig_dir, "10_leydig_target_gene_dotplot_by_group.png"), p_target_dot, width = 10, height = 4, dpi = 300)

module_file <- path.expand(Sys.getenv(
  "HDWGCNA_MODULE_ASSIGNMENT_CSV",
  unset = file.path(
    path.expand(Sys.getenv("HDWGCNA_ROOT", unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT"))),
    "tables", "hdwgcna_module_assignment.csv"
  )
))

if (!file.exists(module_file)) {
  stop("Missing current hdWGCNA module assignment: ", module_file)
}

mods <- read.csv(module_file, check.names = FALSE)

blue_genes <- unique(mods$gene_name[mods$module == "blue"])
brown_genes <- unique(mods$gene_name[mods$module == "brown"])
turquoise_genes <- unique(mods$gene_name[mods$module == "turquoise"])

blue_use <- intersect(blue_genes, rownames(leydig_obj))
brown_use <- intersect(brown_genes, rownames(leydig_obj))
turquoise_use <- intersect(turquoise_genes, rownames(leydig_obj))

module_overlap <- data.frame(
  module = c("blue", "brown", "turquoise"),
  mouse_module_n = c(length(blue_genes), length(brown_genes), length(turquoise_genes)),
  overlap_n = c(length(blue_use), length(brown_use), length(turquoise_use)),
  stringsAsFactors = FALSE
)

write.csv(module_overlap, file.path(table_dir, "9_current_module_overlap_in_gse303193_leydig.csv"), row.names = FALSE)


if (length(blue_use) >= 10) {
  leydig_obj <- AddModuleScore(leydig_obj, features = list(blue_use), name = "blue_module")
}

if (length(brown_use) >= 10) {
  leydig_obj <- AddModuleScore(leydig_obj, features = list(brown_use), name = "brown_module")
}

if (length(turquoise_use) >= 10) {
  leydig_obj <- AddModuleScore(leydig_obj, features = list(turquoise_use), name = "turquoise_module")
}

saveRDS(leydig_obj, file.path(rds_dir, "5_leydig_subset_with_current_module_scores.rds"))

module_score_vars <- intersect(
  c("blue_module1", "brown_module1", "turquoise_module1"),
  colnames(leydig_obj@meta.data)
)

if (length(module_score_vars) > 0) {
  p_module_violin <- VlnPlot(
    leydig_obj,
    features = module_score_vars,
    group.by = "group",
    pt.size = 0,
    ncol = length(module_score_vars)
  )

  ggsave(file.path(fig_dir, "11_current_module_scores_violin_by_group.pdf"), p_module_violin, width = 12, height = 4)
  ggsave(file.path(fig_dir, "11_current_module_scores_violin_by_group.png"), p_module_violin, width = 12, height = 4, dpi = 300)

  score_by_sample <- FetchData(
    leydig_obj,
    vars = c("sample_id", "group", module_score_vars)
  ) %>%
    group_by(sample_id, group) %>%
    summarise(
      across(all_of(module_score_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  write.csv(
    score_by_sample,
    file.path(table_dir, "10_current_module_score_by_sample.csv"),
    row.names = FALSE
  )


  for (v in module_score_vars) {
    p <- ggplot(score_by_sample, aes(x = group, y = .data[[v]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, height = 0, size = 2) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0("GSE303193 ", v, " by sample"),
        x = NULL,
        y = "sample mean score"
      )

    ggsave(
      file.path(fig_dir, paste0("12_sample_box_", v, ".pdf")),
      p,
      width = 5,
      height = 4
    )

    ggsave(
      file.path(fig_dir, paste0("12_sample_box_", v, ".png")),
      p,
      width = 5,
      height = 4,
      dpi = 300
    )
  }
}


counts <- GetAssayData(leydig_obj, assay = "RNA", layer = "counts")
sample_ids <- leydig_obj$sample_id
sample_levels <- unique(sample_ids)

pb_mat <- sapply(sample_levels, function(sid) {
  Matrix::rowSums(counts[, sample_ids == sid, drop = FALSE])
})

pb_mat <- as.matrix(pb_mat)
storage.mode(pb_mat) <- "integer"

sample_info <- unique(leydig_obj@meta.data[, c("sample_id", "group")])
sample_info <- sample_info[match(colnames(pb_mat), sample_info$sample_id), ]
rownames(sample_info) <- sample_info$sample_id
sample_info$group <- factor(sample_info$group, levels = c("Young", "Aged"))


write.csv(pb_mat, file.path(table_dir, "11_leydig_pseudobulk_counts_by_sample.csv"))
write.csv(sample_info, file.path(table_dir, "12_leydig_pseudobulk_sample_info.csv"), row.names = FALSE)

keep_genes <- rowSums(pb_mat >= 10) >= 3

dds <- DESeqDataSetFromMatrix(
  countData = pb_mat[keep_genes, ],
  colData = sample_info,
  design = ~ group
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "Aged", "Young"))
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- res_df[order(is.na(res_df$padj), res_df$padj), ]

write.csv(
  res_df,
  file.path(table_dir, "13_leydig_pseudobulk_DESeq2_Aged_vs_Young.csv"),
  row.names = FALSE
)

target_res <- res_df %>%
  filter(Gene %in% targets_use) %>%
  arrange(match(Gene, targets_use))

write.csv(
  target_res,
  file.path(table_dir, "14_leydig_target_genes_pseudobulk_DESeq2.csv"),
  row.names = FALSE
)


vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE303193 Leydig pseudobulk PCA",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

ggsave(file.path(fig_dir, "13_leydig_pseudobulk_PCA.pdf"), p_pca, width = 6, height = 5)
ggsave(file.path(fig_dir, "13_leydig_pseudobulk_PCA.png"), p_pca, width = 6, height = 5, dpi = 300)

