#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(readr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

base_dir <- path.expand(Sys.getenv(
  "GSE182786_HUMAN_AGING_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE182786_HUMAN_AGING_ROOT")
))
raw_dir <- file.path(base_dir, "raw_data")
table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
rds_dir <- file.path(base_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


read_processed_matrix <- function(file, group_label) {

  if (!file.exists(file)) {
    stop("Missing matrix file: ", file)
  }

  # GEO processed matrix: first row = cell barcodes, first column in data rows = gene names.
  mat_df <- read.table(
    gzfile(file),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )


  mat <- as.matrix(mat_df)
  rm(mat_df)
  gc()

  storage.mode(mat) <- "integer"

  mat <- Matrix(mat, sparse = TRUE)


  mat
}

young_file <- file.path(raw_dir, "GSE182786_Young_matrix.txt.gz")
older1_file <- file.path(raw_dir, "GSE182786_Older_Group1_matrix.txt.gz")
older2_file <- file.path(raw_dir, "GSE182786_Older_Group2_matrix.txt.gz")

young_mat <- read_processed_matrix(young_file, "Young")
older1_mat <- read_processed_matrix(older1_file, "Aged_Group1")
older2_mat <- read_processed_matrix(older2_file, "Aged_Group2")

common_genes <- Reduce(intersect, list(
  rownames(young_mat),
  rownames(older1_mat),
  rownames(older2_mat)
))


young_mat <- young_mat[common_genes, ]
older1_mat <- older1_mat[common_genes, ]
older2_mat <- older2_mat[common_genes, ]

all_mat <- cbind(young_mat, older1_mat, older2_mat)

rm(young_mat, older1_mat, older2_mat)
gc()


obj <- CreateSeuratObject(
  counts = all_mat,
  project = "GSE182786_human_testis",
  min.cells = 3,
  min.features = 200
)

rm(all_mat)
gc()

DefaultAssay(obj) <- "RNA"

obj$cell_id <- colnames(obj)
obj$sample_id <- sub("\\..*$", "", colnames(obj))
obj$group <- ifelse(grepl("^Young", obj$sample_id), "Young", "Aged")

# Map GEO sample names
sample_manifest <- data.frame(
  sample_id = c(
    "Young1", "Young2", "Young3", "Young4",
    "Older1", "Older2", "Older3", "Older4",
    "Older5", "Older6", "Older7", "Older8"
  ),
  geo_title = c(
    "Young_1", "Young_2", "Young_3", "Young_4",
    "Older_1", "Older_2", "Older_3", "Older_4",
    "Older_5", "Older_6", "Older_7", "Older_8"
  ),
  group = c(rep("Young", 4), rep("Aged", 8)),
  age_class = c(rep("17_22", 4), rep("gt60", 8)),
  stringsAsFactors = FALSE
)

write.csv(sample_manifest, file.path(table_dir, "0_sample_manifest.csv"), row.names = FALSE)


saveRDS(obj, file.path(rds_dir, "1_raw_merged_seurat.rds"))

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

qc_raw <- obj@meta.data %>%
  group_by(sample_id, group) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(qc_raw, file.path(table_dir, "1_raw_qc_by_sample.csv"), row.names = FALSE)


p_qc <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample_id",
  ncol = 3,
  pt.size = 0
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "1_raw_qc_violin_by_sample.pdf"), p_qc, width = 12, height = 5)
ggsave(file.path(fig_dir, "1_raw_qc_violin_by_sample.png"), p_qc, width = 12, height = 5, dpi = 300)

obj <- subset(
  obj,
  subset = nFeature_RNA > 300 &
           nFeature_RNA < 7000 &
           percent.mt < 20
)


qc_filtered <- obj@meta.data %>%
  group_by(sample_id, group) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(qc_filtered, file.path(table_dir, "2_filtered_qc_by_sample.csv"), row.names = FALSE)

saveRDS(obj, file.path(rds_dir, "2_filtered_preprocess_input.rds"))


obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, npcs = 30, features = VariableFeatures(obj))

p_elbow <- ElbowPlot(obj, ndims = 30)
ggsave(file.path(fig_dir, "2_elbow.pdf"), p_elbow, width = 6, height = 4)
ggsave(file.path(fig_dir, "2_elbow.png"), p_elbow, width = 6, height = 4, dpi = 300)

obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)
obj <- RunUMAP(obj, dims = 1:20)

saveRDS(obj, file.path(rds_dir, "3_clustered_seurat_res0.4.rds"))


write.csv(
  as.data.frame(table(obj$seurat_clusters)),
  file.path(table_dir, "3_cluster_cell_counts.csv"),
  row.names = FALSE
)

group_cluster <- as.data.frame(table(obj$group, obj$seurat_clusters))
colnames(group_cluster) <- c("group", "cluster", "n_cells")
write.csv(group_cluster, file.path(table_dir, "4_group_by_cluster_counts.csv"), row.names = FALSE)

p_umap_cluster <- DimPlot(obj, reduction = "umap", label = TRUE) +
  ggtitle("GSE182786 UMAP by cluster")

p_umap_group <- DimPlot(obj, reduction = "umap", group.by = "group") +
  ggtitle("GSE182786 UMAP by group")

p_umap_sample <- DimPlot(obj, reduction = "umap", group.by = "sample_id") +
  ggtitle("GSE182786 UMAP by sample")

ggsave(file.path(fig_dir, "3_umap_by_cluster.pdf"), p_umap_cluster, width = 7, height = 6)
ggsave(file.path(fig_dir, "3_umap_by_cluster.png"), p_umap_cluster, width = 7, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "4_umap_by_group.pdf"), p_umap_group, width = 7, height = 6)
ggsave(file.path(fig_dir, "4_umap_by_group.png"), p_umap_group, width = 7, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "5_umap_by_sample.pdf"), p_umap_sample, width = 9, height = 7)
ggsave(file.path(fig_dir, "5_umap_by_sample.png"), p_umap_sample, width = 9, height = 7, dpi = 300)

leydig_markers <- c("INSL3", "CYP11A1", "STAR", "HSD3B2", "HSD3B1", "CYP17A1", "LHCGR", "NR5A1")
target_genes <- c("HMGCS2", "CYP11A1", "CYP17A1", "STAR", "HSD3B1", "HSD3B2", "LHCGR", "NR5A1", "FOXO3", "INSL3")

marker_presence <- data.frame(
  gene = unique(c(leydig_markers, target_genes)),
  present = unique(c(leydig_markers, target_genes)) %in% rownames(obj)
)

write.csv(marker_presence, file.path(table_dir, "5_marker_presence.csv"), row.names = FALSE)


markers_use <- leydig_markers[leydig_markers %in% rownames(obj)]
targets_use <- target_genes[target_genes %in% rownames(obj)]

p_marker_dot <- DotPlot(obj, features = markers_use) + RotatedAxis()
ggsave(file.path(fig_dir, "6_leydig_marker_dotplot.pdf"), p_marker_dot, width = 10, height = 5)
ggsave(file.path(fig_dir, "6_leydig_marker_dotplot.png"), p_marker_dot, width = 10, height = 5, dpi = 300)

p_marker_feature <- FeaturePlot(obj, features = markers_use, ncol = 3)
ggsave(file.path(fig_dir, "7_leydig_marker_featureplot.pdf"), p_marker_feature, width = 12, height = 8)
ggsave(file.path(fig_dir, "7_leydig_marker_featureplot.png"), p_marker_feature, width = 12, height = 8, dpi = 300)

avg_marker <- AverageExpression(
  obj,
  features = markers_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)

avg_marker_mat <- avg_marker$RNA
write.csv(avg_marker_mat, file.path(table_dir, "6_leydig_marker_average_expression_by_cluster.csv"))

marker_score_df <- data.frame(
  cluster = colnames(avg_marker_mat),
  leydig_marker_mean = colMeans(avg_marker_mat[markers_use, , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(leydig_marker_mean))

write.csv(marker_score_df, file.path(table_dir, "7_leydig_marker_score_by_cluster.csv"), row.names = FALSE)


Idents(obj) <- "seurat_clusters"
cluster_ids <- levels(Idents(obj))

module_file <- path.expand(Sys.getenv(
  "HDWGCNA_MODULE_ASSIGNMENT_CSV",
  unset = file.path(
    path.expand(Sys.getenv("HDWGCNA_ROOT", unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT"))),
    "tables", "hdwgcna_module_assignment.csv"
  )
))
if (!file.exists(module_file)) stop("Missing module file: ", module_file)

mods <- read.csv(module_file, check.names = FALSE)
mods$human_symbol_guess <- toupper(mods$gene_name)

blue_genes <- unique(mods$human_symbol_guess[mods$module == "blue"])
brown_genes <- unique(mods$human_symbol_guess[mods$module == "brown"])
turquoise_genes <- unique(mods$human_symbol_guess[mods$module == "turquoise"])

marker_rank <- marker_score_df
marker_rank$cluster_clean <- sub("^g", "", marker_rank$cluster)

old_core <- c("0", "17")
old_ext <- c("0", "14", "17")
auto_top2 <- marker_rank$cluster_clean[1:min(2, nrow(marker_rank))]
auto_top3 <- marker_rank$cluster_clean[1:min(3, nrow(marker_rank))]
union_core_auto <- unique(c(old_core, auto_top3))
union_ext_auto <- unique(c(old_ext, auto_top3))

boundaries <- list(
  old_core_0_17 = old_core,
  old_ext_0_14_17 = old_ext,
  auto_top2_marker = auto_top2,
  auto_top3_marker = auto_top3,
  union_core_auto = union_core_auto,
  union_ext_auto = union_ext_auto
)


run_one_boundary <- function(boundary_name, clusters_use) {

  clusters_use <- clusters_use[clusters_use %in% cluster_ids]

  if (length(clusters_use) == 0) {
    warning("No valid clusters for boundary: ", boundary_name)
    return(NULL)
  }

  sub_obj <- subset(obj, idents = clusters_use)
  DefaultAssay(sub_obj) <- "RNA"


  saveRDS(sub_obj, file.path(rds_dir, paste0("4_boundary_", boundary_name, ".rds")))

  blue_use <- intersect(blue_genes, rownames(sub_obj))
  brown_use <- intersect(brown_genes, rownames(sub_obj))
  turquoise_use <- intersect(turquoise_genes, rownames(sub_obj))

  overlap_df <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    module = c("blue", "brown", "turquoise"),
    module_n = c(length(blue_genes), length(brown_genes), length(turquoise_genes)),
    overlap_n = c(length(blue_use), length(brown_use), length(turquoise_use)),
    stringsAsFactors = FALSE
  )

  if (length(blue_use) >= 10) {
    sub_obj <- AddModuleScore(sub_obj, features = list(blue_use), name = "blue_module")
  }
  if (length(brown_use) >= 10) {
    sub_obj <- AddModuleScore(sub_obj, features = list(brown_use), name = "brown_module")
  }
  if (length(turquoise_use) >= 10) {
    sub_obj <- AddModuleScore(sub_obj, features = list(turquoise_use), name = "turquoise_module")
  }

  module_score_vars <- intersect(
    c("blue_module1", "brown_module1", "turquoise_module1"),
    colnames(sub_obj@meta.data)
  )

  score_by_sample <- FetchData(
    sub_obj,
    vars = c("sample_id", "group", module_score_vars)
  ) %>%
    group_by(sample_id, group) %>%
    summarise(
      across(all_of(module_score_vars), ~ mean(.x, na.rm = TRUE)),
      n_cells = n(),
      .groups = "drop"
    )

  score_by_sample$boundary <- boundary_name
  score_by_sample$clusters <- paste(clusters_use, collapse = "+")

  wilcox_rows <- list()
  for (v in module_score_vars) {
    wt <- suppressWarnings(wilcox.test(as.formula(paste(v, "~ group")), data = score_by_sample))

    mean_df <- score_by_sample %>%
      group_by(group) %>%
      summarise(
        mean_score = mean(.data[[v]], na.rm = TRUE),
        median_score = median(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      )

    young_mean <- mean_df$mean_score[mean_df$group == "Young"]
    aged_mean <- mean_df$mean_score[mean_df$group == "Aged"]
    young_median <- mean_df$median_score[mean_df$group == "Young"]
    aged_median <- mean_df$median_score[mean_df$group == "Aged"]

    wilcox_rows[[v]] <- data.frame(
      boundary = boundary_name,
      clusters = paste(clusters_use, collapse = "+"),
      module_score = v,
      p_value = wt$p.value,
      young_mean = young_mean,
      aged_mean = aged_mean,
      young_median = young_median,
      aged_median = aged_median,
      direction_aged_vs_young = ifelse(aged_mean > young_mean, "Aged_higher", "Aged_lower"),
      stringsAsFactors = FALSE
    )
  }

  wilcox_df <- bind_rows(wilcox_rows)

  # Pseudobulk counts by sample
  counts <- GetAssayData(sub_obj, assay = "RNA", layer = "counts")
  sample_ids <- sub_obj$sample_id
  sample_levels <- unique(sample_ids)

  pb_mat <- sapply(sample_levels, function(sid) {
    Matrix::rowSums(counts[, sample_ids == sid, drop = FALSE])
  })

  pb_mat <- as.matrix(pb_mat)
  storage.mode(pb_mat) <- "integer"

  sample_info <- unique(sub_obj@meta.data[, c("sample_id", "group")])
  sample_info <- sample_info[match(colnames(pb_mat), sample_info$sample_id), ]
  rownames(sample_info) <- sample_info$sample_id
  sample_info$group <- factor(sample_info$group, levels = c("Young", "Aged"))

  write.csv(pb_mat, file.path(table_dir, paste0("boundary_", boundary_name, "_pseudobulk_counts_by_sample.csv")))
  write.csv(sample_info, file.path(table_dir, paste0("boundary_", boundary_name, "_pseudobulk_sample_info.csv")), row.names = FALSE)

  targets_present <- target_genes[target_genes %in% rownames(pb_mat)]
  target_counts <- pb_mat[targets_present, , drop = FALSE]

  target_count_summary <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    Gene = rownames(target_counts),
    total_count = rowSums(target_counts),
    samples_count_ge1 = rowSums(target_counts >= 1),
    samples_count_ge5 = rowSums(target_counts >= 5),
    samples_count_ge10 = rowSums(target_counts >= 10),
    pass_filter_ge10_in_3samples = rowSums(target_counts >= 10) >= 3,
    stringsAsFactors = FALSE
  )

  keep_genes <- rowSums(pb_mat >= 10) >= 3

  target_res_full <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    Gene = targets_present,
    in_DESeq_filter = targets_present %in% rownames(pb_mat)[keep_genes],
    stringsAsFactors = FALSE
  )

  if (sum(keep_genes) >= 100) {
    dds <- DESeqDataSetFromMatrix(
      countData = pb_mat[keep_genes, ],
      colData = sample_info,
      design = ~ group
    )

    dds <- DESeq(dds, quiet = TRUE)

    res <- results(dds, contrast = c("group", "Aged", "Young"))
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)

    target_res <- res_df %>%
      filter(Gene %in% targets_present) %>%
      select(Gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

    target_res_full <- target_res_full %>%
      left_join(target_res, by = "Gene")

    write.csv(
      res_df[order(is.na(res_df$padj), res_df$padj), ],
      file.path(table_dir, paste0("boundary_", boundary_name, "_DESeq2_Aged_vs_Young_all.csv")),
      row.names = FALSE
    )
  }

  cell_counts <- as.data.frame(table(sub_obj$sample_id, sub_obj$group, sub_obj$seurat_clusters))
  colnames(cell_counts) <- c("sample_id", "group", "cluster", "n_cells")
  cell_counts$boundary <- boundary_name
  cell_counts$clusters <- paste(clusters_use, collapse = "+")

  saveRDS(sub_obj, file.path(rds_dir, paste0("5_boundary_", boundary_name, "_with_module_scores.rds")))

  list(
    overlap = overlap_df,
    score_by_sample = score_by_sample,
    wilcox = wilcox_df,
    target_count = target_count_summary,
    target_res = target_res_full,
    cell_counts = cell_counts
  )
}

res_list <- list()

for (nm in names(boundaries)) {
  res_list[[nm]] <- run_one_boundary(nm, boundaries[[nm]])
}

res_list <- res_list[!sapply(res_list, is.null)]

overlap_all <- bind_rows(lapply(res_list, `[[`, "overlap"))
score_all <- bind_rows(lapply(res_list, `[[`, "score_by_sample"))
wilcox_all <- bind_rows(lapply(res_list, `[[`, "wilcox"))
target_count_all <- bind_rows(lapply(res_list, `[[`, "target_count"))
target_res_all <- bind_rows(lapply(res_list, `[[`, "target_res"))
cell_counts_all <- bind_rows(lapply(res_list, `[[`, "cell_counts"))

write.csv(overlap_all, file.path(table_dir, "8_boundary_module_overlap.csv"), row.names = FALSE)
write.csv(score_all, file.path(table_dir, "9_boundary_module_score_by_sample.csv"), row.names = FALSE)
write.csv(wilcox_all, file.path(table_dir, "10_boundary_module_score_young_vs_aged.csv"), row.names = FALSE)
write.csv(target_count_all, file.path(table_dir, "11_boundary_target_count_filter_summary.csv"), row.names = FALSE)
write.csv(target_res_all, file.path(table_dir, "12_boundary_target_DESeq2_Aged_vs_Young.csv"), row.names = FALSE)
write.csv(cell_counts_all, file.path(table_dir, "13_boundary_cell_counts.csv"), row.names = FALSE)




