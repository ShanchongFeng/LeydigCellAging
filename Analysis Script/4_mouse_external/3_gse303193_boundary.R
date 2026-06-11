#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

base_dir <- path.expand(Sys.getenv(
  "GSE303193_EXTERNAL_MOUSE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE303193_EXTERNAL_MOUSE_ROOT")
))
table_dir <- file.path(base_dir, "tables")
rds_dir <- file.path(base_dir, "rds")

obj_file <- file.path(rds_dir, "3_clustered_seurat_res0.4.rds")
marker_rank_file <- file.path(table_dir, "7_leydig_marker_score_by_cluster.csv")
module_file <- path.expand(Sys.getenv(
  "HDWGCNA_MODULE_ASSIGNMENT_CSV",
  unset = file.path(
    path.expand(Sys.getenv("HDWGCNA_ROOT", unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT"))),
    "tables", "hdwgcna_module_assignment.csv"
  )
))

if (!file.exists(obj_file)) stop("Missing object: ", obj_file)
if (!file.exists(marker_rank_file)) stop("Missing marker rank file: ", marker_rank_file)
if (!file.exists(module_file)) stop("Missing module file: ", module_file)


obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
Idents(obj) <- "seurat_clusters"

marker_rank <- read.csv(marker_rank_file, check.names = FALSE)
marker_rank$cluster_clean <- sub("^g", "", marker_rank$cluster)

mods <- read.csv(module_file, check.names = FALSE)

blue_genes <- unique(mods$gene_name[mods$module == "blue"])
brown_genes <- unique(mods$gene_name[mods$module == "brown"])
turquoise_genes <- unique(mods$gene_name[mods$module == "turquoise"])

target_genes <- c(
  "Hmgcs2", "Cyp11a1", "Cyp17a1", "Star", "Hsd3b1",
  "Lhcgr", "Nr5a1", "Foxo3", "Insl3"
)

old_boundary <- c("16", "17", "19")
auto_top3 <- marker_rank$cluster_clean[1:3]
union_boundary <- unique(c(old_boundary, auto_top3))

boundaries <- list(
  old_boundary_16_17_19 = old_boundary,
  auto_top3_marker = auto_top3,
  union_old_auto = union_boundary
)


run_one_boundary <- function(boundary_name, clusters_use) {

  clusters_existing <- levels(Idents(obj))
  clusters_use <- clusters_use[clusters_use %in% clusters_existing]

  if (length(clusters_use) == 0) {
    stop("No valid clusters for boundary: ", boundary_name)
  }

  sub_obj <- subset(obj, idents = clusters_use)
  DefaultAssay(sub_obj) <- "RNA"


  # module overlaps
  blue_use <- intersect(blue_genes, rownames(sub_obj))
  brown_use <- intersect(brown_genes, rownames(sub_obj))
  turquoise_use <- intersect(turquoise_genes, rownames(sub_obj))

  overlap_df <- data.frame(
    boundary = boundary_name,
    module = c("blue", "brown", "turquoise"),
    module_n = c(length(blue_genes), length(brown_genes), length(turquoise_genes)),
    overlap_n = c(length(blue_use), length(brown_use), length(turquoise_use)),
    stringsAsFactors = FALSE
  )

  # AddModuleScore
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
  ) |>
    dplyr::group_by(sample_id, group) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(module_score_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  score_by_sample$boundary <- boundary_name
  score_by_sample$clusters <- paste(clusters_use, collapse = "+")

  # Pseudobulk DESeq2
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

  keep_genes <- rowSums(pb_mat >= 10) >= 3

  dds <- DESeqDataSetFromMatrix(
    countData = pb_mat[keep_genes, ],
    colData = sample_info,
    design = ~ group
  )

  dds <- DESeq(dds, quiet = TRUE)

  res <- results(dds, contrast = c("group", "Aged", "Young"))
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)

  target_res <- res_df |>
    dplyr::filter(Gene %in% target_genes) |>
    dplyr::arrange(match(Gene, target_genes))

  target_res$boundary <- boundary_name
  target_res$clusters <- paste(clusters_use, collapse = "+")
  target_res$n_cells <- ncol(sub_obj)

  # cell counts
  cell_counts <- as.data.frame(table(sub_obj$sample_id, sub_obj$group, sub_obj$seurat_clusters))
  colnames(cell_counts) <- c("sample_id", "group", "cluster", "n_cells")
  cell_counts$boundary <- boundary_name
  cell_counts$clusters <- paste(clusters_use, collapse = "+")

  saveRDS(
    sub_obj,
    file.path(rds_dir, paste0("6_boundary_sensitivity_", boundary_name, ".rds"))
  )

  list(
    overlap = overlap_df,
    score_by_sample = score_by_sample,
    target_res = target_res,
    cell_counts = cell_counts
  )
}

res_list <- list()

for (nm in names(boundaries)) {
  res_list[[nm]] <- run_one_boundary(nm, boundaries[[nm]])
}

overlap_all <- bind_rows(lapply(res_list, `[[`, "overlap"))
score_all <- bind_rows(lapply(res_list, `[[`, "score_by_sample"))
target_all <- bind_rows(lapply(res_list, `[[`, "target_res"))
cell_counts_all <- bind_rows(lapply(res_list, `[[`, "cell_counts"))

write.csv(
  overlap_all,
  file.path(table_dir, "15_boundary_sensitivity_module_overlap.csv"),
  row.names = FALSE
)

write.csv(
  score_all,
  file.path(table_dir, "16_boundary_sensitivity_module_score_by_sample.csv"),
  row.names = FALSE
)

write.csv(
  target_all,
  file.path(table_dir, "17_boundary_sensitivity_target_gene_DESeq2.csv"),
  row.names = FALSE
)

write.csv(
  cell_counts_all,
  file.path(table_dir, "18_boundary_sensitivity_cell_counts.csv"),
  row.names = FALSE
)

# Compact interpretation table
blue_score_summary <- score_all |>
  dplyr::group_by(boundary, group) |>
  dplyr::summarise(
    blue_mean = mean(blue_module1, na.rm = TRUE),
    blue_median = median(blue_module1, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  blue_score_summary,
  file.path(table_dir, "19_boundary_sensitivity_blue_score_group_summary.csv"),
  row.names = FALSE
)

target_compact <- target_all |>
  dplyr::filter(Gene %in% c("Hmgcs2", "Cyp11a1", "Star", "Lhcgr", "Insl3", "Nr5a1")) |>
  dplyr::select(boundary, clusters, Gene, log2FoldChange, pvalue, padj, n_cells)

write.csv(
  target_compact,
  file.path(table_dir, "20_boundary_sensitivity_target_compact.csv"),
  row.names = FALSE
)



