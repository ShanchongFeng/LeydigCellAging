#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

LOCKED_ROOT <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
OUT <- path.expand(Sys.getenv(
  "GSE254315_PSEUDOBULK_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE254315_PSEUDOBULK_ROOT")
))
DIRS <- list(
  rds = file.path(OUT, "rds"),
  tables = file.path(OUT, "tables"),
  metadata = file.path(OUT, "metadata"),
  figures = file.path(OUT, "figures"),
  session = file.path(OUT, "session")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
on.exit(writeLines(capture.output(sessionInfo()), file.path(DIRS$session, "GSE254315_sessionInfo.txt")), add = TRUE)

obj_file <- Sys.glob(file.path(LOCKED_ROOT, "rds/primary_locked/GSE254315/*3_filtered_clustered_human_seurat.rds"))[1]
if (is.na(obj_file) || !file.exists(obj_file)) stop("Missing RDS input for GSE254315: ", obj_file)
obj <- readRDS(obj_file)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"

required_cols <- c("sample_id", "group", "seurat_clusters")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse=", "))

boundaries <- list(
  core_1_6_18 = c('1', '6', '18'),
  ext_1_6_18_7 = c('1', '6', '18', '7')
)
target_genes <- c('HMGCS2', 'CYP11A1', 'CYP17A1', 'HSD3B1', 'HSD3B2', 'LHCGR', 'STAR', 'INSL3', 'NR5A1', 'FOXO3', 'HDAC1')

write.csv(data.frame(role="input_rds", path=obj_file, exists=file.exists(obj_file)),
          file.path(DIRS$metadata, "GSE254315_input_paths.csv"), row.names=FALSE)

sample_meta <- obj@meta.data |>
  as.data.frame() |>
  count(.data[["sample_id"]], .data[["group"]], name = "n_cells_total")
colnames(sample_meta)[1:2] <- c("donor_id", "group")
write.csv(sample_meta, file.path(DIRS$metadata, "GSE254315_sample_metadata.csv"), row.names=FALSE)

run_boundary <- function(boundary_name, clusters_requested) {
  clusters_existing <- levels(Idents(obj))
  clusters_use <- clusters_requested[clusters_requested %in% clusters_existing]
  if (length(clusters_use) == 0) stop("No valid clusters for boundary ", boundary_name)

  sub_obj <- subset(obj, idents = clusters_use)
  sub_obj <- tryCatch(JoinLayers(sub_obj, assay = "RNA"), error = function(e) sub_obj)
  counts <- tryCatch(
    GetAssayData(sub_obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(sub_obj, assay = "RNA", slot = "counts")
  )
  ids <- sub_obj@meta.data[["sample_id"]]
  id_levels <- unique(ids)

  pb_mat <- sapply(id_levels, function(x) Matrix::rowSums(counts[, ids == x, drop = FALSE]))
  pb_mat <- as.matrix(pb_mat)
  storage.mode(pb_mat) <- "integer"

  info <- unique(sub_obj@meta.data[, c("sample_id", "group"), drop = FALSE])
  info <- info[match(colnames(pb_mat), info[["sample_id"]]), , drop = FALSE]
  colnames(info) <- c("donor_id", "group")
  rownames(info) <- info[[ "donor_id" ]]
  info$group <- factor(info$group, levels = c("Young", "Aged"))
  keep_samples <- !is.na(info$group)
  pb_mat <- pb_mat[, keep_samples, drop = FALSE]
  info <- info[keep_samples, , drop = FALSE]

  design <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    donor_id = info[[ "donor_id" ]],
    group = as.character(info$group),
    design_intercept = 1,
    design_Aged = as.integer(info$group == "Aged"),
    n_cells_in_boundary = as.integer(table(factor(ids, levels = info[[ "donor_id" ]]))),
    stringsAsFactors = FALSE
  )

  write.csv(pb_mat, file.path(DIRS$tables, paste0("GSE254315_", boundary_name, "_pseudobulk_counts.csv")), quote=FALSE)
  write.csv(design, file.path(DIRS$metadata, paste0("GSE254315_", boundary_name, "_design_matrix.csv")), row.names=FALSE)

  target_present <- target_genes[target_genes %in% rownames(pb_mat)]
  target_counts <- pb_mat[target_present, , drop = FALSE]
  sparsity <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    gene = rownames(target_counts),
    total_count = rowSums(target_counts),
    samples_ge1 = rowSums(target_counts >= 1),
    samples_ge5 = rowSums(target_counts >= 5),
    samples_ge10 = rowSums(target_counts >= 10),
    stringsAsFactors = FALSE
  )

  keep_genes <- rowSums(pb_mat) > 0
  sparsity$in_DESeq2_filter <- sparsity$gene %in% rownames(pb_mat)[keep_genes]
  write.csv(sparsity, file.path(DIRS$tables, paste0("GSE254315_", boundary_name, "_target_sparsity.csv")), row.names=FALSE)

  target_status <- data.frame(
    gene = target_present,
    boundary = boundary_name,
    clusters = paste(clusters_use, collapse = "+"),
    in_DESeq2_filter = target_present %in% rownames(pb_mat)[keep_genes],
    n_units = ncol(pb_mat),
    n_young = sum(info$group == "Young"),
    n_aged = sum(info$group == "Aged"),
    stringsAsFactors = FALSE
  )

  if (sum(keep_genes) >= 50 && nlevels(droplevels(info$group)) == 2) {
    dds <- DESeqDataSetFromMatrix(countData = pb_mat[keep_genes,,drop=FALSE], colData = info, design = ~ group)
    dds <- DESeq(dds, quiet=TRUE)
    res <- as.data.frame(results(dds, contrast=c("group","Aged","Young")))
    res$gene <- rownames(res)
    res <- res[, c("gene", setdiff(colnames(res), "gene"))]
    write.csv(res, file.path(DIRS$tables, paste0("GSE254315_", boundary_name, "_DESeq2_Aged_vs_Young.csv")), row.names=FALSE)
    target_res <- target_status |>
      left_join(res, by = "gene") |>
      mutate(direction = case_when(is.na(log2FoldChange) ~ "not_tested",
                                   log2FoldChange < 0 ~ "Aged_lower",
                                   log2FoldChange > 0 ~ "Aged_higher",
                                   TRUE ~ "zero"))
  } else {
    target_res <- target_status
    target_res$baseMean <- NA_real_; target_res$log2FoldChange <- NA_real_; target_res$lfcSE <- NA_real_
    target_res$stat <- NA_real_; target_res$pvalue <- NA_real_; target_res$padj <- NA_real_; target_res$direction <- "not_tested"
  }
  write.csv(target_res, file.path(DIRS$tables, paste0("GSE254315_", boundary_name, "_target_genes_DESeq2.csv")), row.names=FALSE)
  target_res
}

all_targets <- bind_rows(lapply(names(boundaries), function(nm) run_boundary(nm, boundaries[[nm]])))
write.csv(all_targets, file.path(DIRS$tables, "GSE254315_all_boundaries_target_genes_DESeq2.csv"), row.names=FALSE)
