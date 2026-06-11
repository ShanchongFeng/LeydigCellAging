suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

base_dir <- path.expand(Sys.getenv(
  "BULK_GSE287203_ROOT",
  unset = file.path(Sys.getenv("HOME"), "BULK_GSE287203_ROOT")
))
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
setwd(base_dir)

dir.create("results/in_vitro", recursive = TRUE, showWarnings = FALSE)
dir.create("results/in_vivo", recursive = TRUE, showWarnings = FALSE)
dir.create("results/shared_qc", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

counts <- read.table(
  "data/GSE287203_raw_counts.txt.gz",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(counts)[1] <- "gene"
rownames(counts) <- counts$gene
counts <- counts[, -1, drop = FALSE]

meta <- read.csv("data/metadata_clean.csv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

sample_order <- meta$sample_id
stopifnot(all(sample_order %in% colnames(counts)))
counts <- counts[, sample_order, drop = FALSE]

write.table(
  data.frame(
    n_genes_raw = nrow(counts),
    n_samples_raw = ncol(counts)
  ),
  file = "results/shared_qc/raw_counts_dimensions.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

run_deseq_pipeline <- function(model_type, group_levels, contrast_vec, prefix) {

  meta_sub <- meta[meta$model_type == model_type, , drop = FALSE]
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]

  meta_sub$group <- factor(meta_sub$group, levels = group_levels)

  min_group_n <- min(table(meta_sub$group))
  keep <- rowSums(counts_sub >= 10) >= min_group_n

  counts_filt <- counts_sub[keep, , drop = FALSE]

  write.csv(
    counts_filt,
    file = file.path("results", prefix, paste0("counts_filtered_", prefix, ".csv")),
    quote = FALSE
  )

  qc_lib <- data.frame(
    sample_id = colnames(counts_sub),
    library_size_raw = colSums(counts_sub),
    library_size_filtered = colSums(counts_filt)
  )
  write.csv(
    qc_lib,
    file = file.path("results", prefix, paste0("library_sizes_", prefix, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts_filt)),
    colData = meta_sub,
    design = ~ group
  )

  dds <- DESeq(dds)

  res <- results(dds, contrast = contrast_vec)
  res <- res[order(res$padj), ]
  res_df <- data.frame(Gene = rownames(res), as.data.frame(res), row.names = NULL)

  coef_name <- resultsNames(dds)[2]
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  res_shrunk <- res_shrunk[order(res_shrunk$padj), ]
  res_shrunk_df <- data.frame(
    Gene = rownames(res_shrunk),
    as.data.frame(res_shrunk),
    row.names = NULL
  )

  write.csv(
    res_df,
    file = file.path("results", prefix, paste0("DEG_", prefix, "_raw.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  write.csv(
    res_shrunk_df,
    file = file.path("results", prefix, paste0("DEG_", prefix, "_shrunk.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  enrichment_name <- switch(
    prefix,
    in_vitro = "GSE287203_invitro_H_vs_C_DESeq2_results.csv",
    in_vivo = "GSE287203_invivo_Old_vs_Young_DESeq2_results.csv",
    stop("Unsupported bulk prefix: ", prefix)
  )
  write.csv(
    res_df,
    file = file.path("tables", enrichment_name),
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(dds, file = file.path("results", prefix, paste0("dds_", prefix, ".rds")))

  vsd <- vst(dds, blind = TRUE)
  saveRDS(vsd, file = file.path("results", prefix, paste0("vsd_", prefix, ".rds")))

  vst_mat <- assay(vsd)
  write.csv(
    vst_mat,
    file = file.path("results", prefix, paste0("vst_", prefix, ".csv")),
    quote = FALSE
  )

  pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  p <- ggplot(pcaData, aes(PC1, PC2, color = group, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3.5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste0("PCA - ", prefix)) +
    theme_bw(base_size = 12)

  ggsave(
    filename = file.path("results", prefix, paste0("PCA_", prefix, ".pdf")),
    plot = p,
    width = 7,
    height = 5
  )

  sample_dists <- dist(t(assay(vsd)))
  sample_dist_mat <- as.matrix(sample_dists)
  rownames(sample_dist_mat) <- colnames(vsd)
  colnames(sample_dist_mat) <- colnames(vsd)

  pdf(file.path("results", prefix, paste0("sample_distance_heatmap_", prefix, ".pdf")), width = 6, height = 5)
  pheatmap(sample_dist_mat, main = paste0("Sample distance - ", prefix))
  dev.off()

  pdf(file.path("results", prefix, paste0("MAplot_", prefix, ".pdf")), width = 6, height = 5)
  plotMA(res_shrunk, main = paste0("MA plot - ", prefix), ylim = c(-5, 5))
  dev.off()

  summary_tab <- data.frame(
    prefix = prefix,
    n_samples = ncol(counts_sub),
    n_genes_raw = nrow(counts_sub),
    n_genes_filtered = nrow(counts_filt),
    n_sig_padj_0.05 = sum(res$padj < 0.05, na.rm = TRUE),
    n_sig_shrunk_padj_0.05 = sum(res_shrunk$padj < 0.05, na.rm = TRUE)
  )

  write.csv(
    summary_tab,
    file = file.path("results", prefix, paste0("summary_", prefix, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
}

run_deseq_pipeline(
  model_type = "in_vitro",
  group_levels = c("C", "H"),
  contrast_vec = c("group", "H", "C"),
  prefix = "in_vitro"
)

run_deseq_pipeline(
  model_type = "in_vivo",
  group_levels = c("young", "old"),
  contrast_vec = c("group", "old", "young"),
  prefix = "in_vivo"
)

