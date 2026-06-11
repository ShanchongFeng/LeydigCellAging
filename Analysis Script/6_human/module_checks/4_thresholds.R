#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv("GSE182786_SCHEME3_ROOT", unset = file.path(Sys.getenv("HOME"), "GSE182786_SCHEME3_ROOT")))
locked_root <- path.expand(Sys.getenv("LOCKED_INPUT_ROOT", unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")))
dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

core_rds <- file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__21_human_leydig_core_0_17_with_module_scores.rds")

old_threshold_core_genes <- file.path(locked_root, "tables/primary_locked/GSE182786_check/GSE182786_check__43_threshold_sensitivity_core_genes.csv")
old_threshold_summary    <- file.path(locked_root, "tables/primary_locked/GSE182786_check/GSE182786_check__44_threshold_sensitivity_summary.csv")

if (!file.exists(path.expand(core_rds))) {
  stop("找不到输入文件：", core_rds)
}


core_obj <- readRDS(path.expand(core_rds))
DefaultAssay(core_obj) <- "RNA"


pb_counts_core <- AggregateExpression(
  core_obj,
  group.by = "sample_id",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

pb_mat <- as.matrix(pb_counts_core$RNA)

sample_info <- unique(core_obj@meta.data[, c("sample_id", "group")])
sample_info <- sample_info[match(colnames(pb_mat), sample_info$sample_id), , drop = FALSE]
rownames(sample_info) <- sample_info$sample_id
sample_info$group <- factor(sample_info$group, levels = c("Young", "Aged"))



target_genes <- c(
  "HMGCS2", "CYP11A1", "NR5A1", "INSL3",
  "STAR", "HSD3B2", "LHCGR", "CYP17A1"
)

target_present <- intersect(target_genes, rownames(pb_mat))
target_missing <- setdiff(target_genes, rownames(pb_mat))



pb_target_counts <- pb_mat[target_present, , drop = FALSE]
write_csv(
  as.data.frame(pb_target_counts) %>% rownames_to_column("Gene"),
  file.path(base_dir, "tables", "17_core_gene_raw_pseudobulk_counts.csv")
)

gene_sparsity <- tibble(
  Gene = target_present,
  mean_count = apply(pb_target_counts, 1, mean),
  median_count = apply(pb_target_counts, 1, median),
  min_count = apply(pb_target_counts, 1, min),
  max_count = apply(pb_target_counts, 1, max),
  n_zero = apply(pb_target_counts, 1, function(x) sum(x == 0)),
  prop_zero = apply(pb_target_counts, 1, function(x) mean(x == 0)),
  n_lt5 = apply(pb_target_counts, 1, function(x) sum(x < 5)),
  prop_lt5 = apply(pb_target_counts, 1, function(x) mean(x < 5)),
  n_lt10 = apply(pb_target_counts, 1, function(x) sum(x < 10)),
  prop_lt10 = apply(pb_target_counts, 1, function(x) mean(x < 10))
) %>%
  arrange(desc(prop_zero), mean_count)

write_csv(gene_sparsity, file.path(base_dir, "tables", "18_core_gene_sparsity_metrics.csv"))

threshold_grid <- tribble(
  ~threshold_label, ~min_count, ~min_samples,
  "count_ge10_in_ge3", 10, 3,
  "count_ge5_in_ge3",   5, 3,
  "count_ge10_in_ge2", 10, 2
)

write_csv(threshold_grid, file.path(base_dir, "tables", "19_threshold_grid_used.csv"))

run_deseq_threshold <- function(count_mat, coldata, min_count = 10, min_samples = 3) {
  keep <- rowSums(count_mat >= min_count) >= min_samples

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat[keep, , drop = FALSE],
    colData = coldata,
    design = ~ group
  )

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "Aged", "Young"))
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)

  list(
    keep = keep,
    res = res_df
  )
}

threshold_gene_long <- map_dfr(seq_len(nrow(threshold_grid)), function(i) {
  one <- threshold_grid[i, ]

  out <- run_deseq_threshold(
    count_mat = pb_mat,
    coldata = sample_info,
    min_count = one$min_count,
    min_samples = one$min_samples
  )

  keep_vec <- out$keep
  kept_genes <- names(keep_vec)[keep_vec]

  res_df <- out$res

  map_dfr(target_genes, function(g) {
    pass_filter <- g %in% kept_genes

    if (g %in% res_df$Gene) {
      row <- res_df %>% filter(Gene == g) %>% slice(1)
      tibble(
        threshold_label = one$threshold_label,
        min_count = one$min_count,
        min_samples = one$min_samples,
        Gene = g,
        present_in_matrix = g %in% rownames(pb_mat),
        pass_filter = pass_filter,
        baseMean = row$baseMean,
        log2FoldChange = row$log2FoldChange,
        lfcSE = row$lfcSE,
        stat = row$stat,
        pvalue = row$pvalue,
        padj = row$padj
      )
    } else {
      tibble(
        threshold_label = one$threshold_label,
        min_count = one$min_count,
        min_samples = one$min_samples,
        Gene = g,
        present_in_matrix = g %in% rownames(pb_mat),
        pass_filter = pass_filter,
        baseMean = NA_real_,
        log2FoldChange = NA_real_,
        lfcSE = NA_real_,
        stat = NA_real_,
        pvalue = NA_real_,
        padj = NA_real_
      )
    }
  })
})

threshold_gene_long <- threshold_gene_long %>%
  mutate(
    direction = case_when(
      is.na(log2FoldChange) ~ "not_tested",
      log2FoldChange > 0 ~ "positive",
      log2FoldChange < 0 ~ "negative",
      TRUE ~ "zero"
    )
  )

write_csv(threshold_gene_long, file.path(base_dir, "tables", "20_threshold_sensitivity_all_thresholds_long.csv"))

gene_stability_summary <- threshold_gene_long %>%
  group_by(Gene) %>%
  summarise(
    present_in_matrix = any(present_in_matrix),
    n_thresholds = n(),
    n_passed = sum(pass_filter, na.rm = TRUE),
    any_passed = any(pass_filter, na.rm = TRUE),
    n_tested_with_lfc = sum(!is.na(log2FoldChange)),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    min_log2FC = suppressWarnings(min(log2FoldChange, na.rm = TRUE)),
    max_log2FC = suppressWarnings(max(log2FoldChange, na.rm = TRUE)),
    range_log2FC = max_log2FC - min_log2FC,
    any_positive = any(log2FoldChange > 0, na.rm = TRUE),
    any_negative = any(log2FoldChange < 0, na.rm = TRUE),
    cross_zero_across_thresholds = any_positive & any_negative,
    direction_set = paste(sort(unique(direction[direction %in% c("positive", "negative", "zero")])), collapse = ";"),
    direction_consistent = !cross_zero_across_thresholds,
    best_pvalue = suppressWarnings(min(pvalue, na.rm = TRUE)),
    best_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(gene_sparsity, by = "Gene") %>%
  arrange(desc(prop_zero), range_log2FC)

gene_stability_summary$best_pvalue[is.infinite(gene_stability_summary$best_pvalue)] <- NA_real_
gene_stability_summary$best_padj[is.infinite(gene_stability_summary$best_padj)] <- NA_real_
gene_stability_summary$min_log2FC[is.infinite(gene_stability_summary$min_log2FC)] <- NA_real_
gene_stability_summary$max_log2FC[is.infinite(gene_stability_summary$max_log2FC)] <- NA_real_

write_csv(gene_stability_summary, file.path(base_dir, "tables", "21_threshold_gene_stability_summary.csv"))

if (file.exists(path.expand(old_threshold_core_genes))) {
  old_core_genes <- suppressWarnings(read_csv(path.expand(old_threshold_core_genes), show_col_types = FALSE))
  write_csv(old_core_genes, file.path(base_dir, "tables", "old_43_threshold_sensitivity_core_genes_backup.csv"))
}

if (file.exists(path.expand(old_threshold_summary))) {
  old_summary <- suppressWarnings(read_csv(path.expand(old_threshold_summary), show_col_types = FALSE))
  write_csv(old_summary, file.path(base_dir, "tables", "old_44_threshold_sensitivity_summary_backup.csv"))
}

plot_df1 <- threshold_gene_long %>%
  filter(Gene %in% target_present) %>%
  mutate(
    threshold_label = factor(
      threshold_label,
      levels = c("count_ge10_in_ge3", "count_ge5_in_ge3", "count_ge10_in_ge2")
    )
  )

p1 <- ggplot(plot_df1, aes(x = threshold_label, y = log2FoldChange, group = 1)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_line(na.rm = TRUE) +
  geom_point(aes(shape = pass_filter), size = 2.8, na.rm = TRUE) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  theme_bw(base_size = 12) +
  labs(
    title = "Threshold sensitivity of donor-level single-gene log2FC",
    x = NULL,
    y = "log2FC (Aged vs Young)"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "7_threshold_sensitivity_lfc_by_gene.png"),
  plot = p1,
  width = 12,
  height = 8,
  dpi = 300
)

sample_info_tbl <- sample_info %>%
  as.data.frame() %>%
  tibble::as_tibble()

count_heat_df <- as.data.frame(pb_target_counts) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample_id",
    values_to = "raw_count"
  ) %>%
  left_join(
    sample_info_tbl,
    by = "sample_id"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = rev(target_present)),
    sample_id = factor(sample_id, levels = sample_info$sample_id),
    log10_count = log10(raw_count + 1)
  )

p2 <- ggplot(count_heat_df, aes(x = sample_id, y = Gene, fill = log10_count)) +
  geom_tile(color = "white") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Raw donor-level pseudobulk counts of key genes",
    x = NULL,
    y = NULL,
    fill = "log10(count+1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "8_core_gene_raw_pseudobulk_heatmap.png"),
  plot = p2,
  width = 10,
  height = 5.5,
  dpi = 300
)

plot_df3 <- gene_sparsity %>%
  mutate(Gene = factor(Gene, levels = Gene[order(prop_zero, decreasing = TRUE)]))

p3 <- ggplot(plot_df3, aes(x = Gene, y = prop_zero)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Proportion of zero counts across donor-level pseudobulk samples",
    x = NULL,
    y = "Proportion zero"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "9_core_gene_zero_proportion.png"),
  plot = p3,
  width = 7,
  height = 4.8,
  dpi = 300
)





