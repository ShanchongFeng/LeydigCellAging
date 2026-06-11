#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

locked_root <- normalizePath(path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
)), mustWork = TRUE)
out_root <- path.expand(Sys.getenv(
  "GSE182786_BOUNDARY_STRESS_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE182786_BOUNDARY_STRESS_ROOT")
))
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), recursive = TRUE, showWarnings = FALSE)

n_iter <- 2000
max_cells_per_sample <- 200
pseudocounts <- c(0.1, 0.5, 1)
set.seed(20260601)

get_counts <- function(obj) {
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = DefaultAssay(obj)), error = function(e) obj)
  counts <- tryCatch(
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
    error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  )
  list(obj = obj, counts = as(counts, "dgCMatrix"))
}

exact_permutation_p <- function(x, group, reference) {
  group <- droplevels(factor(group))
  if (length(levels(group)) != 2 || min(table(group)) < 2) return(NA_real_)
  group <- relevel(group, ref = reference)
  n_ref <- sum(group == levels(group)[1])
  observed <- mean(x[group == levels(group)[2]]) - mean(x[group == levels(group)[1]])
  ref_sets <- combn(seq_along(x), n_ref, simplify = FALSE)
  null_effects <- vapply(ref_sets, function(ref_idx) {
    alt_idx <- setdiff(seq_along(x), ref_idx)
    mean(x[alt_idx]) - mean(x[ref_idx])
  }, numeric(1))
  mean(abs(null_effects) >= abs(observed) - 1e-12)
}

summarize_iterations <- function(x) {
  split_x <- split(x, paste(x$boundary, x$comparison, x$pseudocount, sep = "|||"))
  rows <- lapply(split_x, function(df) {
    data.frame(
      boundary = df$boundary[[1]],
      comparison = df$comparison[[1]],
      pseudocount = df$pseudocount[[1]],
      sampled_cells_per_sample = df$sampled_cells_per_sample[[1]],
      n_iterations = nrow(df),
      median_effect = median(df$effect_alternative_minus_reference),
      mean_effect = mean(df$effect_alternative_minus_reference),
      ci95_low = unname(quantile(df$effect_alternative_minus_reference, 0.025)),
      ci95_high = unname(quantile(df$effect_alternative_minus_reference, 0.975)),
      pct_negative = mean(df$effect_alternative_minus_reference < 0),
      pct_exact_permutation_p_lt_0.05 = mean(df$exact_permutation_p_two_sided < 0.05),
      median_exact_permutation_p_two_sided = median(df$exact_permutation_p_two_sided)
    )
  })
  do.call(rbind, rows)
}

object_path <- file.path(locked_root, "rds", "primary_locked", "GSE182786_check", "GSE182786_check__results_01_basic_clustered.rds")
metadata_path <- file.path(locked_root, "tables", "source_copies", "hmgcs2_gse182786_phase1", "data", "meta", "gse182786_sample_metadata_curated.csv")

obj <- readRDS(object_path)
got <- get_counts(obj)
obj <- got$obj
counts <- got$counts
meta <- obj@meta.data
meta$seurat_clusters <- as.character(meta$seurat_clusters)
meta$sample_id <- as.character(meta$sample_id)
cell_library_sizes <- Matrix::colSums(counts)
hmgcs2_cell_counts <- as.numeric(counts["HMGCS2", ])
names(hmgcs2_cell_counts) <- colnames(counts)

sample_meta <- read.csv(metadata_path, check.names = FALSE)
sample_meta$sample_id <- sample_meta$cell_prefix
sample_meta <- sample_meta[, c("sample_id", "bmi_group_curated")]

boundaries <- list(CMB_0_17 = c("0", "17"), ESB_0_14_17 = c("0", "14", "17"))
comparisons <- list(
  all_older_vs_young = list(include = c("reference_young", "older_normal_BMI", "older_high_BMI"), reference = "reference_young", alternative = "older_all"),
  older_normal_vs_young = list(include = c("reference_young", "older_normal_BMI"), reference = "reference_young", alternative = "older_normal_BMI"),
  older_high_vs_young = list(include = c("reference_young", "older_high_BMI"), reference = "reference_young", alternative = "older_high_BMI"),
  older_high_vs_older_normal = list(include = c("older_normal_BMI", "older_high_BMI"), reference = "older_normal_BMI", alternative = "older_high_BMI")
)

make_groups <- function(sample_ids, comparison) {
  strata <- sample_meta$bmi_group_curated[match(sample_ids, sample_meta$sample_id)]
  keep <- strata %in% comparison$include
  groups <- strata[keep]
  if (comparison$alternative == "older_all") {
    groups[groups %in% c("older_normal_BMI", "older_high_BMI")] <- "older_all"
  }
  list(
    keep = keep,
    group = factor(groups, levels = c(comparison$reference, comparison$alternative))
  )
}

all_rows <- list()
design_rows <- list()

for (boundary_name in names(boundaries)) {
  clusters <- boundaries[[boundary_name]]
  cells <- rownames(meta)[meta$seurat_clusters %in% clusters & meta$sample_id %in% sample_meta$sample_id]
  boundary_meta <- meta[cells, , drop = FALSE]
  sample_to_cells <- split(rownames(boundary_meta), boundary_meta$sample_id)
  sample_ids <- sample_meta$sample_id
  sample_counts <- vapply(sample_ids, function(s) length(sample_to_cells[[s]]), integer(1))
  sample_n <- min(max_cells_per_sample, min(sample_counts))
  design_rows[[boundary_name]] <- data.frame(
    boundary = boundary_name,
    clusters = paste(clusters, collapse = "+"),
    n_samples = length(sample_ids),
    min_observed_cells_per_sample = min(sample_counts),
    max_observed_cells_per_sample = max(sample_counts),
    sampled_cells_per_sample = sample_n,
    n_iterations = n_iter
  )

  for (iter in seq_len(n_iter)) {
    sample_metrics <- lapply(sample_ids, function(sample_id) {
      selected <- sample(sample_to_cells[[sample_id]], sample_n, replace = FALSE)
      c(
        HMGCS2 = sum(hmgcs2_cell_counts[selected]),
        library_size = sum(cell_library_sizes[selected])
      )
    })
    sample_metrics <- do.call(cbind, sample_metrics)
    colnames(sample_metrics) <- sample_ids

    for (pc in pseudocounts) {
      hmgcs2_logcpm <- log2(((sample_metrics["HMGCS2", ] + pc) / (sample_metrics["library_size", ] + 2 * pc)) * 1e6)
      for (comparison_name in names(comparisons)) {
        comparison <- comparisons[[comparison_name]]
        grouped <- make_groups(sample_ids, comparison)
        x <- as.numeric(hmgcs2_logcpm[grouped$keep])
        group <- grouped$group
        effect <- mean(x[group == comparison$alternative]) - mean(x[group == comparison$reference])
        all_rows[[length(all_rows) + 1]] <- data.frame(
          boundary = boundary_name,
          iteration = iter,
          pseudocount = pc,
          comparison = comparison_name,
          sampled_cells_per_sample = sample_n,
          effect_alternative_minus_reference = effect,
          exact_permutation_p_two_sided = exact_permutation_p(x, group, comparison$reference)
        )
      }
    }
    if (iter %% 100 == 0) message(boundary_name, " iteration ", iter, "/", n_iter)
  }
}

iterations <- do.call(rbind, all_rows)
summary <- summarize_iterations(iterations)
design <- do.call(rbind, design_rows)

write.csv(design, file.path(out_root, "tables", "13_BMI_stratified_downsampling_design.csv"), row.names = FALSE)
write.csv(iterations, file.path(out_root, "tables", "14_BMI_stratified_downsampling_iterations.csv"), row.names = FALSE)
write.csv(summary, file.path(out_root, "tables", "15_BMI_stratified_downsampling_summary.csv"), row.names = FALSE)

plot_df <- summary[summary$pseudocount == 0.5 & summary$comparison != "older_high_vs_older_normal", ]
plot_df$comparison <- factor(
  plot_df$comparison,
  levels = c("all_older_vs_young", "older_normal_vs_young", "older_high_vs_young"),
  labels = c("All older vs young", "Older normal BMI vs young", "Older high BMI vs young")
)
p <- ggplot(plot_df, aes(x = median_effect, y = comparison, colour = boundary)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci95_low, xmax = ci95_high), height = 0.12, linewidth = 0.5, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.35)) +
  scale_colour_manual(
    values = c(CMB_0_17 = "#1F5A75", ESB_0_14_17 = "#B85B3F"),
    labels = c(CMB_0_17 = "CMB (0+17)", ESB_0_14_17 = "ESB (0+14+17)")
  ) +
  labs(x = "Median HMGCS2 logCPM effect after equal-cell recovery", y = NULL, colour = "Boundary") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

grDevices::cairo_pdf(file.path(out_root, "figures", "Fig_S4_HMGCS2_BMI_stratified_equal_cell_downsampling.pdf"), width = 170 / 25.4, height = 72 / 25.4, family = "Arial")
print(p)
grDevices::dev.off()
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(file.path(out_root, "figures", "Fig_S4_HMGCS2_BMI_stratified_equal_cell_downsampling.png"), width = 170 / 25.4, height = 72 / 25.4, units = "in", res = 600)
} else {
  grDevices::png(file.path(out_root, "figures", "Fig_S4_HMGCS2_BMI_stratified_equal_cell_downsampling.png"), width = 170 / 25.4, height = 72 / 25.4, units = "in", res = 600, type = "cairo")
}
print(p)
grDevices::dev.off()

capture.output(sessionInfo(), file = file.path(out_root, "logs", "4_bmi_sessionInfo.txt"))
