#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
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
dirs <- file.path(out_root, c("scripts", "tables", "figures", "docs", "logs"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, name) {
  path <- file.path(out_root, "tables", name)
  write.csv(x, path, row.names = FALSE, na = "")
  path
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 4, format = "f")))
}

get_counts <- function(obj) {
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = DefaultAssay(obj)), error = function(e) obj)
  counts <- tryCatch(
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
    error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  )
  list(obj = obj, counts = counts)
}

wilcox_safe <- function(x, group) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(levels(group)) != 2 || min(table(group)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x ~ group, exact = FALSE)$p.value)
}

lm_safe <- function(x, group, reference) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(levels(group)) != 2 || min(table(group)) < 2) return(c(beta = NA_real_, p = NA_real_))
  group <- stats::relevel(group, ref = reference)
  fit <- summary(lm(x ~ group))
  c(beta = unname(coef(fit)[2, 1]), p = unname(coef(fit)[2, 4]))
}

exact_permutation <- function(x, group, reference) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(levels(group)) != 2 || min(table(group)) < 2) {
    return(c(effect = NA_real_, p_two_sided = NA_real_, p_lower_one_sided = NA_real_, n_permutations = NA_real_))
  }
  group <- stats::relevel(group, ref = reference)
  ref <- levels(group)[1]
  alt <- levels(group)[2]
  n_ref <- sum(group == ref)
  observed <- mean(x[group == alt]) - mean(x[group == ref])
  ref_sets <- combn(seq_along(x), n_ref, simplify = FALSE)
  null_effects <- vapply(ref_sets, function(ref_idx) {
    alt_idx <- setdiff(seq_along(x), ref_idx)
    mean(x[alt_idx]) - mean(x[ref_idx])
  }, numeric(1))
  c(
    effect = observed,
    p_two_sided = mean(abs(null_effects) >= abs(observed) - 1e-12),
    p_lower_one_sided = mean(null_effects <= observed + 1e-12),
    n_permutations = length(null_effects)
  )
}

bootstrap_ci <- function(x, group, reference, n_boot = 10000, seed = 20260601) {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(levels(group)) != 2 || min(table(group)) < 2) return(c(low = NA_real_, high = NA_real_))
  group <- stats::relevel(group, ref = reference)
  ref <- levels(group)[1]
  alt <- levels(group)[2]
  ref_x <- x[group == ref]
  alt_x <- x[group == alt]
  set.seed(seed)
  effects <- replicate(n_boot, mean(sample(alt_x, length(alt_x), replace = TRUE)) - mean(sample(ref_x, length(ref_x), replace = TRUE)))
  unname(stats::quantile(effects, c(0.025, 0.975), na.rm = TRUE))
}

save_plot <- function(p, stem, width_mm = 175, height_mm = 90, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(p)
  grDevices::dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
  } else {
    grDevices::png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi, type = "cairo")
  }
  print(p)
  grDevices::dev.off()
}

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

object_path <- file.path(locked_root, "rds", "primary_locked", "GSE182786_check", "GSE182786_check__results_01_basic_clustered.rds")
metadata_path <- file.path(locked_root, "tables", "source_copies", "hmgcs2_gse182786_phase1", "data", "meta", "gse182786_sample_metadata_curated.csv")

obj <- readRDS(object_path)
count_result <- get_counts(obj)
obj <- count_result$obj
counts <- count_result$counts
meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)
meta$seurat_clusters <- as.character(meta$seurat_clusters)
meta$sample_id <- as.character(meta$sample_id)
meta$group <- as.character(meta$group)

sample_meta <- read.csv(metadata_path, check.names = FALSE)
sample_meta$sample_id <- sample_meta$cell_prefix
sample_meta$group <- ifelse(sample_meta$age_group == "young", "Young", "Aged")
sample_meta <- sample_meta[, c("sample_id", "donor_id", "geo_accession", "group", "analysis_group", "bmi_group_curated")]

boundaries <- list(CMB_0_17 = c("0", "17"), ESB_0_14_17 = c("0", "14", "17"))
target_genes <- c("HMGCS2", "CYP11A1", "NR5A1", "INSL3", "STAR", "HSD3B2", "LHCGR", "CYP17A1", "FOXO3")
program_sets <- list(
  steroidogenesis = c("STAR", "CYP11A1", "HSD3B2", "HSD3B1", "CYP17A1", "HSD17B3", "LHCGR", "NR5A1", "SCARB1"),
  ketogenesis_fao = c("HMGCS2", "HMGCL", "BDH1", "OXCT1", "CPT1A", "CPT2", "ACADL", "ACADM", "HADHA", "HADH", "ACOX1"),
  cholesterol_handling = c("HMGCR", "HMGCS1", "SQLE", "LSS", "FDFT1", "DHCR7", "SCARB1", "SOAT1"),
  detox_antioxidant = c("GSTA1", "GSTA2", "GSTM1", "GSTM2", "NQO1", "GPX1", "GPX3", "PRDX1", "SOD2")
)

aggregate_boundary <- function(boundary_name, clusters) {
  cells <- rownames(meta)[meta$seurat_clusters %in% clusters & meta$sample_id %in% sample_meta$sample_id]
  sub_meta <- meta[cells, , drop = FALSE]
  sub_counts <- counts[, cells, drop = FALSE]
  samples <- sample_meta$sample_id
  design_matrix <- sparseMatrix(
    i = seq_along(cells),
    j = match(sub_meta$sample_id, samples),
    x = 1,
    dims = c(length(cells), length(samples)),
    dimnames = list(cells, samples)
  )
  pseudobulk <- sub_counts %*% design_matrix
  storage.mode(pseudobulk@x) <- "numeric"
  cell_counts <- as.data.frame(table(sample_id = sub_meta$sample_id), stringsAsFactors = FALSE)
  names(cell_counts)[2] <- "n_cells"
  merged_meta <- merge(sample_meta, cell_counts, by = "sample_id", all.x = TRUE)
  merged_meta$n_cells[is.na(merged_meta$n_cells)] <- 0
  merged_meta$boundary <- boundary_name
  merged_meta$clusters <- paste(clusters, collapse = "+")
  list(counts = pseudobulk[, merged_meta$sample_id, drop = FALSE], meta = merged_meta)
}

aggregated <- lapply(names(boundaries), function(boundary_name) aggregate_boundary(boundary_name, boundaries[[boundary_name]]))
names(aggregated) <- names(boundaries)

cell_count_audit <- do.call(rbind, lapply(aggregated, function(x) x$meta))
write_csv(cell_count_audit, "1_boundary_sample_cell_counts.csv")

for (boundary_name in names(aggregated)) {
  pb <- as.matrix(aggregated[[boundary_name]]$counts)
  write.csv(data.frame(Gene = rownames(pb), pb, check.names = FALSE), file.path(out_root, "tables", paste0("2_", boundary_name, "_full_pseudobulk_counts.csv")), row.names = FALSE)
}

comparisons <- list(
  all_older_vs_young = list(include = c("reference_young", "older_normal_BMI", "older_high_BMI"), reference = "reference_young", alternative = "older_all"),
  older_normal_vs_young = list(include = c("reference_young", "older_normal_BMI"), reference = "reference_young", alternative = "older_normal_BMI"),
  older_high_vs_young = list(include = c("reference_young", "older_high_BMI"), reference = "reference_young", alternative = "older_high_BMI"),
  older_high_vs_older_normal = list(include = c("older_normal_BMI", "older_high_BMI"), reference = "older_normal_BMI", alternative = "older_high_BMI")
)

make_comparison_group <- function(meta_df, comparison) {
  group <- meta_df$bmi_group_curated
  keep <- group %in% comparison$include
  group[group == "reference_young" & comparison$alternative == "older_all"] <- "reference_young"
  group[group %in% c("older_normal_BMI", "older_high_BMI") & comparison$alternative == "older_all"] <- "older_all"
  factor(group[keep], levels = c(comparison$reference, comparison$alternative))
}

pseudocounts <- c(0.1, 0.5, 1)
hmgcs2_rows <- list()
detection_rows <- list()
detection_test_rows <- list()
raw_target_rows <- list()
program_rows <- list()
edge_rows <- list()
voom_rows <- list()

for (boundary_name in names(aggregated)) {
  pb <- as.matrix(aggregated[[boundary_name]]$counts)
  sm <- aggregated[[boundary_name]]$meta
  lib_size <- colSums(pb)

  target_present <- intersect(target_genes, rownames(pb))
  raw_target <- data.frame(boundary = boundary_name, Gene = target_present, pb[target_present, , drop = FALSE], check.names = FALSE)
  raw_target_rows[[boundary_name]] <- raw_target

  hmgcs2_raw <- as.numeric(pb["HMGCS2", ])
  detection_rows[[boundary_name]] <- do.call(rbind, lapply(unique(sm$bmi_group_curated), function(stratum) {
    idx <- sm$bmi_group_curated == stratum
    data.frame(
      boundary = boundary_name,
      stratum = stratum,
      n_samples = sum(idx),
      mean_raw_count = mean(hmgcs2_raw[idx]),
      median_raw_count = median(hmgcs2_raw[idx]),
      n_detected = sum(hmgcs2_raw[idx] > 0),
      detection_rate = mean(hmgcs2_raw[idx] > 0),
      n_zero = sum(hmgcs2_raw[idx] == 0),
      zero_rate = mean(hmgcs2_raw[idx] == 0)
    )
  }))
  detection_test_rows[[boundary_name]] <- do.call(rbind, lapply(names(comparisons), function(comparison_name) {
    comparison <- comparisons[[comparison_name]]
    keep <- sm$bmi_group_curated %in% comparison$include
    group <- make_comparison_group(sm, comparison)
    detected <- hmgcs2_raw[keep] > 0
    reference_detected <- sum(detected[group == comparison$reference])
    reference_zero <- sum(!detected[group == comparison$reference])
    alternative_detected <- sum(detected[group == comparison$alternative])
    alternative_zero <- sum(!detected[group == comparison$alternative])
    fisher_p <- fisher.test(
      matrix(
        c(reference_detected, reference_zero, alternative_detected, alternative_zero),
        nrow = 2,
        byrow = TRUE
      )
    )$p.value
    data.frame(
      boundary = boundary_name,
      comparison = comparison_name,
      reference_detected = reference_detected,
      reference_zero = reference_zero,
      alternative_detected = alternative_detected,
      alternative_zero = alternative_zero,
      reference_detection_rate = mean(detected[group == comparison$reference]),
      alternative_detection_rate = mean(detected[group == comparison$alternative]),
      fisher_exact_p_two_sided = fisher_p
    )
  }))

  for (pc in pseudocounts) {
    log_cpm <- log2(sweep(pb + pc, 2, lib_size + 2 * pc, "/") * 1e6)
    hmgcs2 <- as.numeric(log_cpm["HMGCS2", ])
    for (comparison_name in names(comparisons)) {
      comparison <- comparisons[[comparison_name]]
      keep <- sm$bmi_group_curated %in% comparison$include
      group <- make_comparison_group(sm, comparison)
      x <- hmgcs2[keep]
      perm <- exact_permutation(x, group, comparison$reference)
      ci <- bootstrap_ci(x, group, comparison$reference)
      lm_result <- lm_safe(x, group, comparison$reference)
      hmgcs2_rows[[paste(boundary_name, pc, comparison_name, sep = "__")]] <- data.frame(
        boundary = boundary_name,
        pseudocount = pc,
        comparison = comparison_name,
        n_reference = sum(group == comparison$reference),
        n_alternative = sum(group == comparison$alternative),
        mean_reference = mean(x[group == comparison$reference]),
        mean_alternative = mean(x[group == comparison$alternative]),
        effect_alternative_minus_reference = unname(perm["effect"]),
        bootstrap_ci_low = ci[1],
        bootstrap_ci_high = ci[2],
        wilcox_p = wilcox_safe(x, group),
        lm_beta = unname(lm_result["beta"]),
        lm_p = unname(lm_result["p"]),
        exact_permutation_p_two_sided = unname(perm["p_two_sided"]),
        exact_permutation_p_lower_one_sided = unname(perm["p_lower_one_sided"]),
        n_exact_label_allocations = unname(perm["n_permutations"])
      )
    }
  }

  dge <- DGEList(counts = pb, samples = sm)
  dge <- calcNormFactors(dge)
  tmm_logcpm <- cpm(dge, log = TRUE, prior.count = 0.5)
  for (comparison_name in names(comparisons)) {
    comparison <- comparisons[[comparison_name]]
    keep_samples <- sm$bmi_group_curated %in% comparison$include
    sm_sub <- sm[keep_samples, , drop = FALSE]
    group <- make_comparison_group(sm, comparison)
    dge_sub <- DGEList(counts = pb[, keep_samples, drop = FALSE], group = group)
    keep_genes <- filterByExpr(dge_sub, group = group)
    dge_sub <- dge_sub[keep_genes, , keep.lib.sizes = FALSE]
    dge_sub <- calcNormFactors(dge_sub)
    design <- model.matrix(~ group)
    dge_sub <- estimateDisp(dge_sub, design)
    edge_fit <- glmQLFit(dge_sub, design, robust = TRUE)
    edge_test <- glmQLFTest(edge_fit, coef = 2)
    edge_tab <- topTags(edge_test, n = Inf, sort.by = "none")$table
    edge_tab$Gene <- rownames(edge_tab)
    edge_tab$boundary <- boundary_name
    edge_tab$comparison <- comparison_name
    edge_tab$passed_filterByExpr <- TRUE
    missing_targets <- setdiff(target_genes, edge_tab$Gene)
    if (length(missing_targets) > 0) {
      edge_tab <- rbind(
        edge_tab,
        data.frame(
          logFC = NA_real_, logCPM = NA_real_, F = NA_real_, PValue = NA_real_, FDR = NA_real_,
          Gene = missing_targets, boundary = boundary_name, comparison = comparison_name, passed_filterByExpr = FALSE
        )
      )
    }
    edge_rows[[paste(boundary_name, comparison_name, sep = "__")]] <- edge_tab

    v <- voom(dge_sub, design, plot = FALSE)
    voom_fit <- eBayes(lmFit(v, design), robust = TRUE)
    voom_tab <- topTable(voom_fit, coef = 2, number = Inf, sort.by = "none")
    voom_tab$Gene <- rownames(voom_tab)
    voom_tab$boundary <- boundary_name
    voom_tab$comparison <- comparison_name
    voom_tab$passed_filterByExpr <- TRUE
    missing_targets <- setdiff(target_genes, voom_tab$Gene)
    if (length(missing_targets) > 0) {
      voom_tab <- rbind(
        voom_tab,
        data.frame(
          logFC = NA_real_, AveExpr = NA_real_, t = NA_real_, P.Value = NA_real_, adj.P.Val = NA_real_, B = NA_real_,
          Gene = missing_targets, boundary = boundary_name, comparison = comparison_name, passed_filterByExpr = FALSE
        )
      )
    }
    voom_rows[[paste(boundary_name, comparison_name, sep = "__")]] <- voom_tab

    for (program in names(program_sets)) {
      genes_use <- intersect(program_sets[[program]], rownames(tmm_logcpm))
      score <- colMeans(tmm_logcpm[genes_use, , drop = FALSE])
      x <- score[keep_samples]
      perm <- exact_permutation(x, group, comparison$reference)
      ci <- bootstrap_ci(x, group, comparison$reference)
      lm_result <- lm_safe(x, group, comparison$reference)
      program_rows[[paste(boundary_name, comparison_name, program, sep = "__")]] <- data.frame(
        boundary = boundary_name,
        comparison = comparison_name,
        program = program,
        n_genes_present = length(genes_use),
        genes_present = paste(genes_use, collapse = ";"),
        effect_alternative_minus_reference = unname(perm["effect"]),
        bootstrap_ci_low = ci[1],
        bootstrap_ci_high = ci[2],
        wilcox_p = wilcox_safe(x, group),
        lm_beta = unname(lm_result["beta"]),
        lm_p = unname(lm_result["p"]),
        exact_permutation_p_two_sided = unname(perm["p_two_sided"]),
        exact_permutation_p_lower_one_sided = unname(perm["p_lower_one_sided"]),
        n_exact_label_allocations = unname(perm["n_permutations"])
      )
    }
  }
}

hmgcs2_table <- do.call(rbind, hmgcs2_rows)
detection_table <- do.call(rbind, detection_rows)
detection_test_table <- do.call(rbind, detection_test_rows)
raw_target_table <- do.call(rbind, raw_target_rows)
program_table <- do.call(rbind, program_rows)
edge_table <- do.call(rbind, edge_rows)
voom_table <- do.call(rbind, voom_rows)

write_csv(hmgcs2_table, "3_HMGCS2_pseudocount_BMI_exact_permutation_summary.csv")
write_csv(detection_table, "4_HMGCS2_detection_zero_rate_by_BMI_stratum.csv")
write_csv(detection_test_table, "04b_HMGCS2_detection_Fisher_tests.csv")
write_csv(raw_target_table, "5_target_gene_raw_counts_by_boundary.csv")
write_csv(program_table, "6_curated_program_TMM_logCPM_BMI_sensitivity.csv")
write_csv(edge_table, "7_edgeR_QLF_all_genes.csv")
write_csv(edge_table[edge_table$Gene %in% target_genes, ], "8_edgeR_QLF_target_genes.csv")
write_csv(voom_table, "9_limma_voom_all_genes.csv")
write_csv(voom_table[voom_table$Gene %in% target_genes, ], "10_limma_voom_target_genes.csv")

plot_values <- list()
for (boundary_name in names(aggregated)) {
  pb <- as.matrix(aggregated[[boundary_name]]$counts)
  sm <- aggregated[[boundary_name]]$meta
  lib_size <- colSums(pb)
  hmgcs2_logcpm <- as.numeric(log2(((pb["HMGCS2", ] + 0.5) / (lib_size + 1)) * 1e6))
  plot_values[[boundary_name]] <- data.frame(
    boundary = boundary_name,
    sample_id = sm$sample_id,
    stratum = sm$bmi_group_curated,
    raw_count = as.numeric(pb["HMGCS2", ]),
    logCPM_pc0.5 = hmgcs2_logcpm
  )
}
plot_values <- do.call(rbind, plot_values)
write_csv(plot_values, "11_HMGCS2_sample_level_plot_values.csv")

stratum_labels <- c(
  reference_young = "Young reference",
  older_normal_BMI = "Older normal BMI",
  older_high_BMI = "Older high BMI"
)
stratum_colours <- c(
  reference_young = "#3B6E8F",
  older_normal_BMI = "#C86B4A",
  older_high_BMI = "#8B4C55"
)
boundary_labels <- c(
  CMB_0_17 = "CMB (0+17)",
  ESB_0_14_17 = "ESB (0+14+17)"
)
plot_values$stratum <- factor(plot_values$stratum, levels = names(stratum_labels), labels = unname(stratum_labels))

p_raw <- ggplot(plot_values, aes(x = stratum, y = raw_count, colour = stratum)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.42) +
  geom_jitter(width = 0.11, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~ boundary, scales = "free_y", labeller = as_labeller(boundary_labels)) +
  scale_colour_manual(values = unname(stratum_colours), guide = "none") +
  labs(x = NULL, y = "HMGCS2 raw pseudobulk counts") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_raw, file.path(out_root, "figures", "Fig_S1_HMGCS2_raw_counts_by_BMI_stratum"), 170, 70)

p_log <- ggplot(plot_values, aes(x = stratum, y = logCPM_pc0.5, colour = stratum)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.42) +
  geom_jitter(width = 0.11, height = 0, size = 1.8, alpha = 0.9) +
  facet_wrap(~ boundary, scales = "free_y", labeller = as_labeller(boundary_labels)) +
  scale_colour_manual(values = unname(stratum_colours), guide = "none") +
  labs(x = NULL, y = "HMGCS2 targeted logCPM") +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_log, file.path(out_root, "figures", "Fig_S2_HMGCS2_logCPM_by_BMI_stratum"), 170, 70)

forest <- hmgcs2_table[hmgcs2_table$pseudocount == 0.5 & hmgcs2_table$comparison != "older_high_vs_older_normal", ]
forest$comparison <- factor(
  forest$comparison,
  levels = c("all_older_vs_young", "older_normal_vs_young", "older_high_vs_young"),
  labels = c("All older vs young", "Older normal BMI vs young", "Older high BMI vs young")
)
p_forest <- ggplot(forest, aes(x = effect_alternative_minus_reference, y = comparison, colour = boundary)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55", linetype = 2) +
  geom_errorbarh(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), height = 0.12, linewidth = 0.5, position = position_dodge(width = 0.35)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c(CMB_0_17 = "#1F5A75", ESB_0_14_17 = "#B85B3F"), labels = boundary_labels) +
  labs(x = "Alternative-minus-reference HMGCS2 logCPM effect", y = NULL, colour = "Boundary") +
  theme_pub(8)
save_plot(p_forest, file.path(out_root, "figures", "Fig_S3_HMGCS2_BMI_sensitivity_forest"), 170, 72)

edge_target <- edge_table[edge_table$Gene %in% target_genes, ]
edge_target$method <- "edgeR QLF"
edge_target$p_value <- edge_target$PValue
edge_target$fdr <- edge_target$FDR
edge_target$effect <- edge_target$logFC
voom_target <- voom_table[voom_table$Gene %in% target_genes, ]
voom_target$method <- "limma-voom"
voom_target$p_value <- voom_target$P.Value
voom_target$fdr <- voom_target$adj.P.Val
voom_target$effect <- voom_target$logFC
alternative_target_summary <- rbind(
  edge_target[, c("boundary", "comparison", "Gene", "method", "passed_filterByExpr", "effect", "p_value", "fdr")],
  voom_target[, c("boundary", "comparison", "Gene", "method", "passed_filterByExpr", "effect", "p_value", "fdr")]
)
write_csv(alternative_target_summary, "12_alternative_model_target_gene_summary.csv")


capture.output(sessionInfo(), file = file.path(out_root, "logs", "sessionInfo.txt"))
