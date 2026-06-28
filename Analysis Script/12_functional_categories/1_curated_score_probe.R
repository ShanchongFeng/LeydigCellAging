suppressPackageStartupMessages({
  library(ggplot2)
})

work_dir <- path.expand(Sys.getenv(
  "FUNCTIONAL_CATEGORY_ROOT",
  unset = file.path(Sys.getenv("HOME"), "FUNCTIONAL_CATEGORY_ROOT")
))
locked_root <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
ortholog_root <- path.expand(Sys.getenv(
  "ORTHOLOG_CONSERVATION_ROOT",
  unset = file.path(Sys.getenv("HOME"), "ORTHOLOG_CONSERVATION_ROOT")
))
pseudobulk_root <- path.expand(Sys.getenv(
  "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT",
  unset = file.path(Sys.getenv("HOME"), "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT")
))
tables_dir <- file.path(work_dir, "tables")
figures_dir <- file.path(work_dir, "figures")
logs_dir <- file.path(work_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260619)
run_leave_one_gene_out <- FALSE

message("Curated category direction probe started: ", Sys.time())

curated_manifest <- file.path(
  locked_root,
  "tables/primary_locked/Hmgcs2_curated_sets/leydig_scRNA_phase2a__phase2b__hmgcs2_curated_sets__csv__00_curated_gene_sets_used.csv"
)
curated_manifest <- path.expand(Sys.getenv("CURATED_GENE_SET_MANIFEST", unset = curated_manifest))
ortholog_file <- path.expand(Sys.getenv(
  "ORTHOLOG_ONE_TO_ONE_FILE",
  unset = file.path(ortholog_root, "homologene_mouse_human_orthologs_one_to_one.csv")
))
pseudobulk_dir <- path.expand(Sys.getenv(
  "PORTABLE_SIGNATURE_PSEUDOBULK_DIR",
  unset = file.path(pseudobulk_root, "pseudobulk")
))
stopifnot(file.exists(curated_manifest), file.exists(ortholog_file), dir.exists(pseudobulk_dir))

read_csv_base <- function(path, row_names = FALSE) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = if (row_names) 1 else NULL)
}

manifest <- read_csv_base(curated_manifest)
orth <- read_csv_base(ortholog_file)
mouse_to_human <- setNames(orth$human_symbol, orth$mouse_symbol)

map_to_human <- function(genes) {
  mapped <- mouse_to_human[genes]
  mapped[is.na(mapped) | mapped == ""] <- toupper(genes[is.na(mapped) | mapped == ""])
  unique(unname(mapped))
}

manifest_sets <- split(manifest$gene, manifest$set)

hmgcs2_axis_mouse <- c(
  "Hmgcs2", "Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl",
  "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a"
)

mouse_sets <- list(
  hmgcs2_axis = hmgcs2_axis_mouse,
  hmgcs2_axis_no_hmgcs2 = setdiff(hmgcs2_axis_mouse, "Hmgcs2"),
  curated_fatty_acid_oxidation_score = manifest_sets$curated_fatty_acid_oxidation_score,
  curated_cholesterol_biosynthesis_score = manifest_sets$curated_cholesterol_biosynthesis_score,
  curated_steroidogenesis_score = manifest_sets$curated_steroidogenesis_score
)

mouse_sets$metabolic_support_union <- unique(c(
  mouse_sets$hmgcs2_axis,
  mouse_sets$curated_fatty_acid_oxidation_score,
  mouse_sets$curated_cholesterol_biosynthesis_score
))
mouse_sets$metabolic_support_union_no_hmgcs2 <- setdiff(mouse_sets$metabolic_support_union, "Hmgcs2")

human_sets <- lapply(mouse_sets, map_to_human)

gene_sets <- list(mouse = mouse_sets, human = human_sets)

write.csv(
  do.call(rbind, lapply(names(mouse_sets), function(set_name) {
    data.frame(
      set = set_name,
      mouse_gene = mouse_sets[[set_name]],
      human_gene = map_to_human(mouse_sets[[set_name]]),
      source = ifelse(grepl("^curated_", set_name), "existing_curated_manifest", "fixed_from_existing_score_sensitivity_script"),
      stringsAsFactors = FALSE
    )
  })),
  file.path(tables_dir, "curated_probe_gene_sets_used_20260619.csv"),
  row.names = FALSE
)

comparisons <- data.frame(
  dataset = c("GSE182786", "GSE182786", "GSE254315", "GSE254315", "GSE303193", "GSE303193"),
  boundary = c("CORE_0_17", "EXT_0_14_17", "CORE_1_6_18", "EXT_1_6_18_7", "CMB_16", "ESB_16_17_19"),
  species = c("human", "human", "human", "human", "mouse", "mouse"),
  stringsAsFactors = FALSE
)

safe_wilcox <- function(x, group) {
  if (!all(c("Aged", "Young") %in% group)) return(NA_real_)
  if (sum(group == "Aged") < 2 || sum(group == "Young") < 2) return(NA_real_)
  tryCatch(
    wilcox.test(x[group == "Aged"], x[group == "Young"], exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

score_gene_z_mean <- function(expr, genes) {
  available <- intersect(genes, colnames(expr))
  if (length(available) == 0) return(rep(NA_real_, nrow(expr)))
  sub <- as.matrix(expr[, available, drop = FALSE])
  z <- as.matrix(scale(sub))
  if (is.null(colnames(z))) colnames(z) <- available
  bad <- apply(z, 2, function(v) all(is.na(v)))
  if (any(bad)) z[, bad] <- NA_real_
  rowMeans(z, na.rm = TRUE)
}

calc_effects_for_comparison <- function(dataset, boundary, species) {
  prefix <- paste0(dataset, "__", boundary)
  expr_file <- file.path(pseudobulk_dir, paste0(prefix, "__sample_logcpm.csv"))
  meta_file <- file.path(pseudobulk_dir, paste0(prefix, "__sample_metadata.csv"))
  expr <- read_csv_base(expr_file, row_names = TRUE)
  meta <- read_csv_base(meta_file)

  common_samples <- intersect(rownames(expr), meta$sample)
  expr <- expr[common_samples, , drop = FALSE]
  meta <- meta[match(common_samples, meta$sample), , drop = FALSE]

  sets <- gene_sets[[species]]
  do.call(rbind, lapply(names(sets), function(set_name) {
    genes <- unique(sets[[set_name]])
    available <- intersect(genes, colnames(expr))
    score <- score_gene_z_mean(expr, genes)
    aged <- score[meta$group == "Aged"]
    young <- score[meta$group == "Young"]
    data.frame(
      analysis_layer = "curated_probe_gene_z_mean",
      dataset = dataset,
      boundary = boundary,
      species = species,
      gene_set = set_name,
      n_genes_defined = length(genes),
      n_genes_available = length(available),
      coverage = length(available) / max(1, length(genes)),
      n_young = length(young),
      n_aged = length(aged),
      young_mean = mean(young, na.rm = TRUE),
      aged_mean = mean(aged, na.rm = TRUE),
      aged_minus_young = mean(aged, na.rm = TRUE) - mean(young, na.rm = TRUE),
      wilcox_p = safe_wilcox(score, meta$group),
      direction = ifelse(mean(aged, na.rm = TRUE) < mean(young, na.rm = TRUE), "aged_lower", "aged_higher"),
      genes_available = paste(available, collapse = ";"),
      genes_missing = paste(setdiff(genes, available), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}

effect_list <- lapply(seq_len(nrow(comparisons)), function(i) {
  calc_effects_for_comparison(comparisons$dataset[i], comparisons$boundary[i], comparisons$species[i])
})
effects <- do.call(rbind, effect_list)
effects$wilcox_fdr_within_probe <- p.adjust(effects$wilcox_p, method = "BH")

write.csv(
  effects,
  file.path(tables_dir, "curated_category_effect_summary_20260619.csv"),
  row.names = FALSE
)

dataset_boundary <- aggregate(
  aged_minus_young ~ dataset + boundary + gene_set,
  effects,
  median,
  na.rm = TRUE
)
dataset_boundary$direction <- ifelse(dataset_boundary$aged_minus_young < 0, "aged_lower", "aged_higher")

dataset_level <- aggregate(
  aged_minus_young ~ dataset + gene_set,
  effects,
  median,
  na.rm = TRUE
)
dataset_level$direction <- ifelse(dataset_level$aged_minus_young < 0, "aged_lower", "aged_higher")

dataset_blocked_summary <- do.call(rbind, lapply(split(dataset_level, dataset_level$gene_set), function(df) {
  data.frame(
    gene_set = df$gene_set[1],
    n_dataset_blocks = nrow(df),
    n_dataset_aged_lower = sum(df$direction == "aged_lower"),
    aged_lower_fraction_dataset_blocked = mean(df$direction == "aged_lower"),
    median_dataset_effect = median(df$aged_minus_young, na.rm = TRUE),
    min_dataset_effect = min(df$aged_minus_young, na.rm = TRUE),
    max_dataset_effect = max(df$aged_minus_young, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  dataset_level,
  file.path(tables_dir, "curated_category_dataset_level_effects_20260619.csv"),
  row.names = FALSE
)
write.csv(
  dataset_blocked_summary,
  file.path(tables_dir, "curated_category_dataset_blocked_summary_20260619.csv"),
  row.names = FALSE
)

calc_leave_one_gene <- function(dataset, boundary, species, set_name = "metabolic_support_union") {
  prefix <- paste0(dataset, "__", boundary)
  expr_file <- file.path(pseudobulk_dir, paste0(prefix, "__sample_logcpm.csv"))
  meta_file <- file.path(pseudobulk_dir, paste0(prefix, "__sample_metadata.csv"))
  expr <- read_csv_base(expr_file, row_names = TRUE)
  meta <- read_csv_base(meta_file)
  common_samples <- intersect(rownames(expr), meta$sample)
  expr <- expr[common_samples, , drop = FALSE]
  meta <- meta[match(common_samples, meta$sample), , drop = FALSE]

  genes <- unique(gene_sets[[species]][[set_name]])
  full_score <- score_gene_z_mean(expr, genes)
  full_effect <- mean(full_score[meta$group == "Aged"], na.rm = TRUE) -
    mean(full_score[meta$group == "Young"], na.rm = TRUE)

  do.call(rbind, lapply(genes, function(g) {
    reduced <- setdiff(genes, g)
    score <- score_gene_z_mean(expr, reduced)
    effect <- mean(score[meta$group == "Aged"], na.rm = TRUE) -
      mean(score[meta$group == "Young"], na.rm = TRUE)
    data.frame(
      dataset = dataset,
      boundary = boundary,
      species = species,
      gene_set = set_name,
      removed_gene = g,
      removed_gene_available = g %in% colnames(expr),
      full_effect = full_effect,
      leave_one_effect = effect,
      delta_from_full = effect - full_effect,
      direction_after_removal = ifelse(effect < 0, "aged_lower", "aged_higher"),
      stringsAsFactors = FALSE
    )
  }))
}

if (run_leave_one_gene_out) {
  loo <- do.call(rbind, lapply(seq_len(nrow(comparisons)), function(i) {
    calc_leave_one_gene(comparisons$dataset[i], comparisons$boundary[i], comparisons$species[i])
  }))

  write.csv(
    loo,
    file.path(tables_dir, "curated_metabolic_support_leave_one_gene_out_20260619.csv"),
    row.names = FALSE
  )
}

plot_df <- effects
plot_df$boundary_label <- paste(plot_df$dataset, plot_df$boundary, sep = "\n")
plot_df$gene_set_label <- factor(
  plot_df$gene_set,
  levels = rev(c(
    "hmgcs2_axis",
    "hmgcs2_axis_no_hmgcs2",
    "curated_fatty_acid_oxidation_score",
    "curated_cholesterol_biosynthesis_score",
    "metabolic_support_union",
    "metabolic_support_union_no_hmgcs2",
    "curated_steroidogenesis_score"
  )),
  labels = rev(c(
    "HMGCS2 axis",
    "HMGCS2 axis, no HMGCS2",
    "FAO",
    "Cholesterol biosynthesis",
    "Metabolic support union",
    "Metabolic support union, no HMGCS2",
    "Steroidogenesis"
  ))
)
plot_df$effect_capped <- pmax(pmin(plot_df$aged_minus_young, 1.5), -1.5)
plot_df$nominal_sig <- !is.na(plot_df$wilcox_p) & plot_df$wilcox_p < 0.05
plot_df$fdr_sig <- !is.na(plot_df$wilcox_fdr_within_probe) & plot_df$wilcox_fdr_within_probe < 0.10
plot_df$alpha_value <- ifelse(plot_df$coverage >= 0.75, 1, 0.35)

p <- ggplot(plot_df, aes(x = boundary_label, y = gene_set_label, fill = effect_capped, alpha = alpha_value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_point(
    data = subset(plot_df, nominal_sig),
    aes(x = boundary_label, y = gene_set_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.0,
    stroke = 0.35,
    fill = "#2E2A27",
    color = "#2E2A27"
  ) +
  geom_point(
    data = subset(plot_df, fdr_sig),
    aes(x = boundary_label, y = gene_set_label),
    inherit.aes = FALSE,
    shape = 8,
    size = 2.0,
    stroke = 0.45,
    color = "#2E2A27"
  ) +
  scale_fill_gradient2(
    low = "#77a1ec",
    mid = "#f7e3db",
    high = "#f48f89",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    name = "Aged - young\nmean gene z"
  ) +
  scale_alpha_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "#2E2A27"),
    axis.text.y = element_text(color = "#2E2A27"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(
  file.path(figures_dir, "curated_category_direction_heatmap_20260619.pdf"),
  p,
  width = 183,
  height = 120,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  file.path(figures_dir, "curated_category_direction_heatmap_20260619.png"),
  p,
  width = 183,
  height = 120,
  units = "mm",
  dpi = 450
)

sink(file.path(logs_dir, "curated_probe_session_info_20260619.txt"))
cat("Curated category direction probe completed: ", as.character(Sys.time()), "\n\n", sep = "")
print(sessionInfo())
sink()

message("Curated category direction probe completed: ", Sys.time())
