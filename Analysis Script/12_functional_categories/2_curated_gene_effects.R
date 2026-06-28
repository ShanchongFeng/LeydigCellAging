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
tables_dir <- file.path(work_dir, "tables")
figures_dir <- file.path(work_dir, "figures")
logs_dir <- file.path(work_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

message("Curated category gene-effect summary started: ", Sys.time())

effect_file <- path.expand(Sys.getenv(
  "ORTHOLOG_GENE_EFFECTS_FILE",
  unset = file.path(ortholog_root, "ortholog_gene_effects.csv")
))
curated_manifest <- file.path(
  locked_root,
  "tables/primary_locked/Hmgcs2_curated_sets/leydig_scRNA_phase2a__phase2b__hmgcs2_curated_sets__csv__00_curated_gene_sets_used.csv"
)
curated_manifest <- path.expand(Sys.getenv("CURATED_GENE_SET_MANIFEST", unset = curated_manifest))
ortholog_file <- path.expand(Sys.getenv(
  "ORTHOLOG_ONE_TO_ONE_FILE",
  unset = file.path(ortholog_root, "homologene_mouse_human_orthologs_one_to_one.csv")
))
stopifnot(file.exists(effect_file), file.exists(curated_manifest), file.exists(ortholog_file))

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

effects <- read_csv(effect_file)
manifest <- read_csv(curated_manifest)
orth <- read_csv(ortholog_file)
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
  hmgcs2_single_gene = "Hmgcs2",
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

gene_set_manifest <- do.call(rbind, lapply(names(mouse_sets), function(set_name) {
  n <- max(length(mouse_sets[[set_name]]), length(map_to_human(mouse_sets[[set_name]])))
  data.frame(
    set = set_name,
    mouse_gene = c(mouse_sets[[set_name]], rep(NA_character_, n - length(mouse_sets[[set_name]]))),
    human_gene = c(map_to_human(mouse_sets[[set_name]]), rep(NA_character_, n - length(map_to_human(mouse_sets[[set_name]])))),
    stringsAsFactors = FALSE
  )
}))
write.csv(gene_set_manifest, file.path(tables_dir, "curated_gene_effect_sets_used_20260619.csv"), row.names = FALSE)

safe_signed_rank_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (all(x == 0)) return(NA_real_)
  tryCatch(wilcox.test(x, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
}

summaries <- list()
k <- 1L
blocks <- unique(effects[, c("dataset", "boundary", "species", "aggregation")])

for (i in seq_len(nrow(blocks))) {
  block <- blocks[i, ]
  block_effects <- effects[
    effects$dataset == block$dataset &
      effects$boundary == block$boundary &
      effects$species == block$species &
      effects$aggregation == block$aggregation,
  ]
  sets <- gene_sets[[block$species]]
  for (set_name in names(sets)) {
    genes <- unique(sets[[set_name]])
    sub <- block_effects[block_effects$assay_symbol %in% genes, ]
    vals <- sub$aged_minus_young
    summaries[[k]] <- data.frame(
      analysis_layer = "curated_probe_gene_effect_table",
      dataset = block$dataset,
      boundary = block$boundary,
      species = block$species,
      aggregation = block$aggregation,
      gene_set = set_name,
      n_genes_defined = length(genes),
      n_genes_available = length(unique(sub$assay_symbol)),
      coverage = length(unique(sub$assay_symbol)) / max(1, length(genes)),
      median_gene_effect = median(vals, na.rm = TRUE),
      mean_gene_effect = mean(vals, na.rm = TRUE),
      fraction_genes_aged_lower = mean(vals < 0, na.rm = TRUE),
      n_genes_aged_lower = sum(vals < 0, na.rm = TRUE),
      n_nominal_p_lt_0_05 = sum(sub$p_value < 0.05, na.rm = TRUE),
      signed_rank_p_against_zero = safe_signed_rank_p(vals),
      direction_by_median = ifelse(median(vals, na.rm = TRUE) < 0, "aged_lower", "aged_higher"),
      genes_available = paste(sort(unique(sub$assay_symbol)), collapse = ";"),
      genes_missing = paste(sort(setdiff(genes, unique(sub$assay_symbol))), collapse = ";"),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

summary_df <- do.call(rbind, summaries)
summary_df$signed_rank_fdr_within_probe <- p.adjust(summary_df$signed_rank_p_against_zero, method = "BH")
write.csv(summary_df, file.path(tables_dir, "curated_category_gene_effect_summary_20260619.csv"), row.names = FALSE)

dataset_level <- aggregate(
  median_gene_effect ~ dataset + gene_set,
  summary_df[summary_df$aggregation == "sample_mean", ],
  median,
  na.rm = TRUE
)
dataset_level$direction <- ifelse(dataset_level$median_gene_effect < 0, "aged_lower", "aged_higher")

dataset_blocked <- do.call(rbind, lapply(split(dataset_level, dataset_level$gene_set), function(df) {
  data.frame(
    gene_set = df$gene_set[1],
    n_dataset_blocks = nrow(df),
    n_dataset_aged_lower = sum(df$direction == "aged_lower"),
    aged_lower_fraction_dataset_blocked = mean(df$direction == "aged_lower"),
    median_dataset_effect = median(df$median_gene_effect, na.rm = TRUE),
    min_dataset_effect = min(df$median_gene_effect, na.rm = TRUE),
    max_dataset_effect = max(df$median_gene_effect, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(dataset_level, file.path(tables_dir, "curated_category_gene_effect_dataset_level_20260619.csv"), row.names = FALSE)
write.csv(dataset_blocked, file.path(tables_dir, "curated_category_gene_effect_dataset_blocked_20260619.csv"), row.names = FALSE)

plot_order <- c(
  "hmgcs2_single_gene",
  "hmgcs2_axis",
  "hmgcs2_axis_no_hmgcs2",
  "curated_fatty_acid_oxidation_score",
  "curated_cholesterol_biosynthesis_score",
  "metabolic_support_union",
  "metabolic_support_union_no_hmgcs2",
  "curated_steroidogenesis_score"
)
plot_labels <- c(
  "HMGCS2",
  "HMGCS2 axis",
  "HMGCS2 axis, no HMGCS2",
  "FAO",
  "Cholesterol biosynthesis",
  "Metabolic support union",
  "Metabolic support union, no HMGCS2",
  "Steroidogenesis"
)
names(plot_labels) <- plot_order

plot_df <- summary_df[summary_df$gene_set %in% plot_order, ]
plot_df$block_label <- paste(plot_df$dataset, plot_df$boundary, sep = "\n")
plot_df$gene_set_label <- factor(plot_labels[plot_df$gene_set], levels = rev(plot_labels))
plot_df$effect_capped <- pmax(pmin(plot_df$median_gene_effect, 1.5), -1.5)
plot_df$alpha_value <- ifelse(plot_df$coverage >= 0.75, 1, 0.35)
plot_df$is_sample_mean <- plot_df$aggregation == "sample_mean"
plot_df$nominal_category_p <- !is.na(plot_df$signed_rank_p_against_zero) & plot_df$signed_rank_p_against_zero < 0.05

p <- ggplot(plot_df, aes(x = block_label, y = gene_set_label, fill = effect_capped, alpha = alpha_value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_point(
    data = subset(plot_df, nominal_category_p),
    aes(x = block_label, y = gene_set_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.7,
    stroke = 0.35,
    fill = "#2E2A27",
    color = "#2E2A27"
  ) +
  geom_point(
    data = subset(plot_df, !is_sample_mean),
    aes(x = block_label, y = gene_set_label),
    inherit.aes = FALSE,
    shape = 4,
    size = 2.0,
    stroke = 0.4,
    color = "#2E2A27"
  ) +
  scale_fill_gradient2(
    low = "#77a1ec",
    mid = "#f7e3db",
    high = "#f48f89",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    name = "Median gene effect\naged - young"
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

ggsave(file.path(figures_dir, "curated_category_gene_effect_heatmap_20260619.pdf"), p, width = 183, height = 120, units = "mm")
ggsave(file.path(figures_dir, "curated_category_gene_effect_heatmap_20260619.png"), p, width = 183, height = 120, units = "mm", dpi = 450)

sink(file.path(logs_dir, "curated_gene_effect_summary_session_info_20260619.txt"))
cat("Curated category gene-effect summary completed: ", as.character(Sys.time()), "\n\n", sep = "")
print(sessionInfo())
sink()

message("Curated category gene-effect summary completed: ", Sys.time())
