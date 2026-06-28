suppressPackageStartupMessages({
  library(ggplot2)
})

work_dir <- path.expand(Sys.getenv(
  "FUNCTIONAL_CATEGORY_ROOT",
  unset = file.path(Sys.getenv("HOME"), "FUNCTIONAL_CATEGORY_ROOT")
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

message("Reactome external category gene-effect summary started: ", Sys.time())

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

effect_file <- path.expand(Sys.getenv(
  "ORTHOLOG_GENE_EFFECTS_FILE",
  unset = file.path(ortholog_root, "ortholog_gene_effects.csv")
))
gene_set_file <- file.path(tables_dir, "reactome_external_gene_sets_20260619.csv")
stopifnot(file.exists(effect_file), file.exists(gene_set_file))
effects <- read_csv(effect_file)
gene_sets <- read_csv(gene_set_file)
gene_sets$entrez_id <- as.character(gene_sets$entrez_id)
effects$mouse_entrez <- as.character(effects$mouse_entrez)
effects$human_entrez <- as.character(effects$human_entrez)

safe_signed_rank_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  if (all(x == 0)) return(NA_real_)
  tryCatch(wilcox.test(x, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
}

blocks <- unique(effects[, c("dataset", "boundary", "species", "aggregation")])
categories <- unique(gene_sets$category)

rows <- list()
k <- 1L
for (i in seq_len(nrow(blocks))) {
  block <- blocks[i, ]
  block_effects <- effects[
    effects$dataset == block$dataset &
      effects$boundary == block$boundary &
      effects$species == block$species &
      effects$aggregation == block$aggregation,
  ]
  entrez_col <- ifelse(block$species == "human", "human_entrez", "mouse_entrez")
  for (category in categories) {
    entrez_ids <- unique(gene_sets$entrez_id[gene_sets$species == block$species & gene_sets$category == category])
    sub <- block_effects[block_effects[[entrez_col]] %in% entrez_ids, ]
    vals <- sub$aged_minus_young
    rows[[k]] <- data.frame(
      analysis_layer = "reactome_external_gene_effect_table",
      dataset = block$dataset,
      boundary = block$boundary,
      species = block$species,
      aggregation = block$aggregation,
      category = category,
      n_genes_defined = length(entrez_ids),
      n_genes_available = length(unique(sub$assay_symbol)),
      coverage = length(unique(sub[[entrez_col]])) / max(1, length(entrez_ids)),
      median_gene_effect = median(vals, na.rm = TRUE),
      mean_gene_effect = mean(vals, na.rm = TRUE),
      fraction_genes_aged_lower = mean(vals < 0, na.rm = TRUE),
      n_genes_aged_lower = sum(vals < 0, na.rm = TRUE),
      n_nominal_gene_p_lt_0_05 = sum(sub$p_value < 0.05, na.rm = TRUE),
      signed_rank_p_against_zero = safe_signed_rank_p(vals),
      direction_by_median = ifelse(median(vals, na.rm = TRUE) < 0, "aged_lower", "aged_higher"),
      genes_available = paste(sort(unique(sub$assay_symbol)), collapse = ";"),
      missing_entrez = paste(sort(setdiff(entrez_ids, unique(sub[[entrez_col]]))), collapse = ";"),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

summary_df <- do.call(rbind, rows)
summary_df$signed_rank_fdr_within_reactome_probe <- p.adjust(summary_df$signed_rank_p_against_zero, method = "BH")
write.csv(summary_df, file.path(tables_dir, "reactome_external_category_gene_effect_summary_20260619.csv"), row.names = FALSE)

sample_summary <- summary_df[summary_df$aggregation == "sample_mean", ]
dataset_level <- aggregate(
  median_gene_effect ~ dataset + category,
  sample_summary,
  median,
  na.rm = TRUE
)
dataset_level$direction <- ifelse(dataset_level$median_gene_effect < 0, "aged_lower", "aged_higher")

dataset_blocked <- do.call(rbind, lapply(split(dataset_level, dataset_level$category), function(df) {
  data.frame(
    category = df$category[1],
    n_dataset_blocks = nrow(df),
    n_dataset_aged_lower = sum(df$direction == "aged_lower"),
    aged_lower_fraction_dataset_blocked = mean(df$direction == "aged_lower"),
    median_dataset_effect = median(df$median_gene_effect, na.rm = TRUE),
    min_dataset_effect = min(df$median_gene_effect, na.rm = TRUE),
    max_dataset_effect = max(df$median_gene_effect, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(dataset_level, file.path(tables_dir, "reactome_external_category_dataset_level_20260619.csv"), row.names = FALSE)
write.csv(dataset_blocked, file.path(tables_dir, "reactome_external_category_dataset_blocked_20260619.csv"), row.names = FALSE)

plot_order <- c(
  "reactome_fatty_acid_beta_oxidation",
  "reactome_ketone_body_metabolism",
  "reactome_mitochondrial_energy",
  "reactome_cholesterol_handling",
  "reactome_oxidative_stress_redox",
  "reactome_metabolic_support_union_prespecified",
  "reactome_steroidogenic_execution",
  "reactome_steroid_metabolism_broad"
)
plot_labels <- c(
  "FA beta-oxidation",
  "Ketone-body metabolism",
  "Mitochondrial energy",
  "Cholesterol handling",
  "Oxidative stress/redox",
  "Metabolic-support union",
  "Steroidogenic execution",
  "Steroid metabolism, broad"
)
names(plot_labels) <- plot_order

plot_df <- summary_df[summary_df$category %in% plot_order, ]
plot_df$block_label <- paste(plot_df$dataset, plot_df$boundary, sep = "\n")
plot_df$category_label <- factor(plot_labels[plot_df$category], levels = rev(plot_labels))
plot_df$effect_capped <- pmax(pmin(plot_df$median_gene_effect, 0.75), -0.75)
plot_df$alpha_value <- ifelse(plot_df$coverage >= 0.25 & plot_df$n_genes_available >= 3, 1, 0.35)
plot_df$is_sample_mean <- plot_df$aggregation == "sample_mean"
plot_df$nominal_category_p <- !is.na(plot_df$signed_rank_p_against_zero) & plot_df$signed_rank_p_against_zero < 0.05

p <- ggplot(plot_df, aes(x = block_label, y = category_label, fill = effect_capped, alpha = alpha_value)) +
  geom_tile(color = "white", size = 0.35) +
  geom_point(
    data = subset(plot_df, nominal_category_p),
    aes(x = block_label, y = category_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.7,
    stroke = 0.35,
    fill = "#2E2A27",
    color = "#2E2A27"
  ) +
  geom_point(
    data = subset(plot_df, !is_sample_mean),
    aes(x = block_label, y = category_label),
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
    limits = c(-0.75, 0.75),
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

ggsave(file.path(figures_dir, "reactome_external_category_gene_effect_heatmap_20260619.pdf"), p, width = 183, height = 120, units = "mm")
ggsave(file.path(figures_dir, "reactome_external_category_gene_effect_heatmap_20260619.png"), p, width = 183, height = 120, units = "mm", dpi = 450)

sink(file.path(logs_dir, "reactome_external_gene_effect_summary_session_info_20260619.txt"))
cat("Reactome external category gene-effect summary completed: ", as.character(Sys.time()), "\n\n", sep = "")
print(sessionInfo())
sink()

message("Reactome external category gene-effect summary completed: ", Sys.time())
