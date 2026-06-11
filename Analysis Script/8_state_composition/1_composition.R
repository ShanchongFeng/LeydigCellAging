#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260602)

LOCKED_ROOT <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
OUT <- path.expand(Sys.getenv(
  "STATE_STRUCTURE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "STATE_STRUCTURE_ROOT")
))
TABLE_DIR <- file.path(OUT, "tables", "1_state_composition")
FIG_DIR <- file.path(OUT, "figures", "exploratory")
DOC_DIR <- file.path(OUT, "docs")
LOG_DIR <- file.path(OUT, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

n_boot <- 5000
n_perm_random <- 10000

configs <- list(
  list(
    dataset = "GSE270931",
    species = "mouse",
    sample_aware = FALSE,
    sample_col = "orig.ident",
    group_col = "group",
    rds = file.path(LOCKED_ROOT, "rds/primary_locked/GSE270931_mouse_scRNA/leydig_scRNA_phase2a__results__qc__04_filtered_clustered_seurat.rds"),
    boundaries = list(CMB_0_4 = c("0", "4"), ESB_0_4_14 = c("0", "4", "14")),
    state_labels = c("0" = "canonical_like", "4" = "aged_remodeled_like", "14" = "sensitivity_only")
  ),
  list(
    dataset = "GSE303193",
    species = "mouse",
    sample_aware = TRUE,
    sample_col = "sample_id",
    group_col = "group",
    rds = file.path(LOCKED_ROOT, "rds/primary_locked/GSE303193/GSE303193_phase1__results__03_filtered_clustered_seurat.rds"),
    boundaries = list(CMB_16 = c("16"), ESB_16_17_19 = c("16", "17", "19")),
    state_labels = c("16" = "canonical_like_core", "17" = "sensitivity_neighbor", "19" = "sensitivity_only")
  ),
  list(
    dataset = "GSE254315",
    species = "human",
    sample_aware = TRUE,
    sample_col = "sample_id",
    group_col = "group",
    rds = file.path(LOCKED_ROOT, "rds/primary_locked/GSE254315/GSE254315_check__results__03_filtered_clustered_human_seurat.rds"),
    boundaries = list(CMB_1_6_18 = c("1", "6", "18"), ESB_1_6_18_7 = c("1", "6", "18", "7")),
    state_labels = c("1" = "core_1", "6" = "core_6", "18" = "core_18", "7" = "sensitivity_only")
  ),
  list(
    dataset = "GSE182786",
    species = "human",
    sample_aware = TRUE,
    sample_col = "sample_id",
    group_col = "group",
    rds = file.path(LOCKED_ROOT, "rds/primary_locked/GSE182786_check/GSE182786_check__results_01_basic_clustered.rds"),
    boundaries = list(CMB_0_17 = c("0", "17"), ESB_0_14_17 = c("0", "14", "17")),
    state_labels = c("0" = "core_0", "17" = "core_17", "14" = "sensitivity_only")
  )
)

safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)

sample_bmi_stratum <- function(dataset, sample_unit, group) {
  if (dataset != "GSE182786") return(rep("not_applicable", length(sample_unit)))
  out <- rep("Young", length(sample_unit))
  older_num <- suppressWarnings(as.integer(sub("^Older", "", sample_unit)))
  out[group == "Aged" & older_num %in% 1:5] <- "Older_normal_BMI"
  out[group == "Aged" & older_num %in% 6:8] <- "Older_high_BMI"
  out[group == "Aged" & is.na(older_num)] <- "Older_BMI_unknown"
  out
}

all_label_effects <- function(values, groups) {
  keep <- is.finite(values) & groups %in% c("Young", "Aged")
  values <- values[keep]
  groups <- groups[keep]
  n_y <- sum(groups == "Young")
  n_a <- sum(groups == "Aged")
  if (n_y < 2 || n_a < 2) return(NULL)
  idx <- seq_along(values)
  observed <- mean(values[groups == "Aged"]) - mean(values[groups == "Young"])
  n_possible <- choose(length(idx), n_y)
  if (n_possible <= 5000) {
    combos <- combn(idx, n_y, simplify = FALSE)
    effects <- vapply(combos, function(young_idx) {
      aged_idx <- setdiff(idx, young_idx)
      mean(values[aged_idx]) - mean(values[young_idx])
    }, numeric(1))
    permutation_type <- "exact_enumeration"
  } else {
    effects <- replicate(n_perm_random, {
      young_idx <- sample(idx, n_y, replace = FALSE)
      aged_idx <- setdiff(idx, young_idx)
      mean(values[aged_idx]) - mean(values[young_idx])
    })
    permutation_type <- "random_permutation"
  }
  p <- mean(abs(effects) >= abs(observed) - 1e-12)
  list(
    observed = observed,
    permutation_p = p,
    n_permutations = length(effects),
    n_possible_labelings = n_possible,
    permutation_type = permutation_type
  )
}

bootstrap_effect <- function(values, groups, n_boot = 5000) {
  keep <- is.finite(values) & groups %in% c("Young", "Aged")
  values <- values[keep]
  groups <- groups[keep]
  young <- values[groups == "Young"]
  aged <- values[groups == "Aged"]
  if (length(young) < 2 || length(aged) < 2) return(c(low = NA_real_, high = NA_real_))
  effects <- replicate(n_boot, mean(sample(aged, replace = TRUE)) - mean(sample(young, replace = TRUE)))
  unname(quantile(effects, c(0.025, 0.975), na.rm = TRUE))
}

summarize_comparison <- function(df, dataset, boundary, denominator, comparison, group_override = NULL) {
  work <- df
  groups <- if (is.null(group_override)) work$group else group_override
  states <- unique(work$cluster)
  rows <- list()

  for (state in states) {
    sub <- work[work$cluster == state, , drop = FALSE]
    sub_groups <- groups[work$cluster == state]
    vals <- sub[[denominator]]
    perm <- all_label_effects(vals, sub_groups)
    if (is.null(perm)) next
    ci <- bootstrap_effect(vals, sub_groups, n_boot)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = dataset,
      boundary = boundary,
      denominator = denominator,
      comparison = comparison,
      cluster = state,
      state_label = sub$state_label[[1]],
      n_young = sum(sub_groups == "Young"),
      n_aged = sum(sub_groups == "Aged"),
      mean_young = safe_mean(vals[sub_groups == "Young"]),
      mean_aged = safe_mean(vals[sub_groups == "Aged"]),
      aged_minus_young = perm$observed,
      ci95_low = ci[[1]],
      ci95_high = ci[[2]],
      permutation_p = perm$permutation_p,
      permutation_type = perm$permutation_type,
      n_permutations = perm$n_permutations,
      n_possible_labelings = perm$n_possible_labelings,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

composition_rows <- list()
effect_rows <- list()
eligibility_rows <- list()

for (cfg in configs) {
  obj <- readRDS(cfg$rds)
  md <- obj@meta.data
  md$cluster <- as.character(md$seurat_clusters)
  md$sample_unit <- as.character(md[[cfg$sample_col]])
  md$group <- as.character(md[[cfg$group_col]])
  md$bmi_stratum <- sample_bmi_stratum(cfg$dataset, md$sample_unit, md$group)

  all_samples <- unique(md[, c("sample_unit", "group", "bmi_stratum"), drop = FALSE])
  all_samples$n_cells_whole_testis <- as.integer(table(factor(md$sample_unit, levels = all_samples$sample_unit)))

  for (boundary_name in names(cfg$boundaries)) {
    boundary_clusters <- cfg$boundaries[[boundary_name]]
    boundary_clusters <- boundary_clusters[boundary_clusters %in% unique(md$cluster)]
    if (length(boundary_clusters) == 0) next

    boundary_md <- md[md$cluster %in% boundary_clusters, , drop = FALSE]
    boundary_counts <- as.data.frame(table(
      sample_unit = factor(boundary_md$sample_unit, levels = all_samples$sample_unit),
      cluster = factor(boundary_md$cluster, levels = boundary_clusters)
    ), stringsAsFactors = FALSE)
    colnames(boundary_counts)[3] <- "n_cells_state"

    totals <- aggregate(n_cells_state ~ sample_unit, boundary_counts, sum)
    colnames(totals)[2] <- "n_cells_boundary"
    comp <- merge(boundary_counts, totals, by = "sample_unit", all.x = TRUE)
    comp <- merge(comp, all_samples, by = "sample_unit", all.x = TRUE)
    comp$dataset <- cfg$dataset
    comp$species <- cfg$species
    comp$sample_aware <- cfg$sample_aware
    comp$boundary <- boundary_name
    comp$state_label <- unname(cfg$state_labels[as.character(comp$cluster)])
    comp$state_label[is.na(comp$state_label)] <- paste0("cluster_", comp$cluster[is.na(comp$state_label)])
    comp$proportion_within_leydig <- ifelse(comp$n_cells_boundary > 0, comp$n_cells_state / comp$n_cells_boundary, NA_real_)
    comp$proportion_whole_testis <- ifelse(comp$n_cells_whole_testis > 0, comp$n_cells_state / comp$n_cells_whole_testis, NA_real_)
    comp$eligible_ge50 <- comp$n_cells_boundary >= 50
    comp$eligible_ge100 <- comp$n_cells_boundary >= 100
    comp <- comp[, c(
      "dataset", "species", "sample_aware", "boundary", "sample_unit", "group", "bmi_stratum",
      "cluster", "state_label", "n_cells_state", "n_cells_boundary", "n_cells_whole_testis",
      "proportion_within_leydig", "proportion_whole_testis", "eligible_ge50", "eligible_ge100"
    )]
    composition_rows[[length(composition_rows) + 1]] <- comp

    elig <- unique(comp[, c("dataset", "species", "sample_aware", "boundary", "sample_unit", "group", "bmi_stratum", "n_cells_boundary", "n_cells_whole_testis", "eligible_ge50", "eligible_ge100")])
    eligibility_rows[[length(eligibility_rows) + 1]] <- elig

    if (cfg$sample_aware) {
      primary <- comp[comp$eligible_ge50, , drop = FALSE]
      if (nrow(primary) > 0) {
        for (denom in c("proportion_within_leydig", "proportion_whole_testis")) {
          res <- summarize_comparison(primary, cfg$dataset, boundary_name, denom, "all_aged_vs_young")
          if (!is.null(res)) effect_rows[[length(effect_rows) + 1]] <- res
        }
      }

      strict <- comp[comp$eligible_ge100, , drop = FALSE]
      if (nrow(strict) > 0) {
        for (denom in c("proportion_within_leydig", "proportion_whole_testis")) {
          res <- summarize_comparison(strict, cfg$dataset, boundary_name, denom, "all_aged_vs_young_ge100_sensitivity")
          if (!is.null(res)) effect_rows[[length(effect_rows) + 1]] <- res
        }
      }

      if (cfg$dataset == "GSE182786") {
        bmi <- comp[comp$eligible_ge50 & comp$bmi_stratum != "Older_high_BMI", , drop = FALSE]
        if (nrow(bmi) > 0) {
          bmi_group <- ifelse(bmi$bmi_stratum == "Older_normal_BMI", "Aged", bmi$group)
          for (denom in c("proportion_within_leydig", "proportion_whole_testis")) {
            res <- summarize_comparison(bmi, cfg$dataset, boundary_name, denom, "older_normal_BMI_vs_young", bmi_group)
            if (!is.null(res)) effect_rows[[length(effect_rows) + 1]] <- res
          }
        }
      }
    }
  }

  rm(obj, md)
  invisible(gc())
}

composition_df <- do.call(rbind, composition_rows)
eligibility_df <- do.call(rbind, eligibility_rows)
effects_df <- if (length(effect_rows) > 0) do.call(rbind, effect_rows) else data.frame()

write.csv(composition_df, file.path(TABLE_DIR, "substate_composition_by_sample.csv"), row.names = FALSE)
write.csv(eligibility_df, file.path(TABLE_DIR, "sample_boundary_eligibility.csv"), row.names = FALSE)
write.csv(effects_df, file.path(TABLE_DIR, "exact_permutation_bootstrap_effects.csv"), row.names = FALSE)

primary_plot <- composition_df[
  composition_df$sample_aware &
    composition_df$eligible_ge50 &
    grepl("^ESB_", composition_df$boundary),
  ,
  drop = FALSE
]

if (nrow(primary_plot) > 0) {
  primary_plot$group <- factor(primary_plot$group, levels = c("Young", "Aged"))
  p <- ggplot(primary_plot, aes(x = group, y = proportion_within_leydig, colour = group)) +
    geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.16) +
    geom_jitter(width = 0.08, height = 0, size = 1.8, alpha = 0.9) +
    facet_grid(dataset ~ state_label, scales = "free_y", space = "free_x") +
    scale_colour_manual(values = c(Young = "#287D8E", Aged = "#C55A3D")) +
    labs(x = NULL, y = "Proportion within Leydig candidate compartment") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
  ggsave(file.path(FIG_DIR, "WP1_ESB_substate_composition_sample_level.pdf"), p, width = 11, height = 7)
  ggsave(file.path(FIG_DIR, "WP1_ESB_substate_composition_sample_level.png"), p, width = 11, height = 7, dpi = 300)
}

if (nrow(effects_df) > 0) {
  forest <- effects_df[
    effects_df$comparison == "all_aged_vs_young" &
      effects_df$denominator == "proportion_within_leydig" &
      grepl("^ESB_", effects_df$boundary),
    ,
    drop = FALSE
  ]
  forest$label <- paste(forest$dataset, forest$state_label, sep = " | ")
  forest$label <- factor(forest$label, levels = rev(unique(forest$label)))
  p2 <- ggplot(forest, aes(x = aged_minus_young, y = label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "#888888") +
    geom_errorbarh(aes(xmin = ci95_low, xmax = ci95_high), height = 0.18, linewidth = 0.55, colour = "#4D4D4D") +
    geom_point(size = 2.2, colour = "#1E5A6D") +
    labs(x = "Aged minus Young proportion", y = NULL) +
    theme_classic(base_size = 10)
  ggsave(file.path(FIG_DIR, "WP1_ESB_substate_composition_effect_forest.pdf"), p2, width = 7.2, height = 5.8)
  ggsave(file.path(FIG_DIR, "WP1_ESB_substate_composition_effect_forest.png"), p2, width = 7.2, height = 5.8, dpi = 300)
}

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "2_sessionInfo.txt"))

