#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260605)

STATE_ROOT <- path.expand(Sys.getenv(
  "STATE_STRUCTURE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "STATE_STRUCTURE_ROOT")
))
PB_DIR <- file.path(path.expand(Sys.getenv(
  "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT",
  unset = file.path(Sys.getenv("HOME"), "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT")
)), "pseudobulk")
TARGET_MAP_FILE <- file.path(
  path.expand(Sys.getenv("ENSEMBL_ORTHOLOG_ROOT", unset = file.path(Sys.getenv("HOME"), "ENSEMBL_ORTHOLOG_ROOT"))),
  "ensembl_tan_mouse_human_one_to_one_map.csv"
)
ORTHOLOG_POOL_FILE <- file.path(
  path.expand(Sys.getenv("ORTHOLOG_CONSERVATION_ROOT", unset = file.path(Sys.getenv("HOME"), "ORTHOLOG_CONSERVATION_ROOT"))), "homologene_mouse_human_orthologs_one_to_one.csv"
)
TABLE_DIR <- path.expand(Sys.getenv("RANDOM_SIGNATURE_NULL_ROOT", unset = file.path(Sys.getenv("HOME"), "RANDOM_SIGNATURE_NULL_ROOT")))
FIG_DIR <- file.path(STATE_ROOT, "figures", "exploratory")
DOC_DIR <- file.path(STATE_ROOT, "docs")
LOG_DIR <- file.path(STATE_ROOT, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

n_null <- as.integer(Sys.getenv("WP4_RANDOM_NULL_N", unset = "10000"))
match_k <- as.integer(Sys.getenv("WP4_RANDOM_NULL_MATCH_K", unset = "100"))
if (!is.finite(n_null) || n_null < 1000) stop("WP4_RANDOM_NULL_N must be >=1000")
if (!is.finite(match_k) || match_k < 20) stop("WP4_RANDOM_NULL_MATCH_K must be >=20")

configs <- list(
  list(dataset = "GSE303193", boundary = "CMB_16", species = "mouse"),
  list(dataset = "GSE303193", boundary = "ESB_16_17_19", species = "mouse"),
  list(dataset = "GSE254315", boundary = "CORE_1_6_18", species = "human"),
  list(dataset = "GSE254315", boundary = "EXT_1_6_18_7", species = "human"),
  list(dataset = "GSE182786", boundary = "CORE_0_17", species = "human"),
  list(dataset = "GSE182786", boundary = "EXT_0_14_17", species = "human")
)

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

safe_gene_z <- function(expr) {
  means <- colMeans(expr, na.rm = TRUE)
  sds <- apply(expr, 2, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1
  sweep(sweep(expr, 2, means, "-"), 2, sds, "/")
}

sample_percentile_matrix <- function(expr) {
  out <- t(vapply(seq_len(nrow(expr)), function(i) {
    rank(expr[i, ], ties.method = "average", na.last = "keep") / ncol(expr)
  }, numeric(ncol(expr))))
  dimnames(out) <- dimnames(expr)
  out
}

resolve_gene_indices <- function(expr, genes) {
  match(toupper(genes), toupper(colnames(expr)))
}

effect_by_gene <- function(mat, group) {
  keep <- group %in% c("Young", "Aged")
  mat <- mat[keep, , drop = FALSE]
  group <- group[keep]
  if (sum(group == "Young") < 2 || sum(group == "Aged") < 2) {
    return(rep(NA_real_, ncol(mat)))
  }
  colMeans(mat[group == "Aged", , drop = FALSE], na.rm = TRUE) -
    colMeans(mat[group == "Young", , drop = FALSE], na.rm = TRUE)
}

empirical_p_lower <- function(null, observed) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) == 0) return(NA_real_)
  (1 + sum(null <= observed)) / (1 + length(null))
}

empirical_p_two_sided <- function(null, observed) {
  null <- null[is.finite(null)]
  if (!is.finite(observed) || length(null) == 0) return(NA_real_)
  (1 + sum(abs(null) >= abs(observed))) / (1 + length(null))
}

matrix_row_means_from_sets <- function(values, set_index) {
  selected <- matrix(
    values[set_index],
    nrow = nrow(set_index),
    ncol = ncol(set_index)
  )
  rowMeans(selected, na.rm = TRUE)
}

target_map <- read.csv(TARGET_MAP_FILE, check.names = FALSE)
target_map$mapped_one_to_one <- as_bool(target_map$mapped_one_to_one)
target_map <- target_map[
  target_map$module == "hdwgcna_tan" &
    target_map$mapped_one_to_one &
    !is.na(target_map$human_symbol) &
    toupper(target_map$human_symbol) != "HMGCS2",
  ,
  drop = FALSE
]
target_map$mouse_gene <- target_map$mouse_symbol
use_current <- !is.na(target_map$mouse_current_symbol) &
  nzchar(target_map$mouse_current_symbol)
target_map$mouse_gene[use_current] <- target_map$mouse_current_symbol[use_current]
target_map <- unique(target_map[, c(
  "mouse_gene", "human_symbol", "mouse_ensembl_gene_id",
  "human_ensembl_gene_id", "orthology_type"
)])
names(target_map)[names(target_map) == "human_symbol"] <- "human_gene"
if (nrow(target_map) < 10) {
  stop(
    "Too few Ensembl one-to-one genes remain after HMGCS2 removal: ",
    nrow(target_map)
  )
}

ortholog_pool <- read.csv(ORTHOLOG_POOL_FILE, check.names = FALSE)
ortholog_pool <- ortholog_pool[
  as_bool(ortholog_pool$one_to_one) &
    !is.na(ortholog_pool$mouse_symbol) &
    !is.na(ortholog_pool$human_symbol),
  c("homologene_id", "mouse_symbol", "human_symbol"),
  drop = FALSE
]
names(ortholog_pool)[names(ortholog_pool) == "mouse_symbol"] <- "mouse_gene"
names(ortholog_pool)[names(ortholog_pool) == "human_symbol"] <- "human_gene"
ortholog_pool <- ortholog_pool[
  !duplicated(toupper(ortholog_pool$mouse_gene)) &
    !duplicated(toupper(ortholog_pool$human_gene)),
  ,
  drop = FALSE
]

loaded <- list()
for (cfg in configs) {
  prefix <- paste0(cfg$dataset, "__", cfg$boundary)
  expr <- as.matrix(read.csv(
    file.path(PB_DIR, paste0(prefix, "__sample_logcpm.csv")),
    row.names = 1,
    check.names = FALSE
  ))
  storage.mode(expr) <- "double"
  meta <- read.csv(
    file.path(PB_DIR, paste0(prefix, "__sample_metadata.csv")),
    check.names = FALSE
  )
  meta <- meta[match(rownames(expr), meta$sample), , drop = FALSE]
  if (any(is.na(meta$sample))) stop("Metadata/sample mismatch for ", prefix)
  loaded[[prefix]] <- list(cfg = cfg, expr = expr, meta = meta)
}

pool_keep <- rep(TRUE, nrow(ortholog_pool))
for (entry in loaded) {
  genes <- if (entry$cfg$species == "mouse") {
    ortholog_pool$mouse_gene
  } else {
    ortholog_pool$human_gene
  }
  pool_keep <- pool_keep & !is.na(resolve_gene_indices(entry$expr, genes))
}
ortholog_pool <- ortholog_pool[pool_keep, , drop = FALSE]

target_mouse <- toupper(target_map$mouse_gene)
target_human <- toupper(target_map$human_gene)
ortholog_pool <- ortholog_pool[
  !toupper(ortholog_pool$mouse_gene) %in% c(target_mouse, "HMGCS2") &
    !toupper(ortholog_pool$human_gene) %in% c(target_human, "HMGCS2"),
  ,
  drop = FALSE
]
if (nrow(ortholog_pool) < 1000) {
  stop("Too few common expressed one-to-one ortholog pairs: ", nrow(ortholog_pool))
}

for (entry in loaded) {
  target_genes <- if (entry$cfg$species == "mouse") {
    target_map$mouse_gene
  } else {
    target_map$human_gene
  }
  if (any(is.na(resolve_gene_indices(entry$expr, target_genes)))) {
    stop("Target signature is not fully represented in ", entry$cfg$dataset, " ", entry$cfg$boundary)
  }
}

profile_pool <- matrix(
  NA_real_,
  nrow = nrow(ortholog_pool),
  ncol = length(loaded),
  dimnames = list(NULL, names(loaded))
)
profile_target <- matrix(
  NA_real_,
  nrow = nrow(target_map),
  ncol = length(loaded),
  dimnames = list(target_map$human_gene, names(loaded))
)

for (j in seq_along(loaded)) {
  entry <- loaded[[j]]
  pool_genes <- if (entry$cfg$species == "mouse") {
    ortholog_pool$mouse_gene
  } else {
    ortholog_pool$human_gene
  }
  target_genes <- if (entry$cfg$species == "mouse") {
    target_map$mouse_gene
  } else {
    target_map$human_gene
  }
  pool_mean <- colMeans(
    entry$expr[, resolve_gene_indices(entry$expr, pool_genes), drop = FALSE],
    na.rm = TRUE
  )
  target_mean <- colMeans(
    entry$expr[, resolve_gene_indices(entry$expr, target_genes), drop = FALSE],
    na.rm = TRUE
  )
  profile_pool[, j] <- rank(pool_mean, ties.method = "average") / length(pool_mean)
  profile_target[, j] <- vapply(
    target_mean,
    function(x) mean(pool_mean <= x),
    numeric(1)
  )
}

nearest_candidates <- lapply(seq_len(nrow(target_map)), function(i) {
  distance <- rowSums((profile_pool -
    matrix(profile_target[i, ], nrow(profile_pool), ncol(profile_pool), byrow = TRUE))^2)
  order(distance)[seq_len(min(match_k, length(distance)))]
})

unmatched_sets <- t(replicate(
  n_null,
  sample.int(nrow(ortholog_pool), nrow(target_map), replace = FALSE)
))
matched_sets <- matrix(NA_integer_, nrow = n_null, ncol = nrow(target_map))
for (b in seq_len(n_null)) {
  used <- integer()
  target_order <- sample.int(nrow(target_map))
  for (target_i in target_order) {
    candidates <- setdiff(nearest_candidates[[target_i]], used)
    if (length(candidates) == 0) {
      candidates <- setdiff(seq_len(nrow(ortholog_pool)), used)
    }
    selected <- sample(candidates, 1)
    matched_sets[b, target_i] <- selected
    used <- c(used, selected)
  }
}

saveRDS(
  list(
    seed = 20260605,
    n_null = n_null,
    match_k = match_k,
    unmatched_set_index = unmatched_sets,
    expression_matched_set_index = matched_sets
  ),
  file.path(TABLE_DIR, "random_signature_set_indices.rds")
)
write.csv(target_map, file.path(TABLE_DIR, "target_signature_ensembl_release115.csv"), row.names = FALSE)
write.csv(ortholog_pool, file.path(TABLE_DIR, "random_signature_sampling_pool.csv"), row.names = FALSE)

membership_file <- gzfile(
  file.path(TABLE_DIR, "random_signature_set_membership.csv.gz"),
  open = "wt"
)
writeLines("null_type,replicate,position,mouse_gene,human_gene", membership_file)
for (null_type in c("unmatched", "expression_matched")) {
  set_index <- if (null_type == "unmatched") unmatched_sets else matched_sets
  for (position in seq_len(ncol(set_index))) {
    rows <- data.frame(
      null_type = null_type,
      replicate = seq_len(nrow(set_index)),
      position = position,
      mouse_gene = ortholog_pool$mouse_gene[set_index[, position]],
      human_gene = ortholog_pool$human_gene[set_index[, position]]
    )
    write.table(
      rows, membership_file, sep = ",", row.names = FALSE,
      col.names = FALSE, quote = TRUE
    )
  }
}
close(membership_file)

effect_summary <- list()
null_effects <- list()
matching_diagnostics <- list()
observed_store <- list()
null_store <- list()

for (prefix in names(loaded)) {
  entry <- loaded[[prefix]]
  expr <- entry$expr
  meta <- entry$meta
  pool_genes <- if (entry$cfg$species == "mouse") {
    ortholog_pool$mouse_gene
  } else {
    ortholog_pool$human_gene
  }
  target_genes <- if (entry$cfg$species == "mouse") {
    target_map$mouse_gene
  } else {
    target_map$human_gene
  }
  pool_idx <- resolve_gene_indices(expr, pool_genes)
  target_idx <- resolve_gene_indices(expr, target_genes)

  method_mats <- list(
    mean_logcpm = expr,
    gene_z_mean = safe_gene_z(expr),
    sample_expression_percentile_mean = sample_percentile_matrix(expr)
  )
  comparisons <- list(all_aged_vs_young = seq_len(nrow(expr)))
  if (entry$cfg$dataset == "GSE182786") {
    older_num <- suppressWarnings(as.integer(sub("^Older", "", meta$sample)))
    comparisons$older_normal_BMI_vs_young <- which(
      meta$group == "Young" |
        (meta$group == "Aged" & older_num %in% 1:5)
    )
  }

  for (method in names(method_mats)) {
    mat <- method_mats[[method]]
    for (comparison in names(comparisons)) {
      rows <- comparisons[[comparison]]
      group <- meta$group[rows]
      gene_effect <- effect_by_gene(mat[rows, , drop = FALSE], group)
      target_effects <- gene_effect[target_idx]
      pool_effects <- gene_effect[pool_idx]
      observed <- mean(target_effects, na.rm = TRUE)

      for (null_type in c("unmatched", "expression_matched")) {
        set_index <- if (null_type == "unmatched") unmatched_sets else matched_sets
        null <- matrix_row_means_from_sets(pool_effects, set_index)
        key <- paste(prefix, comparison, method, null_type, sep = "||")
        observed_store[[key]] <- observed
        null_store[[key]] <- null

        effect_summary[[length(effect_summary) + 1]] <- data.frame(
          dataset = entry$cfg$dataset,
          boundary = entry$cfg$boundary,
          species = entry$cfg$species,
          comparison = comparison,
          method = method,
          null_type = null_type,
          signature_size = nrow(target_map),
          random_pool_size = nrow(ortholog_pool),
          n_null = length(null),
          observed_aged_minus_young = observed,
          null_mean = mean(null),
          null_sd = sd(null),
          null_q025 = unname(quantile(null, 0.025)),
          null_median = median(null),
          null_q975 = unname(quantile(null, 0.975)),
          empirical_p_aged_lower = empirical_p_lower(null, observed),
          empirical_p_two_sided = empirical_p_two_sided(null, observed),
          null_percentile = mean(null <= observed),
          stringsAsFactors = FALSE
        )
        null_effects[[length(null_effects) + 1]] <- data.frame(
          dataset = entry$cfg$dataset,
          boundary = entry$cfg$boundary,
          comparison = comparison,
          method = method,
          null_type = null_type,
          replicate = seq_along(null),
          random_aged_minus_young = null,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  target_profile_mean <- mean(profile_target[, prefix])
  for (null_type in c("unmatched", "expression_matched")) {
    set_index <- if (null_type == "unmatched") unmatched_sets else matched_sets
    null_profile <- matrix_row_means_from_sets(profile_pool[, prefix], set_index)
    matching_diagnostics[[length(matching_diagnostics) + 1]] <- data.frame(
      dataset = entry$cfg$dataset,
      boundary = entry$cfg$boundary,
      null_type = null_type,
      target_mean_expression_percentile = target_profile_mean,
      null_mean_expression_percentile = mean(null_profile),
      null_q025_expression_percentile = unname(quantile(null_profile, 0.025)),
      null_q975_expression_percentile = unname(quantile(null_profile, 0.975)),
      absolute_target_minus_null_mean = abs(target_profile_mean - mean(null_profile)),
      stringsAsFactors = FALSE
    )
  }
}

effect_summary <- do.call(rbind, effect_summary)
null_effects <- do.call(rbind, null_effects)
matching_diagnostics <- do.call(rbind, matching_diagnostics)
effect_family <- interaction(
  effect_summary$method,
  effect_summary$null_type,
  drop = TRUE
)
effect_summary$empirical_fdr_aged_lower_within_method <- ave(
  effect_summary$empirical_p_aged_lower,
  effect_family,
  FUN = function(x) p.adjust(x, method = "BH")
)
effect_summary$empirical_fdr_two_sided_within_method <- ave(
  effect_summary$empirical_p_two_sided,
  effect_family,
  FUN = function(x) p.adjust(x, method = "BH")
)

cross_summary <- list()
cross_null_plot <- list()
scope_defs <- list(
  six_boundary_sensitivity = c(
    "GSE303193__CMB_16", "GSE303193__ESB_16_17_19",
    "GSE254315__CORE_1_6_18", "GSE254315__EXT_1_6_18_7",
    "GSE182786__CORE_0_17", "GSE182786__EXT_0_14_17"
  ),
  three_canonical_boundaries = c(
    "GSE303193__CMB_16",
    "GSE254315__CORE_1_6_18",
    "GSE182786__CORE_0_17"
  )
)

for (scope in names(scope_defs)) {
  prefixes <- scope_defs[[scope]]
  for (method in c(
    "mean_logcpm", "gene_z_mean", "sample_expression_percentile_mean"
  )) {
    for (null_type in c("unmatched", "expression_matched")) {
      keys <- paste(
        prefixes, "all_aged_vs_young", method, null_type,
        sep = "||"
      )
      observed_vec <- vapply(observed_store[keys], identity, numeric(1))
      null_matrix <- do.call(cbind, null_store[keys])
      observed_n_lower <- sum(observed_vec < 0)
      observed_mean <- mean(observed_vec)
      null_n_lower <- rowSums(null_matrix < 0)
      null_mean <- rowMeans(null_matrix)

      cross_summary[[length(cross_summary) + 1]] <- data.frame(
        scope = scope,
        method = method,
        null_type = null_type,
        n_comparisons = length(prefixes),
        observed_n_aged_lower = observed_n_lower,
        observed_mean_effect = observed_mean,
        null_mean_n_aged_lower = mean(null_n_lower),
        null_q025_n_aged_lower = unname(quantile(null_n_lower, 0.025)),
        null_q975_n_aged_lower = unname(quantile(null_n_lower, 0.975)),
        null_mean_of_mean_effect = mean(null_mean),
        null_q025_mean_effect = unname(quantile(null_mean, 0.025)),
        null_q975_mean_effect = unname(quantile(null_mean, 0.975)),
        empirical_p_direction_count = (
          1 + sum(null_n_lower >= observed_n_lower)
        ) / (1 + length(null_n_lower)),
        empirical_p_mean_aged_lower = empirical_p_lower(
          null_mean, observed_mean
        ),
        empirical_p_joint = (
          1 + sum(
            null_n_lower >= observed_n_lower &
              null_mean <= observed_mean
          )
        ) / (1 + length(null_n_lower)),
        stringsAsFactors = FALSE
      )
      if (scope == "six_boundary_sensitivity") {
        cross_null_plot[[length(cross_null_plot) + 1]] <- data.frame(
          method = method,
          null_type = null_type,
          replicate = seq_along(null_mean),
          null_mean_effect = null_mean,
          observed_mean_effect = observed_mean,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
cross_summary <- do.call(rbind, cross_summary)
cross_null_plot <- do.call(rbind, cross_null_plot)
cross_family <- interaction(
  cross_summary$scope,
  cross_summary$null_type,
  drop = TRUE
)
cross_summary$empirical_fdr_direction_count_across_methods <- ave(
  cross_summary$empirical_p_direction_count,
  cross_family,
  FUN = function(x) p.adjust(x, method = "BH")
)
cross_summary$empirical_fdr_mean_aged_lower_across_methods <- ave(
  cross_summary$empirical_p_mean_aged_lower,
  cross_family,
  FUN = function(x) p.adjust(x, method = "BH")
)
cross_summary$empirical_fdr_joint_across_methods <- ave(
  cross_summary$empirical_p_joint,
  cross_family,
  FUN = function(x) p.adjust(x, method = "BH")
)

write.csv(
  effect_summary,
  file.path(TABLE_DIR, "null_effect_summary.csv"),
  row.names = FALSE
)
write.csv(
  cross_summary,
  file.path(TABLE_DIR, "null_cross_cohort_summary.csv"),
  row.names = FALSE
)
write.csv(
  matching_diagnostics,
  file.path(TABLE_DIR, "random_signature_expression_matching_diagnostics.csv"),
  row.names = FALSE
)
null_connection <- gzfile(
  file.path(TABLE_DIR, "null_effects.csv.gz"),
  open = "wt"
)
write.csv(null_effects, null_connection, row.names = FALSE)
close(null_connection)

plot_df <- cross_null_plot[
  cross_null_plot$null_type == "expression_matched",
  ,
  drop = FALSE
]
method_labels <- c(
  mean_logcpm = "Mean logCPM",
  gene_z_mean = "Gene-wise z mean",
  sample_expression_percentile_mean = "Expression percentile"
)
plot_df$method_label <- factor(
  method_labels[plot_df$method],
  levels = unname(method_labels)
)
observed_lines <- unique(plot_df[, c("method_label", "observed_mean_effect")])

p_global <- ggplot(plot_df, aes(x = null_mean_effect)) +
  geom_histogram(
    bins = 60, fill = "#BFD7F2", color = "white", linewidth = 0.15
  ) +
  geom_vline(
    data = observed_lines,
    aes(xintercept = observed_mean_effect),
    color = "#DC879F", linewidth = 0.8
  ) +
  facet_wrap(~ method_label, scales = "free", ncol = 3) +
  labs(
    x = "Mean Aged - Young effect across six boundary-level comparisons",
    y = "Expression-matched random signatures"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#F5F1EA", color = "#B8B2AA"),
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(FIG_DIR, "WP4_size_matched_random_signature_global_null.pdf"),
  p_global, width = 9.5, height = 3.8
)
ggsave(
  file.path(FIG_DIR, "WP4_size_matched_random_signature_global_null.png"),
  p_global, width = 9.5, height = 3.8, dpi = 300
)

forest_df <- effect_summary[
  effect_summary$comparison == "all_aged_vs_young" &
    effect_summary$method == "gene_z_mean" &
    effect_summary$null_type == "expression_matched",
  ,
  drop = FALSE
]
forest_df$dataset_boundary <- paste(
  forest_df$dataset, forest_df$boundary, sep = " | "
)
forest_df$dataset_boundary <- factor(
  forest_df$dataset_boundary,
  levels = rev(forest_df$dataset_boundary)
)
p_forest <- ggplot(
  forest_df,
  aes(x = observed_aged_minus_young, y = dataset_boundary)
) +
  geom_vline(
    xintercept = 0, color = "#6B7280", linewidth = 0.35, linetype = 2
  ) +
  geom_errorbarh(
    aes(xmin = null_q025, xmax = null_q975),
    height = 0.18, color = "#84B1DA", linewidth = 0.7
  ) +
  geom_point(size = 2.2, color = "#DC879F") +
  labs(
    x = "Observed gene-z effect; blue interval = matched random-set 95% null",
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(FIG_DIR, "WP4_size_matched_random_signature_by_comparison.pdf"),
  p_forest, width = 6.2, height = 4.5
)
ggsave(
  file.path(FIG_DIR, "WP4_size_matched_random_signature_by_comparison.png"),
  p_forest, width = 6.2, height = 4.5, dpi = 300
)

capture.output(
  sessionInfo(),
  file = file.path(LOG_DIR, "2_random_null_sessionInfo.txt")
)

