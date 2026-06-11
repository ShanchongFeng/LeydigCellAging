#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260602)

STATE_ROOT <- path.expand(Sys.getenv(
  "STATE_STRUCTURE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "STATE_STRUCTURE_ROOT")
))
PB_DIR <- file.path(path.expand(Sys.getenv(
  "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT",
  unset = file.path(Sys.getenv("HOME"), "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT")
)), "pseudobulk")
ORTHOLOG_FILE <- Sys.getenv(
  "WP4_ORTHOLOG_FILE",
  unset = file.path(
    path.expand(Sys.getenv("ORTHOLOG_CONSERVATION_ROOT", unset = file.path(Sys.getenv("HOME"), "ORTHOLOG_CONSERVATION_ROOT"))),
    "ortholog_module_gene_map_one_to_one.csv"
  )
)
TABLE_DIR <- Sys.getenv(
  "WP4_TABLE_DIR",
  unset = path.expand(Sys.getenv("PORTABLE_SIGNATURE_AUDIT_ROOT", unset = file.path(Sys.getenv("HOME"), "PORTABLE_SIGNATURE_AUDIT_ROOT")))
)
FIG_DIR <- file.path(STATE_ROOT, "figures", "exploratory")
DOC_DIR <- file.path(STATE_ROOT, "docs")
LOG_DIR <- file.path(STATE_ROOT, "logs")
FIG_PREFIX <- Sys.getenv("WP4_FIG_PREFIX", unset = "WP4_signature_gene_z_effect_forest")
DEFINITION_FILE <- Sys.getenv("WP4_DEFINITION_FILE", unset = "WP4_analysis_definition.md")
SESSION_FILE <- Sys.getenv("WP4_SESSION_FILE", unset = "3_homologene_sessionInfo.txt")
MAPPING_LABEL <- Sys.getenv("WP4_MAPPING_LABEL", unset = "HomoloGene build68 one-to-one")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

n_boot <- 5000
n_perm_random <- 10000

configs <- list(
  list(dataset = "GSE303193", boundary = "CMB_16", species = "mouse"),
  list(dataset = "GSE303193", boundary = "ESB_16_17_19", species = "mouse"),
  list(dataset = "GSE254315", boundary = "CORE_1_6_18", species = "human"),
  list(dataset = "GSE254315", boundary = "EXT_1_6_18_7", species = "human"),
  list(dataset = "GSE182786", boundary = "CORE_0_17", species = "human"),
  list(dataset = "GSE182786", boundary = "EXT_0_14_17", species = "human")
)

ortholog_map <- read.csv(ORTHOLOG_FILE, check.names = FALSE)
as_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}
ortholog_map$mapped_one_to_one <- as_bool(ortholog_map$mapped_one_to_one)
tan_map <- ortholog_map[
  ortholog_map$module == "hdwgcna_tan" &
    ortholog_map$mapped_one_to_one &
    !is.na(ortholog_map$human_symbol),
  ,
  drop = FALSE
]
if (nrow(tan_map) == 0) {
  stop("No hdwgcna_tan one-to-one orthologs were loaded from: ", ORTHOLOG_FILE)
}

compact_metabolic <- c(
  "HMGCS2", "HMGCL", "BDH1", "ACAT1", "ACADM", "ACADL", "ACADVL",
  "CPT1A", "CPT2", "HADH", "HADHA", "PPARA", "PPARGC1A"
)
steroidogenesis <- c(
  "INSL3", "STAR", "CYP11A1", "HSD3B1", "HSD3B2",
  "CYP17A1", "LHCGR", "NR5A1", "SCARB1", "TSPO"
)

species_sets <- function(species) {
  if (species == "human") {
    tan <- unique(tan_map$human_symbol)
    hmgcs2 <- "HMGCS2"
    compact <- compact_metabolic
    steroid <- steroidogenesis
  } else {
    mouse_symbols <- tan_map$mouse_symbol
    if ("mouse_current_symbol" %in% colnames(tan_map)) {
      use_current <- !is.na(tan_map$mouse_current_symbol) &
        nzchar(tan_map$mouse_current_symbol)
      mouse_symbols[use_current] <- tan_map$mouse_current_symbol[use_current]
    }
    tan <- unique(mouse_symbols)
    hmgcs2 <- "Hmgcs2"
    compact <- c(
      "Hmgcs2", "Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl",
      "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a"
    )
    steroid <- c(
      "Insl3", "Star", "Cyp11a1", "Hsd3b1",
      "Cyp17a1", "Lhcgr", "Nr5a1", "Scarb1", "Tspo"
    )
  }
  list(
    hmgcs2_single_gene = hmgcs2,
    tan_ortholog_portable = tan,
    tan_ortholog_portable_leaveout_hmgcs2 = setdiff(tan, hmgcs2),
    compact_metabolic_support = compact,
    compact_metabolic_support_leaveout_hmgcs2 = setdiff(compact, hmgcs2),
    steroidogenesis_comparator = steroid
  )
}

bh_adjust <- function(p) {
  ifelse(is.finite(p), p.adjust(p, method = "BH"), NA_real_)
}

bmi_stratum <- function(dataset, sample_name, group) {
  if (dataset[[1]] != "GSE182786") return(rep("not_applicable", length(sample_name)))
  out <- rep("Young", length(sample_name))
  older_num <- suppressWarnings(as.integer(sub("^Older", "", sample_name)))
  out[group == "Aged" & older_num %in% 1:5] <- "Older_normal_BMI"
  out[group == "Aged" & older_num %in% 6:8] <- "Older_high_BMI"
  out[group == "Aged" & is.na(older_num)] <- "Older_BMI_unknown"
  out
}

resolve_indices <- function(expr, genes) {
  match_upper <- match(toupper(genes), toupper(colnames(expr)))
  match_upper[!is.na(match_upper)]
}

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

score_signatures <- function(expr, dataset, boundary, species) {
  sets <- species_sets(species)
  gene_z <- safe_gene_z(expr)
  percentile <- sample_percentile_matrix(expr)
  score_rows <- list()
  coverage_rows <- list()
  for (signature in names(sets)) {
    genes <- sets[[signature]]
    idx <- resolve_indices(expr, genes)
    present <- colnames(expr)[idx]
    missing <- genes[!toupper(genes) %in% toupper(present)]
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      dataset = dataset,
      boundary = boundary,
      species = species,
      signature = signature,
      n_requested = length(genes),
      n_present = length(idx),
      fraction_present = length(idx) / length(genes),
      present_genes = paste(present, collapse = ";"),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
    method_mats <- list(
      mean_logcpm = expr,
      gene_z_mean = gene_z,
      sample_expression_percentile_mean = percentile
    )
    for (method in names(method_mats)) {
      values <- if (length(idx) == 0) {
        rep(NA_real_, nrow(expr))
      } else {
        rowMeans(method_mats[[method]][, idx, drop = FALSE], na.rm = TRUE)
      }
      score_rows[[length(score_rows) + 1]] <- data.frame(
        dataset = dataset,
        boundary = boundary,
        species = species,
        sample = rownames(expr),
        signature = signature,
        method = method,
        score = values,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    scores = do.call(rbind, score_rows),
    coverage = do.call(rbind, coverage_rows)
  )
}

bootstrap_mean_effect <- function(values, group, n_boot = 5000) {
  young <- values[group == "Young" & is.finite(values)]
  aged <- values[group == "Aged" & is.finite(values)]
  if (length(young) < 2 || length(aged) < 2) return(c(low = NA_real_, high = NA_real_))
  effects <- replicate(n_boot, mean(sample(aged, replace = TRUE)) - mean(sample(young, replace = TRUE)))
  unname(quantile(effects, c(0.025, 0.975), na.rm = TRUE))
}

permutation_mean_effect <- function(values, group, n_random = 10000) {
  keep <- is.finite(values) & group %in% c("Young", "Aged")
  values <- values[keep]
  group <- group[keep]
  idx <- seq_along(values)
  n_young <- sum(group == "Young")
  n_aged <- sum(group == "Aged")
  if (n_young < 2 || n_aged < 2) {
    return(list(effect = NA_real_, p = NA_real_, type = "insufficient_group_size", n = 0, n_possible = NA_real_))
  }
  observed <- mean(values[group == "Aged"]) - mean(values[group == "Young"])
  n_possible <- choose(length(idx), n_young)
  one_effect <- function(young_idx) {
    aged_idx <- setdiff(idx, young_idx)
    mean(values[aged_idx]) - mean(values[young_idx])
  }
  if (n_possible <= 5000) {
    combos <- combn(idx, n_young, simplify = FALSE)
    null <- vapply(combos, one_effect, numeric(1))
    p <- mean(abs(null) >= abs(observed) - 1e-12)
    type <- "exact_enumeration"
  } else {
    null <- replicate(n_random, one_effect(sample(idx, n_young, replace = FALSE)))
    p <- (1 + sum(abs(null) >= abs(observed) - 1e-12)) / (1 + length(null))
    type <- "random_permutation"
  }
  list(effect = observed, p = p, type = type, n = length(null), n_possible = n_possible)
}

summarize_effect <- function(score_df, dataset, boundary, comparison, signature, method) {
  values <- score_df$score
  group <- score_df$analysis_group
  perm <- permutation_mean_effect(values, group, n_perm_random)
  ci <- bootstrap_mean_effect(values, group, n_boot)
  data.frame(
    dataset = dataset,
    boundary = boundary,
    comparison = comparison,
    signature = signature,
    method = method,
    n_young = sum(group == "Young"),
    n_aged = sum(group == "Aged"),
    mean_young = mean(values[group == "Young"], na.rm = TRUE),
    mean_aged = mean(values[group == "Aged"], na.rm = TRUE),
    aged_minus_young = perm$effect,
    bootstrap_ci95_low = ci[[1]],
    bootstrap_ci95_high = ci[[2]],
    permutation_p = perm$p,
    permutation_type = perm$type,
    n_permutations = perm$n,
    n_possible_labelings = perm$n_possible,
    stringsAsFactors = FALSE
  )
}

summarize_loo <- function(score_df, dataset, boundary, comparison, signature, method) {
  full <- mean(score_df$score[score_df$analysis_group == "Aged"], na.rm = TRUE) -
    mean(score_df$score[score_df$analysis_group == "Young"], na.rm = TRUE)
  rows <- lapply(score_df$sample, function(omitted) {
    work <- score_df[score_df$sample != omitted, , drop = FALSE]
    effect <- mean(work$score[work$analysis_group == "Aged"], na.rm = TRUE) -
      mean(work$score[work$analysis_group == "Young"], na.rm = TRUE)
    data.frame(
      dataset = dataset,
      boundary = boundary,
      comparison = comparison,
      signature = signature,
      method = method,
      omitted_sample = omitted,
      full_aged_minus_young = full,
      leave_one_out_aged_minus_young = effect,
      sign_consistent = is.finite(full) && is.finite(effect) && sign(full) == sign(effect),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

all_scores <- list()
all_coverage <- list()
all_effects <- list()
all_loo <- list()

for (cfg in configs) {
  prefix <- paste0(cfg$dataset, "__", cfg$boundary)
  expr <- read.csv(
    file.path(PB_DIR, paste0(prefix, "__sample_logcpm.csv")),
    row.names = 1,
    check.names = FALSE
  )
  meta <- read.csv(
    file.path(PB_DIR, paste0(prefix, "__sample_metadata.csv")),
    check.names = FALSE
  )
  scored <- score_signatures(expr, cfg$dataset, cfg$boundary, cfg$species)
  scores <- merge(scored$scores, meta, by = c("dataset", "boundary", "species", "sample"), all.x = TRUE, sort = FALSE)
  scores$bmi_stratum <- bmi_stratum(scores$dataset, scores$sample, scores$group)
  all_scores[[length(all_scores) + 1]] <- scores
  all_coverage[[length(all_coverage) + 1]] <- scored$coverage

  comparisons <- list(all_aged_vs_young = scores)
  if (cfg$dataset == "GSE182786") {
    comparisons$older_normal_BMI_vs_young <- scores[scores$bmi_stratum != "Older_high_BMI", , drop = FALSE]
  }
  for (comparison in names(comparisons)) {
    work <- comparisons[[comparison]]
    work$analysis_group <- work$group
    if (comparison == "older_normal_BMI_vs_young") {
      work$analysis_group <- ifelse(work$bmi_stratum == "Older_normal_BMI", "Aged", work$group)
    }
    for (signature in unique(work$signature)) {
      for (method in unique(work$method)) {
        sub <- work[work$signature == signature & work$method == method, , drop = FALSE]
        all_effects[[length(all_effects) + 1]] <- summarize_effect(
          sub, cfg$dataset, cfg$boundary, comparison, signature, method
        )
        all_loo[[length(all_loo) + 1]] <- summarize_loo(
          sub, cfg$dataset, cfg$boundary, comparison, signature, method
        )
      }
    }
  }
  rm(expr)
  invisible(gc())
}

score_table <- do.call(rbind, all_scores)
coverage_table <- do.call(rbind, all_coverage)
effect_table <- do.call(rbind, all_effects)
loo_table <- do.call(rbind, all_loo)

family_key <- interaction(
  effect_table$dataset,
  effect_table$boundary,
  effect_table$comparison,
  effect_table$method,
  drop = TRUE
)
effect_table$permutation_fdr_within_test_family <- ave(
  effect_table$permutation_p,
  family_key,
  FUN = bh_adjust
)

loo_key <- interaction(
  loo_table$dataset,
  loo_table$boundary,
  loo_table$comparison,
  loo_table$signature,
  loo_table$method,
  drop = TRUE
)
loo_summary <- do.call(rbind, lapply(split(loo_table, loo_key), function(x) {
  valid <- is.finite(x$leave_one_out_aged_minus_young)
  data.frame(
    dataset = x$dataset[[1]],
    boundary = x$boundary[[1]],
    comparison = x$comparison[[1]],
    signature = x$signature[[1]],
    method = x$method[[1]],
    n_leave_one_out = nrow(x),
    n_valid = sum(valid),
    fraction_sign_consistent = if (sum(valid) > 0) mean(x$sign_consistent[valid]) else NA_real_,
    min_leave_one_out_effect = if (sum(valid) > 0) min(x$leave_one_out_aged_minus_young[valid]) else NA_real_,
    max_leave_one_out_effect = if (sum(valid) > 0) max(x$leave_one_out_aged_minus_young[valid]) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

direction_summary <- do.call(rbind, lapply(
  split(effect_table, interaction(effect_table$signature, effect_table$method, effect_table$comparison, drop = TRUE)),
  function(x) {
    data.frame(
      signature = x$signature[[1]],
      method = x$method[[1]],
      comparison = x$comparison[[1]],
      n_dataset_boundaries = nrow(x),
      n_aged_lower = sum(x$aged_minus_young < 0, na.rm = TRUE),
      fraction_aged_lower = mean(x$aged_minus_young < 0, na.rm = TRUE),
      n_nominal_p_lt_0_05 = sum(x$permutation_p < 0.05, na.rm = TRUE),
      n_fdr_lt_0_10 = sum(x$permutation_fdr_within_test_family < 0.10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

write.csv(coverage_table, file.path(TABLE_DIR, "signature_coverage.csv"), row.names = FALSE)
write.csv(score_table, file.path(TABLE_DIR, "signature_scores_by_sample.csv"), row.names = FALSE)
write.csv(effect_table, file.path(TABLE_DIR, "signature_effects.csv"), row.names = FALSE)
write.csv(loo_table, file.path(TABLE_DIR, "signature_leave_one_out.csv"), row.names = FALSE)
write.csv(loo_summary, file.path(TABLE_DIR, "signature_leave_one_out_summary.csv"), row.names = FALSE)
write.csv(direction_summary, file.path(TABLE_DIR, "signature_direction_summary.csv"), row.names = FALSE)

plot_df <- effect_table[
  effect_table$method == "gene_z_mean" &
    effect_table$comparison == "all_aged_vs_young" &
    effect_table$signature %in% c(
      "hmgcs2_single_gene",
      "tan_ortholog_portable_leaveout_hmgcs2",
      "compact_metabolic_support_leaveout_hmgcs2",
      "steroidogenesis_comparator"
    ),
  ,
  drop = FALSE
]
plot_df$dataset_boundary <- paste(plot_df$dataset, plot_df$boundary, sep = " | ")
plot_df$dataset_boundary <- factor(plot_df$dataset_boundary, levels = rev(unique(plot_df$dataset_boundary)))

p <- ggplot(plot_df, aes(
  x = aged_minus_young,
  y = dataset_boundary,
  color = dataset
)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "#6B7280", linetype = 2) +
  geom_errorbarh(aes(
    xmin = bootstrap_ci95_low,
    xmax = bootstrap_ci95_high
  ), height = 0.16, linewidth = 0.55, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE) +
  facet_wrap(~ signature, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c(
    "GSE303193" = "#4F7D5C",
    "GSE254315" = "#176B87",
    "GSE182786" = "#B8573F"
  )) +
  labs(x = "Aged - Young pseudobulk score", y = NULL, color = NULL) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 7)
  )

ggsave(file.path(FIG_DIR, paste0(FIG_PREFIX, ".pdf")), p, width = 9.5, height = 7.5)
ggsave(file.path(FIG_DIR, paste0(FIG_PREFIX, ".png")), p, width = 9.5, height = 7.5, dpi = 300)

capture.output(sessionInfo(), file = file.path(LOG_DIR, SESSION_FILE))

