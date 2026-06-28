#!/usr/bin/env Rscript

# Composition-adjusted donor-level sensitivity regression. This script uses
# only the first-step output tables and writes additional tables in the same
# isolated analysis directory.

options(stringsAsFactors = FALSE)

OUT <- path.expand(Sys.getenv(
  "GSE254315_COMPOSITION_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE254315_COMPOSITION_ROOT")
))
TABLE_DIR <- file.path(OUT, "tables")
LOG_DIR <- file.path(OUT, "logs")
DOC_DIR <- file.path(OUT, "docs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)
on.exit(writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "composition_adjusted_regression_sessionInfo.txt")), add = TRUE)

boundaries <- list(
  CMB_1_6_18 = c("1", "6", "18"),
  ESB_1_6_18_7 = c("1", "6", "18", "7")
)

read_csv <- function(name) read.csv(file.path(TABLE_DIR, name), stringsAsFactors = FALSE, check.names = FALSE)
expression <- read_csv("04_donor_cluster_target_logCPM.csv")
composition <- read_csv("02_donor_cluster_composition.csv")
membership <- read_csv("00_readout_gene_membership.csv")

required_expression <- c("sample_id", "group", "cluster", "gene", "raw_count", "library_size")
required_composition <- c("sample_id", "group", "boundary", "cluster", "proportion_within_boundary")
stopifnot(all(required_expression %in% colnames(expression)))
stopifnot(all(required_composition %in% colnames(composition)))

extract_group_term <- function(fit) {
  sm <- summary(fit)$coefficients
  term <- grep("^groupAged$", rownames(sm), value = TRUE)
  if (length(term) != 1L) stop("Could not identify the Aged-vs-Young coefficient.")
  ci <- suppressMessages(confint(fit, parm = term))
  data.frame(
    age_beta = unname(sm[term, "Estimate"]),
    age_se = unname(sm[term, "Std. Error"]),
    age_t = unname(sm[term, "t value"]),
    age_p = unname(sm[term, "Pr(>|t|)"]),
    age_ci95_low = unname(ci[1]),
    age_ci95_high = unname(ci[2]),
    model_r_squared = summary(fit)$r.squared,
    model_adj_r_squared = summary(fit)$adj.r.squared,
    residual_df = df.residual(fit),
    stringsAsFactors = FALSE
  )
}

boundary_score_rows <- list()
for (boundary_name in names(boundaries)) {
  clusters <- boundaries[[boundary_name]]
  sub <- expression[expression$cluster %in% clusters, , drop = FALSE]
  counts_by_gene <- aggregate(raw_count ~ sample_id + group + gene, sub, sum, na.rm = TRUE)
  library_by_cluster <- unique(sub[, c("sample_id", "group", "cluster", "library_size"), drop = FALSE])
  library_by_sample <- aggregate(library_size ~ sample_id + group, library_by_cluster, sum, na.rm = TRUE)
  merged <- merge(counts_by_gene, library_by_sample, by = c("sample_id", "group"), all.x = TRUE)
  merged$logCPM_pc0.5 <- log2(((merged$raw_count + 0.5) / (merged$library_size + 1)) * 1e6)

  for (readout in unique(membership$readout)) {
    genes <- membership$gene[membership$readout == readout & membership$present_in_assay]
    defined <- membership$gene[membership$readout == readout]
    coverage <- length(unique(genes)) / length(unique(defined))
    min_genes <- if (length(unique(defined)) == 1L) 1L else 3L
    use_score <- length(unique(genes)) >= min_genes && coverage >= 0.5
    score_sub <- merged[merged$gene %in% genes, , drop = FALSE]
    score <- aggregate(logCPM_pc0.5 ~ sample_id + group, score_sub, median, na.rm = TRUE)
    names(score)[3] <- "boundary_score_median_logCPM"
    score$boundary <- boundary_name
    score$readout <- readout
    score$n_genes_defined <- length(unique(defined))
    score$n_genes_available <- length(unique(genes))
    score$coverage <- coverage
    score$sufficient_coverage <- use_score
    if (!use_score) score$boundary_score_median_logCPM <- NA_real_
    boundary_score_rows[[length(boundary_score_rows) + 1L]] <- score
  }
}
boundary_scores <- do.call(rbind, boundary_score_rows)
write.csv(boundary_scores, file.path(TABLE_DIR, "08_boundary_pseudobulk_readout_scores.csv"), row.names = FALSE)

model_rows <- list()
for (boundary_name in names(boundaries)) {
  clusters <- boundaries[[boundary_name]]
  comp <- composition[
    composition$boundary == boundary_name & composition$cluster %in% clusters,
    c("sample_id", "group", "cluster", "proportion_within_boundary"), drop = FALSE
  ]
  comp_wide <- reshape(
    comp,
    idvar = c("sample_id", "group"),
    timevar = "cluster",
    direction = "wide"
  )
  names(comp_wide) <- sub("proportion_within_boundary\\.", "p_cluster_", names(comp_wide))
  reference_cluster <- clusters[[1]]
  adjustment_clusters <- setdiff(clusters, reference_cluster)
  adjustment_vars <- paste0("p_cluster_", adjustment_clusters)
  if (!all(adjustment_vars %in% names(comp_wide))) stop("Missing composition covariates for ", boundary_name)

  scores <- boundary_scores[boundary_scores$boundary == boundary_name & boundary_scores$sufficient_coverage, , drop = FALSE]
  for (readout in unique(scores$readout)) {
    df <- merge(
      scores[scores$readout == readout, c("sample_id", "group", "boundary_score_median_logCPM"), drop = FALSE],
      comp_wide,
      by = c("sample_id", "group"),
      all = FALSE
    )
    df <- df[is.finite(df$boundary_score_median_logCPM), , drop = FALSE]
    df$group <- factor(df$group, levels = c("Young", "Aged"))
    if (nrow(df) < length(adjustment_vars) + 6L || any(table(df$group) < 3L)) {
      next
    }
    unadjusted <- lm(boundary_score_median_logCPM ~ group, data = df)
    adjusted <- lm(
      as.formula(paste("boundary_score_median_logCPM ~ group +", paste(adjustment_vars, collapse = " + "))),
      data = df
    )
    unadj_stats <- extract_group_term(unadjusted)
    adj_stats <- extract_group_term(adjusted)
    model_rows[[length(model_rows) + 1L]] <- cbind(
      data.frame(
        boundary = boundary_name,
        readout = readout,
        model = "unadjusted_age_only",
        n_donors = nrow(df),
        n_young = sum(df$group == "Young"),
        n_aged = sum(df$group == "Aged"),
        reference_cluster = reference_cluster,
        adjustment_clusters = paste(adjustment_clusters, collapse = ";"),
        stringsAsFactors = FALSE
      ),
      unadj_stats
    )
    model_rows[[length(model_rows) + 1L]] <- cbind(
      data.frame(
        boundary = boundary_name,
        readout = readout,
        model = "age_plus_cluster_composition",
        n_donors = nrow(df),
        n_young = sum(df$group == "Young"),
        n_aged = sum(df$group == "Aged"),
        reference_cluster = reference_cluster,
        adjustment_clusters = paste(adjustment_clusters, collapse = ";"),
        stringsAsFactors = FALSE
      ),
      adj_stats
    )
  }
}

models <- do.call(rbind, model_rows)
models$age_p_fdr_within_regression_sensitivity <- p.adjust(models$age_p, method = "BH")
write.csv(models, file.path(TABLE_DIR, "09_composition_adjusted_age_effects.csv"), row.names = FALSE)

writeLines(c(
  "# Composition-Adjusted Boundary Regression",
  "",
  "Each model uses all 23 GSE254315 donor/sample units.",
  "The adjusted model includes the donor-level proportions of all non-reference",
  "clusters in the fixed boundary; cluster 1 is the reference component.",
  "This is an exclusionary sensitivity analysis. It can show whether the age",
  "coefficient is strongly attenuated after accounting for cluster proportions,",
  "but it does not distinguish confounding from a composition change caused by age",
  "and does not establish a cell-intrinsic causal mechanism."
), file.path(DOC_DIR, "composition_adjusted_regression_README.md"))

message("Composition-adjusted boundary regression completed: ", Sys.time())
