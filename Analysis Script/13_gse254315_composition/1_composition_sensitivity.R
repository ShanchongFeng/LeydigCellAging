#!/usr/bin/env Rscript

# Donor-aware composition and within-cluster expression sensitivity analysis.
# This script writes only under GSE254315_COMPOSITION_ROOT and never changes
# locked inputs.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20260621)

LOCKED_ROOT <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
FUNCTIONAL_ROOT <- path.expand(Sys.getenv(
  "FUNCTIONAL_CATEGORY_ROOT",
  unset = file.path(Sys.getenv("HOME"), "FUNCTIONAL_CATEGORY_ROOT")
))
OUT <- path.expand(Sys.getenv(
  "GSE254315_COMPOSITION_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE254315_COMPOSITION_ROOT")
))

DIRS <- list(
  tables = file.path(OUT, "tables"),
  logs = file.path(OUT, "logs"),
  docs = file.path(OUT, "docs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))
on.exit(writeLines(capture.output(sessionInfo()), file.path(DIRS$logs, "sessionInfo.txt")), add = TRUE)

INPUT_RDS <- path.expand(Sys.getenv(
  "GSE254315_LOCKED_RDS",
  unset = file.path(LOCKED_ROOT, "rds/primary_locked/GSE254315/GSE254315_check__results__03_filtered_clustered_human_seurat.rds")
))
CURATED_SETS <- path.expand(Sys.getenv(
  "FUNCTIONAL_CURATED_SETS_FILE",
  unset = file.path(FUNCTIONAL_ROOT, "tables/curated_gene_effect_sets_used_20260619.csv")
))
GO_SETS <- path.expand(Sys.getenv(
  "FUNCTIONAL_GO_SETS_FILE",
  unset = file.path(FUNCTIONAL_ROOT, "tables/go_external_gene_sets_20260619.csv")
))

stopifnot(file.exists(INPUT_RDS), file.exists(CURATED_SETS), file.exists(GO_SETS))

boundaries <- list(
  CMB_1_6_18 = c("1", "6", "18"),
  ESB_1_6_18_7 = c("1", "6", "18", "7")
)
cell_thresholds <- c(1L, 25L, 50L, 100L)
primary_cell_threshold <- 25L
bootstrap_iterations <- 5000L
random_permutations <- 100000L
exact_permutation_limit <- 250000L

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

get_counts <- function(obj) {
  if ("RNA" %in% names(obj@assays) && inherits(obj[["RNA"]], "Assay5")) {
    count_layers <- grep("^counts", Layers(obj[["RNA"]]), value = TRUE)
    if (length(count_layers) > 1L) obj <- JoinLayers(obj, assay = "RNA")
  }
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  list(object = obj, counts = as(counts, "dgCMatrix"))
}

mean_difference <- function(values, groups) {
  mean(values[groups == "Aged"]) - mean(values[groups == "Young"])
}

bootstrap_ci <- function(values, groups, n_boot = bootstrap_iterations) {
  values <- as.numeric(values)
  groups <- as.character(groups)
  young <- values[groups == "Young"]
  aged <- values[groups == "Aged"]
  if (length(young) < 2L || length(aged) < 2L) return(c(NA_real_, NA_real_))
  draws <- replicate(n_boot, mean(sample(aged, replace = TRUE)) - mean(sample(young, replace = TRUE)))
  unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE))
}

permutation_test <- local({
  combo_cache <- new.env(parent = emptyenv())

  function(values, groups) {
    keep <- is.finite(values) & groups %in% c("Young", "Aged")
    values <- as.numeric(values[keep])
    groups <- as.character(groups[keep])
    n_young <- sum(groups == "Young")
    n_aged <- sum(groups == "Aged")
    n_total <- length(values)
    if (n_young < 2L || n_aged < 2L) {
      return(list(p_value = NA_real_, method = "not_testable", n_labelings = NA_real_))
    }

    observed <- mean_difference(values, groups)
    n_labelings <- choose(n_total, n_young)
    key <- paste(n_total, n_young, sep = "_")
    if (!exists(key, envir = combo_cache, inherits = FALSE)) {
      if (is.finite(n_labelings) && n_labelings <= exact_permutation_limit) {
        cached <- list(
          young_index = combn(seq_len(n_total), n_young),
          method = "exact_label_permutation"
        )
      } else {
        # The same fixed-seed random labelings are reused across readouts with
        # the same donor count and group split. This is both faster and makes
        # sensitivity-grid p values directly comparable.
        cached <- list(
          young_index = replicate(
            random_permutations,
            sample.int(n_total, n_young, replace = FALSE)
          ),
          method = paste0("Monte_Carlo_label_permutation_", random_permutations)
        )
      }
      assign(key, cached, envir = combo_cache)
    }
    cached <- get(key, envir = combo_cache, inherits = FALSE)
    young_sum <- colSums(matrix(values[cached$young_index], nrow = n_young))
    perm_effect <- (sum(values) - young_sum) / n_aged - young_sum / n_young
    list(
      p_value = mean(abs(perm_effect) >= abs(observed) - 1e-12),
      method = cached$method,
      n_labelings = n_labelings
    )
  }
})

summarize_group_effect <- function(values, groups) {
  keep <- is.finite(values) & groups %in% c("Young", "Aged")
  values <- as.numeric(values[keep])
  groups <- as.character(groups[keep])
  n_young <- sum(groups == "Young")
  n_aged <- sum(groups == "Aged")
  if (n_young < 2L || n_aged < 2L) {
    return(data.frame(
      n_young = n_young, n_aged = n_aged,
      mean_young = NA_real_, mean_aged = NA_real_,
      aged_minus_young = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_,
      permutation_p = NA_real_, permutation_method = "not_testable", n_labelings = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  perm <- permutation_test(values, groups)
  ci <- bootstrap_ci(values, groups)
  data.frame(
    n_young = n_young,
    n_aged = n_aged,
    mean_young = mean(values[groups == "Young"]),
    mean_aged = mean(values[groups == "Aged"]),
    aged_minus_young = mean_difference(values, groups),
    ci95_low = ci[[1]],
    ci95_high = ci[[2]],
    permutation_p = perm$p_value,
    permutation_method = perm$method,
    n_labelings = perm$n_labelings,
    stringsAsFactors = FALSE
  )
}

summarize_descriptive_effect <- function(values, groups) {
  keep <- is.finite(values) & groups %in% c("Young", "Aged")
  values <- as.numeric(values[keep])
  groups <- as.character(groups[keep])
  n_young <- sum(groups == "Young")
  n_aged <- sum(groups == "Aged")
  if (n_young < 2L || n_aged < 2L) {
    return(data.frame(
      n_young = n_young, n_aged = n_aged,
      mean_young = NA_real_, mean_aged = NA_real_,
      aged_minus_young = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_,
      permutation_p = NA_real_, permutation_method = "descriptive_sensitivity", n_labelings = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    n_young = n_young,
    n_aged = n_aged,
    mean_young = mean(values[groups == "Young"]),
    mean_aged = mean(values[groups == "Aged"]),
    aged_minus_young = mean_difference(values, groups),
    ci95_low = NA_real_,
    ci95_high = NA_real_,
    permutation_p = NA_real_,
    permutation_method = "descriptive_sensitivity",
    n_labelings = NA_real_,
    stringsAsFactors = FALSE
  )
}

obj <- readRDS(INPUT_RDS)
got <- get_counts(obj)
obj <- got$object
counts <- got$counts
metadata <- as.data.frame(obj@meta.data)
metadata$cell_id <- rownames(metadata)

required_columns <- c("sample_id", "group", "seurat_clusters")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
}
metadata$sample_id <- as.character(metadata$sample_id)
metadata$group <- as.character(metadata$group)
metadata$cluster <- as.character(metadata$seurat_clusters)

sample_group_n <- aggregate(group ~ sample_id, metadata, function(x) length(unique(x)))
if (any(sample_group_n$group != 1L)) stop("At least one sample_id maps to multiple groups.")
sample_metadata <- unique(metadata[, c("sample_id", "group"), drop = FALSE])
sample_metadata <- sample_metadata[order(sample_metadata$sample_id), , drop = FALSE]
if (!all(sample_metadata$group %in% c("Young", "Aged"))) {
  stop("Unexpected group labels: ", paste(sort(unique(sample_metadata$group)), collapse = ", "))
}
if (any(table(sample_metadata$group) < 2L)) stop("Both groups need at least two donor/sample units.")

curated <- read_csv(CURATED_SETS)
go_sets <- read_csv(GO_SETS)
curated_sets <- split(curated$human_gene, curated$set)
go_human <- go_sets[go_sets$species == "human", c("category", "gene"), drop = FALSE]
go_readouts <- split(go_human$gene, go_human$category)

selected_readouts <- list(
  HMGCS2 = "HMGCS2",
  hmgcs2_axis_no_hmgcs2 = curated_sets[["hmgcs2_axis_no_hmgcs2"]],
  curated_fatty_acid_oxidation = curated_sets[["curated_fatty_acid_oxidation_score"]],
  curated_steroidogenesis = curated_sets[["curated_steroidogenesis_score"]],
  go_fatty_acid_beta_oxidation = go_readouts[["go_fatty_acid_beta_oxidation"]],
  go_ketone_body_metabolism = go_readouts[["go_ketone_body_metabolism"]],
  go_mitochondrial_energy = go_readouts[["go_mitochondrial_energy"]],
  go_steroidogenic_execution = go_readouts[["go_steroidogenic_execution"]]
)
selected_readouts <- lapply(selected_readouts, function(x) sort(unique(as.character(x))))
if (any(vapply(selected_readouts, length, integer(1)) == 0L)) stop("One or more frozen readouts are empty.")

target_genes <- sort(unique(unlist(selected_readouts, use.names = FALSE)))
target_genes <- intersect(target_genes, rownames(counts))
if (length(target_genes) == 0L) stop("No selected readout genes are present in the RNA assay.")
target_counts <- counts[target_genes, , drop = FALSE]

membership <- do.call(rbind, lapply(names(selected_readouts), function(readout) {
  genes <- selected_readouts[[readout]]
  data.frame(
    readout = readout,
    gene = genes,
    present_in_assay = genes %in% rownames(counts),
    stringsAsFactors = FALSE
  )
}))
write.csv(membership, file.path(DIRS$tables, "00_readout_gene_membership.csv"), row.names = FALSE)

analysis_plan <- data.frame(
  parameter = c(
    "input_rds", "boundaries", "cell_threshold_grid", "primary_cell_threshold", "readout_score",
    "composition_estimator", "standardized_estimator", "permutation_rule",
    "bootstrap_iterations", "random_seed"
  ),
  value = c(
    INPUT_RDS,
    "CMB=1+6+18; ESB=1+6+18+7",
    paste(cell_thresholds, collapse = ";"),
    as.character(primary_cell_threshold),
    "median log2(CPM+0.5) across available genes in each frozen readout",
    "within-boundary donor cluster proportion",
    "equal mean of donor-specific cluster readout scores across fixed clusters",
    paste0("exact label permutation if <=", exact_permutation_limit, " labelings; otherwise ", random_permutations, " fixed-seed Monte Carlo permutations"),
    as.character(bootstrap_iterations),
    "20260621"
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_plan, file.path(DIRS$tables, "00_analysis_plan.csv"), row.names = FALSE)

all_clusters <- sort(unique(unlist(boundaries, use.names = FALSE)), method = "radix")
cell_count_grid <- expand.grid(
  sample_id = sample_metadata$sample_id,
  cluster = all_clusters,
  stringsAsFactors = FALSE
)
cell_count_grid$group <- sample_metadata$group[match(cell_count_grid$sample_id, sample_metadata$sample_id)]
observed_counts <- as.data.frame(table(
  factor(metadata$sample_id, levels = sample_metadata$sample_id),
  factor(metadata$cluster, levels = all_clusters)
), stringsAsFactors = FALSE)
colnames(observed_counts) <- c("sample_id", "cluster", "n_cells")
cell_count_grid <- merge(cell_count_grid, observed_counts, by = c("sample_id", "cluster"), all.x = TRUE, sort = FALSE)
cell_count_grid$n_cells[is.na(cell_count_grid$n_cells)] <- 0L
cell_count_grid <- cell_count_grid[order(cell_count_grid$sample_id, cell_count_grid$cluster), , drop = FALSE]
write.csv(cell_count_grid, file.path(DIRS$tables, "01_donor_cluster_cell_counts.csv"), row.names = FALSE)

composition_rows <- list()
composition_effect_rows <- list()
for (boundary_name in names(boundaries)) {
  boundary_clusters <- boundaries[[boundary_name]]
  sub <- cell_count_grid[cell_count_grid$cluster %in% boundary_clusters, , drop = FALSE]
  boundary_total <- aggregate(n_cells ~ sample_id, sub, sum)
  colnames(boundary_total)[2] <- "n_cells_boundary"
  sub <- merge(sub, boundary_total, by = "sample_id", all.x = TRUE, sort = FALSE)
  sub$boundary <- boundary_name
  sub$proportion_within_boundary <- ifelse(sub$n_cells_boundary > 0, sub$n_cells / sub$n_cells_boundary, NA_real_)
  composition_rows[[boundary_name]] <- sub

  for (cluster_id in boundary_clusters) {
    cluster_sub <- sub[sub$cluster == cluster_id, , drop = FALSE]
    effect <- summarize_group_effect(cluster_sub$proportion_within_boundary, cluster_sub$group)
    composition_effect_rows[[length(composition_effect_rows) + 1L]] <- cbind(
      data.frame(boundary = boundary_name, cluster = cluster_id, stringsAsFactors = FALSE), effect
    )
  }
}
composition <- do.call(rbind, composition_rows)
composition_effects <- do.call(rbind, composition_effect_rows)
composition_effects$permutation_fdr_within_composition <- p.adjust(composition_effects$permutation_p, method = "BH")
write.csv(composition, file.path(DIRS$tables, "02_donor_cluster_composition.csv"), row.names = FALSE)
write.csv(composition_effects, file.path(DIRS$tables, "03_donor_cluster_composition_effects.csv"), row.names = FALSE)

cell_lookup <- split(metadata$cell_id, paste(metadata$sample_id, metadata$cluster, sep = "||"))
expression_rows <- list()
for (i in seq_len(nrow(cell_count_grid))) {
  sample_id <- cell_count_grid$sample_id[[i]]
  cluster_id <- cell_count_grid$cluster[[i]]
  key <- paste(sample_id, cluster_id, sep = "||")
  cells <- cell_lookup[[key]]
  n_cells <- cell_count_grid$n_cells[[i]]
  if (is.null(cells) || length(cells) == 0L) {
    expression_rows[[i]] <- data.frame(
      sample_id = sample_id,
      group = cell_count_grid$group[[i]],
      cluster = cluster_id,
      n_cells = n_cells,
      gene = target_genes,
      raw_count = NA_real_,
      library_size = NA_real_,
      logCPM_pc0.5 = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  raw_count <- Matrix::rowSums(target_counts[, cells, drop = FALSE])
  library_size <- sum(Matrix::colSums(counts[, cells, drop = FALSE]))
  logcpm <- log2(((raw_count + 0.5) / (library_size + 1)) * 1e6)
  expression_rows[[i]] <- data.frame(
    sample_id = sample_id,
    group = cell_count_grid$group[[i]],
    cluster = cluster_id,
    n_cells = n_cells,
    gene = names(raw_count),
    raw_count = as.numeric(raw_count),
    library_size = library_size,
    logCPM_pc0.5 = as.numeric(logcpm),
    stringsAsFactors = FALSE
  )
}
expression <- do.call(rbind, expression_rows)
write.csv(expression, file.path(DIRS$tables, "04_donor_cluster_target_logCPM.csv"), row.names = FALSE)

score_rows <- list()
for (readout in names(selected_readouts)) {
  defined_genes <- selected_readouts[[readout]]
  available_genes <- intersect(defined_genes, target_genes)
  coverage <- length(available_genes) / length(defined_genes)
  min_genes <- if (length(defined_genes) == 1L) 1L else 3L
  sufficient_coverage <- length(available_genes) >= min_genes && coverage >= 0.5
  sub <- expression[expression$gene %in% available_genes, , drop = FALSE]
  unit_key <- interaction(sub$sample_id, sub$group, sub$cluster, sub$n_cells, drop = TRUE, lex.order = TRUE)
  unit_rows <- lapply(split(sub, unit_key), function(df) {
    data.frame(
      sample_id = df$sample_id[[1]],
      group = df$group[[1]],
      cluster = df$cluster[[1]],
      n_cells = df$n_cells[[1]],
      readout = readout,
      n_genes_defined = length(defined_genes),
      n_genes_available = length(available_genes),
      coverage = coverage,
      sufficient_coverage = sufficient_coverage,
      score_median_logCPM = if (sufficient_coverage) median(df$logCPM_pc0.5, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  score_rows[[readout]] <- do.call(rbind, unit_rows)
}
scores <- do.call(rbind, score_rows)
write.csv(scores, file.path(DIRS$tables, "05_donor_cluster_readout_scores.csv"), row.names = FALSE)

cluster_effect_rows <- list()
for (threshold in cell_thresholds) {
  for (readout in unique(scores$readout)) {
    for (cluster_id in all_clusters) {
      sub <- scores[
        scores$readout == readout & scores$cluster == cluster_id & scores$n_cells >= threshold,
        , drop = FALSE
      ]
      effect <- summarize_descriptive_effect(sub$score_median_logCPM, sub$group)
      cluster_effect_rows[[length(cluster_effect_rows) + 1L]] <- cbind(
        data.frame(cell_threshold = threshold, readout = readout, cluster = cluster_id, stringsAsFactors = FALSE), effect
      )
    }
  }
}
cluster_effects <- do.call(rbind, cluster_effect_rows)
write.csv(cluster_effects, file.path(DIRS$tables, "06_within_cluster_readout_effects.csv"), row.names = FALSE)

standardized_rows <- list()
for (boundary_name in names(boundaries)) {
  boundary_clusters <- boundaries[[boundary_name]]
  for (threshold in cell_thresholds) {
    for (readout in unique(scores$readout)) {
      sub <- scores[
        scores$readout == readout & scores$cluster %in% boundary_clusters & scores$n_cells >= threshold,
        c("sample_id", "group", "cluster", "score_median_logCPM"), drop = FALSE
      ]
      sample_cluster_n <- aggregate(cluster ~ sample_id + group, sub, function(x) length(unique(x)))
      complete_samples <- sample_cluster_n$sample_id[sample_cluster_n$cluster == length(boundary_clusters)]
      complete <- sub[sub$sample_id %in% complete_samples, , drop = FALSE]
      if (nrow(complete) > 0L) {
        donor_scores <- aggregate(score_median_logCPM ~ sample_id + group, complete, mean)
        colnames(donor_scores)[3] <- "equal_cluster_score"
      } else {
        donor_scores <- data.frame(sample_id = character(), group = character(), equal_cluster_score = numeric())
      }
      effect <- if (threshold == primary_cell_threshold) {
        summarize_group_effect(donor_scores$equal_cluster_score, donor_scores$group)
      } else {
        summarize_descriptive_effect(donor_scores$equal_cluster_score, donor_scores$group)
      }
      standardized_rows[[length(standardized_rows) + 1L]] <- cbind(
        data.frame(
          boundary = boundary_name,
          cell_threshold = threshold,
          inference_tier = ifelse(threshold == primary_cell_threshold, "primary_inference", "descriptive_sensitivity"),
          readout = readout,
          n_fixed_clusters = length(boundary_clusters),
          n_complete_donors = nrow(donor_scores),
          stringsAsFactors = FALSE
        ), effect
      )
      if (nrow(donor_scores) > 0L) {
        donor_scores$boundary <- boundary_name
        donor_scores$cell_threshold <- threshold
        donor_scores$readout <- readout
        donor_scores$n_fixed_clusters <- length(boundary_clusters)
        write.csv(
          donor_scores,
          file.path(DIRS$tables, paste0("07a_equal_cluster_scores_", boundary_name, "_", readout, "_ge", threshold, ".csv")),
          row.names = FALSE
        )
      }
    }
  }
}
standardized <- do.call(rbind, standardized_rows)
standardized$permutation_fdr_within_primary_threshold <- NA_real_
primary_rows <- standardized$cell_threshold == primary_cell_threshold
standardized$permutation_fdr_within_primary_threshold[primary_rows] <- p.adjust(
  standardized$permutation_p[primary_rows], method = "BH"
)
write.csv(standardized, file.path(DIRS$tables, "07_equal_cluster_composition_standardized_effects.csv"), row.names = FALSE)

summary_lines <- c(
  "# GSE254315 Composition-Sensitivity Output",
  "",
  "This directory contains a frozen, donor-aware composition sensitivity analysis.",
  "The main standardized table uses equal cluster weights within a donor and therefore",
  "tests whether a group effect remains after the predefined cluster composition is held fixed.",
  "Do not interpret persistence as causal evidence or as validation of post-hoc subclusters.",
  "",
  paste0("Input: `", INPUT_RDS, "`"),
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(summary_lines, file.path(DIRS$docs, "README.md"))

message("GSE254315 composition sensitivity analysis completed: ", Sys.time())
