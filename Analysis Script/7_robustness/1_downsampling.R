suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

locked_root <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
out_dir <- path.expand(Sys.getenv(
  "DONOR_SAMPLE_DOWNSAMPLING_ROOT",
  unset = file.path(Sys.getenv("HOME"), "DONOR_SAMPLE_DOWNSAMPLING_ROOT")
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260529)
n_iter <- 500
min_cells_per_sample <- 100
max_cells_per_sample <- 200
pseudocount <- 0.5

datasets <- list(
  list(
    dataset = "GSE303193",
    boundary = "CMB_16",
    species = "mouse",
    rds = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__00_leydig_cluster16.rds")
  ),
  list(
    dataset = "GSE303193",
    boundary = "ESB_16_17_19",
    species = "mouse",
    rds = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__04_leydig_cluster16_17_19.rds")
  ),
  list(
    dataset = "GSE254315",
    boundary = "CORE_1_6_18",
    species = "human",
    rds = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__11_human_leydig_core_1_6_18_with_module_scores.rds")
  ),
  list(
    dataset = "GSE254315",
    boundary = "EXT_1_6_18_7",
    species = "human",
    rds = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__12_human_leydig_ext_1_6_18_7.rds")
  ),
  list(
    dataset = "GSE182786",
    boundary = "CORE_0_17",
    species = "human",
    rds = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__21_human_leydig_core_0_17_with_module_scores.rds")
  ),
  list(
    dataset = "GSE182786",
    boundary = "EXT_0_14_17",
    species = "human",
    rds = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__22_human_leydig_ext_0_14_17_with_module_scores.rds")
  )
)

target_genes <- list(
  mouse = c("Hmgcs2", "Cyp11a1", "Insl3", "Star", "Hsd3b1", "Hsd17b3", "Cyp17a1", "Lhcgr", "Foxo3", "Nr5a1", "Scarb1", "Ppara", "Pparg", "Rxra", "Ar"),
  human = c("HMGCS2", "CYP11A1", "INSL3", "STAR", "HSD3B1", "HSD3B2", "HSD17B3", "CYP17A1", "LHCGR", "FOXO3", "NR5A1", "SCARB1", "PPARA", "PPARG", "RXRA", "AR")
)

gene_sets <- list(
  mouse = list(
    hmgcs2_axis = c("Hmgcs2", "Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl", "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a"),
    hmgcs2_axis_no_hmgcs2 = c("Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl", "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a"),
    steroidogenesis = c("Star", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Hsd17b3", "Lhcgr", "Nr5a1", "Scarb1", "Tspo"),
    leydig_identity = c("Nr5a1", "Insl3", "Lhcgr", "Star", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Hsd17b3"),
    lipid_cholesterol = c("Scarb1", "Star", "Tspo", "Cyp11a1", "Cyp17a1", "Ppara", "Pparg", "Rxra"),
    stress_senescence = c("Cdkn1a", "Cdkn2a", "Trp53", "Fos", "Jun", "Atf3", "Ddit3", "Hmox1", "Gadd45a", "Gadd45b", "Btg2")
  ),
  human = list(
    hmgcs2_axis = c("HMGCS2", "HMGCL", "BDH1", "ACAT1", "ACADM", "ACADL", "ACADVL", "CPT1A", "CPT2", "HADH", "HADHA", "PPARA", "PPARGC1A"),
    hmgcs2_axis_no_hmgcs2 = c("HMGCL", "BDH1", "ACAT1", "ACADM", "ACADL", "ACADVL", "CPT1A", "CPT2", "HADH", "HADHA", "PPARA", "PPARGC1A"),
    steroidogenesis = c("STAR", "CYP11A1", "CYP17A1", "HSD3B1", "HSD3B2", "HSD17B3", "LHCGR", "NR5A1", "SCARB1", "TSPO"),
    leydig_identity = c("NR5A1", "INSL3", "LHCGR", "STAR", "CYP11A1", "CYP17A1", "HSD3B1", "HSD3B2", "HSD17B3"),
    lipid_cholesterol = c("SCARB1", "STAR", "TSPO", "CYP11A1", "CYP17A1", "PPARA", "PPARG", "RXRA"),
    stress_senescence = c("CDKN1A", "CDKN2A", "TP53", "FOS", "JUN", "ATF3", "DDIT3", "HMOX1", "GADD45A", "GADD45B", "BTG2")
  )
)

get_counts <- function(obj) {
  if (inherits(obj[["RNA"]], "Assay5") && length(grep("^counts\\.", Layers(obj[["RNA"]]), value = TRUE)) > 1) {
    obj <- JoinLayers(obj, assay = "RNA")
  }
  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  list(obj = obj, counts = as(counts, "dgCMatrix"))
}

pick_col <- function(md, candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

safe_wilcox <- function(x, g) {
  if (length(unique(g)) != 2) return(NA_real_)
  if (sum(g == "Aged") < 2 || sum(g == "Young") < 2) return(NA_real_)
  tryCatch(wilcox.test(x[g == "Aged"], x[g == "Young"], exact = FALSE)$p.value, error = function(e) NA_real_)
}

compute_sample_readouts <- function(counts, md, cells, species, score_cols) {
  sub_counts <- counts[, cells, drop = FALSE]
  lib_size <- sum(Matrix::colSums(sub_counts))
  gene_rows <- intersect(target_genes[[species]], rownames(counts))
  gene_counts <- Matrix::rowSums(sub_counts[gene_rows, , drop = FALSE])
  gene_logcpm <- log2(((as.numeric(gene_counts) + pseudocount) / (lib_size + 1)) * 1e6)
  names(gene_logcpm) <- paste0("gene_logCPM:", gene_rows)

  set_vals <- c()
  for (set_name in names(gene_sets[[species]])) {
    genes <- intersect(gene_sets[[species]][[set_name]], rownames(counts))
    if (length(genes) >= 2) {
      counts_vec <- Matrix::rowSums(sub_counts[genes, , drop = FALSE])
      vals <- log2(((as.numeric(counts_vec) + pseudocount) / (lib_size + 1)) * 1e6)
      set_vals[paste0("geneset_logCPM_mean:", set_name)] <- mean(vals, na.rm = TRUE)
      set_vals[paste0("geneset_detected_n:", set_name)] <- sum(counts_vec > 0)
    }
  }

  score_vals <- c()
  if (length(score_cols) > 0) {
    cell_md <- md[cells, , drop = FALSE]
    for (sc in score_cols) {
      score_vals[paste0("metadata_score:", sc)] <- mean(cell_md[[sc]], na.rm = TRUE)
    }
  }
  c(gene_logcpm, set_vals, score_vals)
}

summarize_effects <- function(iter_df) {
  split_df <- split(iter_df, paste(iter_df$dataset, iter_df$boundary, iter_df$readout, sep = "|||"))
  rows <- lapply(split_df, function(df) {
    effects <- df$aged_minus_young
    pvals <- df$p_value
    data.frame(
      dataset = df$dataset[[1]],
      boundary = df$boundary[[1]],
      readout = df$readout[[1]],
      n_iterations = nrow(df),
      median_effect = median(effects, na.rm = TRUE),
      mean_effect = mean(effects, na.rm = TRUE),
      ci95_low = unname(quantile(effects, 0.025, na.rm = TRUE)),
      ci95_high = unname(quantile(effects, 0.975, na.rm = TRUE)),
      pct_aged_lower = mean(effects < 0, na.rm = TRUE),
      pct_aged_higher = mean(effects > 0, na.rm = TRUE),
      direction_consistency = max(mean(effects < 0, na.rm = TRUE), mean(effects > 0, na.rm = TRUE)),
      median_p_value = median(pvals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$direction <- ifelse(out$pct_aged_lower >= 0.9, "aged_lower_consistent",
                          ifelse(out$pct_aged_higher >= 0.9, "aged_higher_consistent",
                                 ifelse(out$pct_aged_lower >= 0.7, "aged_lower_moderate",
                                        ifelse(out$pct_aged_higher >= 0.7, "aged_higher_moderate", "mixed"))))
  out[order(out$dataset, out$boundary, out$readout), ]
}

all_effects <- list()
all_design <- list()
all_sample_counts <- list()

for (cfg in datasets) {
  obj <- readRDS(cfg$rds)
  got <- get_counts(obj)
  obj <- got$obj
  counts <- got$counts
  md <- obj@meta.data
  md$cell_barcode <- rownames(md)
  rownames(md) <- colnames(obj)

  sample_col <- pick_col(md, c("sample_id", "donor_id", "donor", "sample", "orig.ident"))
  group_col <- pick_col(md, c("group", "age_group", "condition"))
  if (is.na(sample_col) || is.na(group_col)) {
    warning("Skipping ", cfg$dataset, " ", cfg$boundary, ": no sample/group column")
    next
  }
  score_cols <- grep("^(tan|greenyellow|green)_module|module1$|Module|Score|score", colnames(md), value = TRUE)
  score_cols <- score_cols[vapply(md[score_cols], is.numeric, logical(1))]

  sample_group <- unique(md[, c(sample_col, group_col), drop = FALSE])
  colnames(sample_group) <- c("sample", "group")
  sample_group$sample <- as.character(sample_group$sample)
  sample_group$group <- as.character(sample_group$group)
  sample_group <- sample_group[sample_group$group %in% c("Young", "Aged"), , drop = FALSE]
  sample_counts <- as.data.frame(table(sample = as.character(md[[sample_col]])), stringsAsFactors = FALSE)
  sample_counts <- merge(sample_counts, sample_group, by = "sample", all.x = TRUE)
  sample_counts$dataset <- cfg$dataset
  sample_counts$boundary <- cfg$boundary
  sample_counts$species <- cfg$species
  sample_counts$eligible <- sample_counts$Freq >= min_cells_per_sample & sample_counts$group %in% c("Young", "Aged")
  sample_counts <- sample_counts[order(sample_counts$group, sample_counts$sample), ]
  all_sample_counts[[length(all_sample_counts) + 1]] <- sample_counts

  eligible <- sample_counts[sample_counts$eligible, , drop = FALSE]
  if (length(unique(eligible$group)) < 2 || any(table(eligible$group) < 2)) {
    warning("Skipping ", cfg$dataset, " ", cfg$boundary, ": insufficient eligible samples")
    next
  }
  sample_n <- min(max_cells_per_sample, min(eligible$Freq))
  design <- data.frame(
    dataset = cfg$dataset,
    boundary = cfg$boundary,
    species = cfg$species,
    rds = cfg$rds,
    sample_col = sample_col,
    group_col = group_col,
    n_total_cells = ncol(obj),
    n_total_genes = nrow(obj),
    n_eligible_samples = nrow(eligible),
    n_young_samples = sum(eligible$group == "Young"),
    n_aged_samples = sum(eligible$group == "Aged"),
    min_cells_per_sample = min_cells_per_sample,
    max_cells_per_sample = max_cells_per_sample,
    sampled_cells_per_sample = sample_n,
    n_iterations = n_iter,
    score_cols = paste(score_cols, collapse = ";"),
    stringsAsFactors = FALSE
  )
  all_design[[length(all_design) + 1]] <- design

  sample_to_cells <- split(rownames(md), as.character(md[[sample_col]]))
  iter_rows <- vector("list", n_iter)
  for (iter in seq_len(n_iter)) {
    sample_metric_list <- list()
    for (i in seq_len(nrow(eligible))) {
      smp <- eligible$sample[[i]]
      group <- eligible$group[[i]]
      cells <- sample(sample_to_cells[[smp]], size = sample_n, replace = FALSE)
      vals <- compute_sample_readouts(counts, md, cells, cfg$species, score_cols)
      sample_metric_list[[i]] <- data.frame(sample = smp, group = group, readout = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)
    }
    sample_metric <- do.call(rbind, sample_metric_list)
    readouts <- split(sample_metric, sample_metric$readout)
    iter_rows[[iter]] <- do.call(rbind, lapply(readouts, function(df) {
      aged <- df$value[df$group == "Aged"]
      young <- df$value[df$group == "Young"]
      data.frame(
        dataset = cfg$dataset,
        boundary = cfg$boundary,
        species = cfg$species,
        iteration = iter,
        sampled_cells_per_sample = sample_n,
        n_aged_samples = length(aged),
        n_young_samples = length(young),
        readout = df$readout[[1]],
        aged_mean = mean(aged, na.rm = TRUE),
        young_mean = mean(young, na.rm = TRUE),
        aged_minus_young = mean(aged, na.rm = TRUE) - mean(young, na.rm = TRUE),
        p_value = safe_wilcox(df$value, df$group),
        stringsAsFactors = FALSE
      )
    }))
    if (iter %% 50 == 0) message(cfg$dataset, " ", cfg$boundary, " iteration ", iter, "/", n_iter)
  }
  all_effects[[length(all_effects) + 1]] <- do.call(rbind, iter_rows)
}

sample_counts_out <- do.call(rbind, all_sample_counts)
design_out <- do.call(rbind, all_design)
effects_out <- do.call(rbind, all_effects)
summary_out <- summarize_effects(effects_out)

write.csv(sample_counts_out, file.path(out_dir, "downsampling_sample_counts_and_eligibility.csv"), row.names = FALSE)
write.csv(design_out, file.path(out_dir, "downsampling_design.csv"), row.names = FALSE)
write.csv(effects_out, file.path(out_dir, "downsampling_iteration_effects.csv"), row.names = FALSE)
write.csv(summary_out, file.path(out_dir, "downsampling_readout_summary.csv"), row.names = FALSE)

priority_patterns <- c(
  "Hmgcs2", "HMGCS2", "Cyp11a1", "CYP11A1", "Insl3", "INSL3", "Star", "STAR",
  "hmgcs2_axis", "steroidogenesis", "leydig_identity", "tan_module1", "greenyellow_module1", "green_module1"
)
priority <- summary_out[Reduce(`|`, lapply(priority_patterns, function(x) grepl(x, summary_out$readout, fixed = TRUE))), ]
priority <- priority[order(priority$dataset, priority$boundary, -priority$direction_consistency, priority$readout), ]
write.csv(priority, file.path(out_dir, "downsampling_priority_readouts_summary.csv"), row.names = FALSE)


