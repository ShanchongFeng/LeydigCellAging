suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

locked_root <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
out_dir <- path.expand(Sys.getenv(
  "ORTHOLOG_CONSERVATION_ROOT",
  unset = file.path(Sys.getenv("HOME"), "ORTHOLOG_CONSERVATION_ROOT")
))
ref_dir <- file.path(out_dir, "reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

homologene_url <- "https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data"
homologene_file <- file.path(ref_dir, "homologene_build68.data")
if (!file.exists(homologene_file) || file.info(homologene_file)$size < 1e6) {
  ok <- tryCatch({
    download.file(homologene_url, homologene_file, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(homologene_file)) {
    stop("Could not download HomoloGene data from ", homologene_url)
  }
}

module_files <- list(
  hdwgcna_tan = file.path(locked_root, "tables/primary_locked/hdWGCNA_enrichment/leydig_scRNA_phase2a__phase2b__hdwgcna__enrich__csv__module_tan_genes.csv"),
  hdwgcna_greenyellow = file.path(locked_root, "tables/primary_locked/hdWGCNA_enrichment/leydig_scRNA_phase2a__phase2b__hdwgcna__enrich__csv__module_greenyellow_genes.csv"),
  hdwgcna_green = file.path(locked_root, "tables/primary_locked/hdWGCNA_enrichment/leydig_scRNA_phase2a__phase2b__hdwgcna__enrich__csv__module_green_genes.csv"),
  hdwgcna_pink = file.path(locked_root, "tables/primary_locked/hdWGCNA_enrichment/leydig_scRNA_phase2a__phase2b__hdwgcna__enrich__csv__module_pink_genes.csv")
)

read_module_genes <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("gene_name" %in% colnames(df)) return(unique(df$gene_name))
  if ("gene" %in% colnames(df)) return(unique(df$gene))
  unique(df[[1]])
}

mouse_modules <- lapply(module_files, read_module_genes)
mouse_modules$hdwgcna_tan_no_hmgcs2 <- setdiff(mouse_modules$hdwgcna_tan, "Hmgcs2")
mouse_modules$hmgcs2_axis <- c("Hmgcs2", "Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl", "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a")
mouse_modules$hmgcs2_axis_no_hmgcs2 <- setdiff(mouse_modules$hmgcs2_axis, "Hmgcs2")
mouse_modules$curated_steroidogenesis <- c("Insl3", "Star", "Cyp11a1", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1", "Scarb1", "Tspo")
mouse_modules$leydig_identity <- c("Nr5a1", "Insl3", "Lhcgr", "Star", "Cyp11a1", "Cyp17a1", "Hsd3b1", "Hsd17b3")

hom <- read.delim(
  homologene_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("homologene_id", "tax_id", "gene_id", "symbol", "protein_gi", "protein_accession")
)
mm <- hom[hom$tax_id == 10090, c("homologene_id", "gene_id", "symbol")]
hs <- hom[hom$tax_id == 9606, c("homologene_id", "gene_id", "symbol")]
colnames(mm) <- c("homologene_id", "mouse_entrez", "mouse_symbol")
colnames(hs) <- c("homologene_id", "human_entrez", "human_symbol")

mm_counts <- as.data.frame(table(mm$homologene_id), stringsAsFactors = FALSE)
hs_counts <- as.data.frame(table(hs$homologene_id), stringsAsFactors = FALSE)
colnames(mm_counts) <- c("homologene_id", "n_mouse_in_group")
colnames(hs_counts) <- c("homologene_id", "n_human_in_group")
orth_all <- merge(mm, hs, by = "homologene_id")
orth_all <- merge(orth_all, mm_counts, by = "homologene_id")
orth_all <- merge(orth_all, hs_counts, by = "homologene_id")
orth_all$one_to_one <- orth_all$n_mouse_in_group == 1 & orth_all$n_human_in_group == 1
orth_1to1 <- orth_all[orth_all$one_to_one, ]
orth_1to1 <- orth_1to1[!duplicated(orth_1to1$mouse_symbol), ]

write.csv(orth_all, file.path(out_dir, "homologene_mouse_human_orthologs_all.csv"), row.names = FALSE)
write.csv(orth_1to1, file.path(out_dir, "homologene_mouse_human_orthologs_one_to_one.csv"), row.names = FALSE)

module_map <- do.call(rbind, lapply(names(mouse_modules), function(module_name) {
  genes <- unique(mouse_modules[[module_name]])
  df <- data.frame(module = module_name, mouse_symbol = genes, stringsAsFactors = FALSE)
  out <- merge(df, orth_1to1[, c("homologene_id", "mouse_symbol", "human_symbol", "mouse_entrez", "human_entrez")], by = "mouse_symbol", all.x = TRUE)
  out$mapped_one_to_one <- !is.na(out$human_symbol)
  out
}))
write.csv(module_map, file.path(out_dir, "ortholog_module_gene_map_one_to_one.csv"), row.names = FALSE)

module_coverage <- do.call(rbind, lapply(split(module_map, module_map$module), function(df) {
  data.frame(
    module = df$module[[1]],
    n_mouse_genes = nrow(df),
    n_one_to_one_orthologs = sum(df$mapped_one_to_one),
    ortholog_coverage = mean(df$mapped_one_to_one),
    mouse_genes_missing_one_to_one = paste(df$mouse_symbol[!df$mapped_one_to_one], collapse = ";"),
    human_orthologs = paste(na.omit(df$human_symbol), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(module_coverage, file.path(out_dir, "ortholog_module_coverage.csv"), row.names = FALSE)

datasets <- list(
  list(
    dataset = "GSE270931",
    boundary = "CMB_0_4",
    species = "mouse",
    sample_aware = FALSE,
    rds = file.path(locked_root, "rds/primary_locked/GSE270931_mouse_scRNA/leydig_scRNA_phase2a__results__target_genes__00_leydig_subset_0_4.rds")
  ),
  list(
    dataset = "GSE270931",
    boundary = "ESB_0_4_14",
    species = "mouse",
    sample_aware = FALSE,
    rds = file.path(locked_root, "rds/primary_locked/GSE270931_mouse_scRNA/leydig_scRNA_phase2a__results__target_genes__00_leydig_subset_0_4_14.rds")
  ),
  list(
    dataset = "GSE303193",
    boundary = "CMB_16",
    species = "mouse",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__00_leydig_cluster16.rds")
  ),
  list(
    dataset = "GSE303193",
    boundary = "ESB_16_17_19",
    species = "mouse",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__04_leydig_cluster16_17_19.rds")
  ),
  list(
    dataset = "GSE254315",
    boundary = "CORE_1_6_18",
    species = "human",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__11_human_leydig_core_1_6_18_with_module_scores.rds")
  ),
  list(
    dataset = "GSE254315",
    boundary = "EXT_1_6_18_7",
    species = "human",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__12_human_leydig_ext_1_6_18_7.rds")
  ),
  list(
    dataset = "GSE182786",
    boundary = "CORE_0_17",
    species = "human",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__21_human_leydig_core_0_17_with_module_scores.rds")
  ),
  list(
    dataset = "GSE182786",
    boundary = "EXT_0_14_17",
    species = "human",
    sample_aware = TRUE,
    rds = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__22_human_leydig_ext_0_14_17_with_module_scores.rds")
  )
)

get_counts <- function(obj) {
  if ("RNA" %in% names(obj@assays) && inherits(obj[["RNA"]], "Assay5")) {
    count_layers <- grep("^counts", Layers(obj[["RNA"]]), value = TRUE)
    if (length(count_layers) > 1) obj <- JoinLayers(obj, assay = "RNA")
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

normalize_group <- function(x) {
  y <- as.character(x)
  ifelse(grepl("young|^y$|m5|5m", y, ignore.case = TRUE), "Young",
         ifelse(grepl("aged|old|^a$|m20|20m", y, ignore.case = TRUE), "Aged", y))
}

safe_wilcox <- function(x, g) {
  if (length(unique(g)) != 2) return(NA_real_)
  if (sum(g == "Aged") < 2 || sum(g == "Young") < 2) return(NA_real_)
  tryCatch(wilcox.test(x[g == "Aged"], x[g == "Young"], exact = FALSE)$p.value, error = function(e) NA_real_)
}

compute_effects <- function(cfg, orth_1to1) {
  obj <- readRDS(cfg$rds)
  got <- get_counts(obj)
  obj <- got$obj
  counts <- got$counts
  md <- obj@meta.data
  rownames(md) <- colnames(obj)
  group_col <- pick_col(md, c("group", "age_group", "condition", "orig.ident"))
  sample_col <- pick_col(md, c("sample_id", "donor_id", "donor", "sample"))
  if (is.na(group_col)) stop("No group column in ", cfg$rds)
  md$age_group_norm <- normalize_group(md[[group_col]])

  if (cfg$species == "mouse") {
    genes <- intersect(unique(orth_1to1$mouse_symbol), rownames(counts))
    map <- orth_1to1[match(genes, orth_1to1$mouse_symbol), ]
    row_gene <- map$mouse_symbol
  } else {
    genes <- intersect(unique(orth_1to1$human_symbol), rownames(counts))
    map <- orth_1to1[match(genes, orth_1to1$human_symbol), ]
    row_gene <- map$human_symbol
  }

  counts_sub <- counts[row_gene, , drop = FALSE]
  rownames(counts_sub) <- row_gene

  use_sample <- cfg$sample_aware && !is.na(sample_col)
  if (use_sample) {
    sg <- unique(data.frame(sample = as.character(md[[sample_col]]), group = md$age_group_norm, stringsAsFactors = FALSE))
    use_sample <- length(unique(sg$group)) == 2 && all(table(sg$group) >= 2) && length(unique(sg$sample)) > length(unique(sg$group))
  }

  if (use_sample) {
    samples <- unique(as.character(md[[sample_col]]))
    sample_counts <- do.call(cbind, lapply(samples, function(smp) {
      cells <- rownames(md)[as.character(md[[sample_col]]) == smp]
      Matrix::rowSums(counts_sub[, cells, drop = FALSE])
    }))
    colnames(sample_counts) <- samples
    sample_lib <- vapply(samples, function(smp) {
      cells <- rownames(md)[as.character(md[[sample_col]]) == smp]
      sum(Matrix::colSums(counts[, cells, drop = FALSE]))
    }, numeric(1))
    sample_group <- vapply(samples, function(smp) unique(md$age_group_norm[as.character(md[[sample_col]]) == smp])[[1]], character(1))
    aggregation <- "sample_mean"
  } else {
    samples <- c("Young", "Aged")
    sample_counts <- do.call(cbind, lapply(samples, function(grp) {
      cells <- rownames(md)[md$age_group_norm == grp]
      Matrix::rowSums(counts_sub[, cells, drop = FALSE])
    }))
    colnames(sample_counts) <- samples
    sample_lib <- vapply(samples, function(grp) {
      cells <- rownames(md)[md$age_group_norm == grp]
      sum(Matrix::colSums(counts[, cells, drop = FALSE]))
    }, numeric(1))
    sample_group <- samples
    aggregation <- "group_pseudobulk_no_sample"
  }

  logcpm <- sweep(as.matrix(sample_counts) + 0.5, 2, sample_lib + 1, "/")
  logcpm <- log2(logcpm * 1e6)
  aged_cols <- sample_group == "Aged"
  young_cols <- sample_group == "Young"
  effects <- rowMeans(logcpm[, aged_cols, drop = FALSE], na.rm = TRUE) - rowMeans(logcpm[, young_cols, drop = FALSE], na.rm = TRUE)
  pvals <- apply(logcpm, 1, function(x) safe_wilcox(x, sample_group))

  out <- data.frame(
    dataset = cfg$dataset,
    boundary = cfg$boundary,
    species = cfg$species,
    aggregation = aggregation,
    assay_symbol = rownames(logcpm),
    aged_minus_young = as.numeric(effects),
    p_value = as.numeric(pvals),
    n_aged_units = sum(aged_cols),
    n_young_units = sum(young_cols),
    stringsAsFactors = FALSE
  )
  if (cfg$species == "mouse") {
    out <- cbind(out, map[, c("homologene_id", "mouse_symbol", "human_symbol", "mouse_entrez", "human_entrez")])
  } else {
    out <- cbind(out, map[, c("homologene_id", "mouse_symbol", "human_symbol", "mouse_entrez", "human_entrez")])
  }
  out
}

all_effects <- do.call(rbind, lapply(datasets, compute_effects, orth_1to1 = orth_1to1))
write.csv(all_effects, file.path(out_dir, "ortholog_gene_effects.csv"), row.names = FALSE)

module_membership <- unique(module_map[module_map$mapped_one_to_one, c("module", "mouse_symbol", "human_symbol", "homologene_id")])
effects_module <- merge(all_effects, module_membership, by = c("homologene_id", "mouse_symbol", "human_symbol"), all.x = FALSE)
write.csv(effects_module, file.path(out_dir, "ortholog_gene_effects_by_module.csv"), row.names = FALSE)

summarize_module <- function(df, all_df) {
  background <- all_df[all_df$dataset == df$dataset[[1]] & all_df$boundary == df$boundary[[1]] & all_df$species == df$species[[1]], ]
  n_lower <- sum(df$aged_minus_young < 0, na.rm = TRUE)
  n_total <- sum(!is.na(df$aged_minus_young))
  data.frame(
    dataset = df$dataset[[1]],
    boundary = df$boundary[[1]],
    species = df$species[[1]],
    aggregation = df$aggregation[[1]],
    module = df$module[[1]],
    n_genes = n_total,
    median_effect = median(df$aged_minus_young, na.rm = TRUE),
    mean_effect = mean(df$aged_minus_young, na.rm = TRUE),
    frac_aged_lower = n_lower / max(1, n_total),
    binom_p_aged_lower = tryCatch(binom.test(n_lower, n_total, p = 0.5, alternative = "greater")$p.value, error = function(e) NA_real_),
    wilcox_vs_background_p = tryCatch(wilcox.test(df$aged_minus_young, background$aged_minus_young, alternative = "less", exact = FALSE)$p.value, error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}

module_summary <- do.call(rbind, lapply(split(effects_module, paste(effects_module$dataset, effects_module$boundary, effects_module$module, sep = "|||")), summarize_module, all_df = all_effects))
module_summary <- module_summary[order(module_summary$dataset, module_summary$boundary, module_summary$module), ]
write.csv(module_summary, file.path(out_dir, "ortholog_module_effect_summary.csv"), row.names = FALSE)

mouse_effects <- all_effects[all_effects$species == "mouse", ]
human_effects <- all_effects[all_effects$species == "human", ]
concordance_rows <- list()
row_i <- 1
for (mkey in unique(paste(mouse_effects$dataset, mouse_effects$boundary, sep = "|||"))) {
  mdf <- mouse_effects[paste(mouse_effects$dataset, mouse_effects$boundary, sep = "|||") == mkey, ]
  for (hkey in unique(paste(human_effects$dataset, human_effects$boundary, sep = "|||"))) {
    hdf <- human_effects[paste(human_effects$dataset, human_effects$boundary, sep = "|||") == hkey, ]
    pair <- merge(
      mdf[, c("homologene_id", "mouse_symbol", "human_symbol", "aged_minus_young")],
      hdf[, c("homologene_id", "aged_minus_young")],
      by = "homologene_id",
      suffixes = c("_mouse", "_human")
    )
    pair_modules <- merge(pair, module_membership[, c("module", "homologene_id")], by = "homologene_id")
    for (module_name in unique(pair_modules$module)) {
      pdf <- pair_modules[pair_modules$module == module_name, ]
      if (nrow(pdf) < 3) next
      concordance_rows[[row_i]] <- data.frame(
        mouse_dataset = unique(mdf$dataset),
        mouse_boundary = unique(mdf$boundary),
        human_dataset = unique(hdf$dataset),
        human_boundary = unique(hdf$boundary),
        module = module_name,
        n_genes = nrow(pdf),
        spearman_r = suppressWarnings(cor(pdf$aged_minus_young_mouse, pdf$aged_minus_young_human, method = "spearman", use = "pairwise.complete.obs")),
        pearson_r = suppressWarnings(cor(pdf$aged_minus_young_mouse, pdf$aged_minus_young_human, method = "pearson", use = "pairwise.complete.obs")),
        sign_concordance = mean(sign(pdf$aged_minus_young_mouse) == sign(pdf$aged_minus_young_human), na.rm = TRUE),
        both_aged_lower_fraction = mean(pdf$aged_minus_young_mouse < 0 & pdf$aged_minus_young_human < 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1
    }
  }
}
concordance <- do.call(rbind, concordance_rows)
write.csv(concordance, file.path(out_dir, "ortholog_cross_dataset_concordance.csv"), row.names = FALSE)

priority_modules <- c("hdwgcna_tan", "hdwgcna_tan_no_hmgcs2", "hmgcs2_axis", "hmgcs2_axis_no_hmgcs2", "curated_steroidogenesis", "hdwgcna_greenyellow", "hdwgcna_green")
plot_df <- module_summary[module_summary$module %in% priority_modules, ]
p <- ggplot(plot_df, aes(x = module, y = paste(dataset, boundary, sep = " | "), fill = median_effect)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2F6DB3", mid = "white", high = "#B54A42", midpoint = 0, na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "Median gene effect", title = "Ortholog-constrained module effects") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 10)
  )
ggsave(file.path(out_dir, "ortholog_module_effect_heatmap.pdf"), p, width = 9, height = 4.8, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "ortholog_module_effect_heatmap.png"), p, width = 9, height = 4.8, units = "in", dpi = 300)


