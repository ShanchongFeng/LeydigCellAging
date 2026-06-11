suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

locked_root <- path.expand(Sys.getenv("LOCKED_INPUT_ROOT", unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")))
out_root <- path.expand(Sys.getenv("CLUSTER_BOUNDARY_AUDIT_ROOT", unset = file.path(Sys.getenv("HOME"), "CLUSTER_BOUNDARY_AUDIT_ROOT")))
dir.create(file.path(out_root, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "docs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.5)
    )
}

save_plot <- function(plot, stem, width_mm = 160, height_mm = 95, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi)
    print(plot)
    grDevices::dev.off()
    ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi)
    print(plot)
    grDevices::dev.off()
  } else {
    grDevices::png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = dpi, type = "cairo")
    print(plot)
    grDevices::dev.off()
    grDevices::tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw", type = "cairo")
    print(plot)
    grDevices::dev.off()
  }
}

z_safe <- function(x) {
  if (all(!is.finite(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

wilcox_safe <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  if (length(unique(g)) != 2 || min(table(g)) < 2) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x ~ g)$p.value)
}

find_sample_col <- function(meta) {
  for (nm in c("sample_id", "donor_id", "orig.ident")) {
    if (nm %in% colnames(meta)) return(nm)
  }
  stop("No sample-like column found.")
}

get_counts <- function(obj) {
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = DefaultAssay(obj)), error = function(e) obj)
  counts <- tryCatch(
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
    error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  )
  list(obj = obj, counts = counts)
}

mouse_gene_sets <- list(
  steroidogenesis = c("Star", "Cyp11a1", "Hsd3b1", "Cyp17a1", "Hsd17b3", "Lhcgr", "Nr5a1", "Scarb1"),
  ketogenesis_fao = c("Hmgcs2", "Hmgcl", "Bdh1", "Oxct1", "Cpt1a", "Cpt2", "Acadl", "Acadm", "Hadha", "Hadh", "Acox1"),
  cholesterol_handling = c("Hmgcr", "Hmgcs1", "Sqle", "Lss", "Fdft1", "Dhcr7", "Scarb1", "Soat1"),
  detox_antioxidant = c("Gsta1", "Gsta2", "Gstm1", "Gstm2", "Nqo1", "Gpx1", "Gpx3", "Prdx1", "Sod2")
)

human_gene_sets <- list(
  steroidogenesis = c("STAR", "CYP11A1", "HSD3B2", "HSD3B1", "CYP17A1", "HSD17B3", "LHCGR", "NR5A1", "SCARB1"),
  ketogenesis_fao = c("HMGCS2", "HMGCL", "BDH1", "OXCT1", "CPT1A", "CPT2", "ACADL", "ACADM", "HADHA", "HADH", "ACOX1"),
  cholesterol_handling = c("HMGCR", "HMGCS1", "SQLE", "LSS", "FDFT1", "DHCR7", "SCARB1", "SOAT1"),
  detox_antioxidant = c("GSTA1", "GSTA2", "GSTM1", "GSTM2", "NQO1", "GPX1", "GPX3", "PRDX1", "SOD2")
)

mouse_celltype_sets <- list(
  Leydig = c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1"),
  Sertoli = c("Sox9", "Amh", "Clu", "Wt1", "Fshr"),
  Myoid_Fibro = c("Acta2", "Tagln", "Myh11", "Col1a1", "Col1a2", "Dcn"),
  Endothelial = c("Pecam1", "Kdr", "Cdh5", "Vwf", "Eng"),
  Immune = c("Ptprc", "Lst1", "C1qa", "Lyz2", "Adgre1"),
  Germ = c("Ddx4", "Dazl", "Sycp3", "Prm1", "Tnp1")
)

human_celltype_sets <- list(
  Leydig = c("INSL3", "CYP11A1", "STAR", "HSD3B2", "HSD3B1", "CYP17A1", "LHCGR", "NR5A1"),
  Sertoli = c("SOX9", "AMH", "CLU", "WT1", "FSHR"),
  Myoid_Fibro = c("ACTA2", "TAGLN", "MYH11", "COL1A1", "COL1A2", "DCN"),
  Endothelial = c("PECAM1", "KDR", "CDH5", "VWF", "ENG"),
  Immune = c("PTPRC", "LST1", "C1QA", "LYZ", "CD68"),
  Germ = c("DDX4", "DAZL", "SYCP3", "PRM1", "TNP1")
)

concept_map <- data.frame(
  concept = c("INSL3", "CYP11A1", "STAR", "HSD3B", "CYP17A1", "LHCGR", "NR5A1", "HMGCS2"),
  mouse = c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1", "Hmgcs2"),
  human = c("INSL3", "CYP11A1", "STAR", "HSD3B2", "CYP17A1", "LHCGR", "NR5A1", "HMGCS2"),
  stringsAsFactors = FALSE
)

dataset_specs <- list(
  list(
    dataset = "GSE270931_mouse_main",
    file = file.path(locked_root, "rds/primary_locked/GSE270931_mouse_scRNA/leydig_scRNA_phase2a__results__qc__04_filtered_clustered_seurat.rds"),
    species = "mouse",
    core = c("0", "4"),
    ext = c("0", "4", "14")
  ),
  list(
    dataset = "GSE270931_mouse_harmony",
    file = file.path(locked_root, "rds/primary_locked/GSE270931_harmony/leydig_scRNA_phase2a__phase2b__rds__01_harmony_clustered.rds"),
    species = "mouse",
    core = c("0", "15"),
    ext = c("0", "15")
  ),
  list(
    dataset = "GSE303193_mouse_external",
    file = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__03_filtered_clustered_seurat.rds"),
    species = "mouse",
    core = c("16"),
    ext = c("16", "17", "19")
  ),
  list(
    dataset = "GSE254315_human_targeted",
    file = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__03_filtered_clustered_human_seurat.rds"),
    species = "human",
    core = c("1", "6", "18"),
    ext = c("1", "6", "18", "7")
  ),
  list(
    dataset = "GSE182786_human_validation",
    file = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__results_01_basic_clustered.rds"),
    species = "human",
    core = c("0", "17"),
    ext = c("0", "14", "17")
  )
)

cluster_scores_path <- file.path(out_root, "tables", "cluster_audit_cluster_scores.csv")
if (!file.exists(cluster_scores_path)) {
  stop("Missing cluster audit scores. Run 3_cluster_boundary_audit.R first.")
}
cluster_scores <- read.csv(cluster_scores_path, stringsAsFactors = FALSE)
cluster_scores$cluster <- as.character(cluster_scores$cluster)

make_boundaries <- function(spec) {
  ds_scores <- cluster_scores %>%
    filter(dataset == spec$dataset) %>%
    arrange(marker_rank)
  top_unselected <- ds_scores %>%
    filter(selection == "not_selected", marker_rank <= 6) %>%
    pull(cluster)
  top_unselected <- head(top_unselected, 4)

  variants <- list(
    CMB = spec$core,
    ESB = spec$ext
  )
  if (length(spec$core) > 1) {
    for (cl in spec$core) {
      remaining <- setdiff(spec$core, cl)
      if (length(remaining) > 0) {
        variants[[paste0("CMB_without_", cl)]] <- remaining
      }
      variants[[paste0("core_only_", cl)]] <- cl
    }
  } else {
    variants[[paste0("core_only_", spec$core[1])]] <- spec$core[1]
  }
  if (length(top_unselected) >= 1) {
    variants[[paste0("CMB_plus_", top_unselected[1])]] <- unique(c(spec$core, top_unselected[1]))
  }
  if (length(top_unselected) >= 2) {
    variants[[paste0("CMB_plus_", paste(top_unselected[1:2], collapse = "_"))]] <- unique(c(spec$core, top_unselected[1:2]))
  }
  if (length(top_unselected) >= 3) {
    variants[[paste0("CMB_plus_top", length(top_unselected))]] <- unique(c(spec$core, top_unselected))
  }
  variants
}

aggregate_boundary <- function(spec, obj, counts, boundary_name, clusters) {
  meta <- obj@meta.data
  meta$cell_barcode_tmp <- rownames(meta)
  meta$cluster_tmp <- as.character(meta$seurat_clusters)
  sample_col <- find_sample_col(meta)
  if (!"group" %in% colnames(meta)) stop("Missing group column in ", spec$dataset)

  cells <- rownames(meta)[meta$cluster_tmp %in% clusters & !is.na(meta$group)]
  if (length(cells) == 0) return(NULL)
  sub_meta <- meta[cells, , drop = FALSE]
  sub_counts <- counts[, cells, drop = FALSE]
  samples <- unique(as.character(sub_meta[[sample_col]]))

  genes_target <- if (spec$species == "mouse") {
    unique(c("Hmgcs2", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1", "Foxo3", "Insl3"))
  } else {
    unique(c("HMGCS2", "CYP11A1", "STAR", "HSD3B2", "HSD3B1", "CYP17A1", "LHCGR", "NR5A1", "FOXO3", "INSL3"))
  }
  gene_sets <- if (spec$species == "mouse") mouse_gene_sets else human_gene_sets
  genes_all <- unique(c(genes_target, unlist(gene_sets, use.names = FALSE)))
  genes_present <- genes_all[genes_all %in% rownames(sub_counts)]
  if (length(genes_present) == 0) return(NULL)

  pb_counts <- sapply(samples, function(s) {
    idx <- as.character(sub_meta[[sample_col]]) == s
    Matrix::rowSums(sub_counts[genes_present, idx, drop = FALSE])
  })
  if (is.null(dim(pb_counts))) {
    pb_counts <- matrix(pb_counts, ncol = 1)
    rownames(pb_counts) <- genes_present
    colnames(pb_counts) <- samples
  }
  lib_sizes <- sapply(samples, function(s) {
    idx <- as.character(sub_meta[[sample_col]]) == s
    sum(Matrix::colSums(sub_counts[, idx, drop = FALSE]))
  })
  log_cpm <- sweep(pb_counts + 0.5, 2, lib_sizes + 1, "/")
  log_cpm <- log2(log_cpm * 1e6)

  sample_info <- unique(sub_meta[, c(sample_col, "group"), drop = FALSE])
  colnames(sample_info) <- c("sample_id", "group")
  sample_info$sample_id <- as.character(sample_info$sample_id)
  sample_info <- sample_info[match(colnames(log_cpm), sample_info$sample_id), , drop = FALSE]
  sample_info$group <- factor(as.character(sample_info$group), levels = c("Young", "Aged"))

  gene_rows <- list()
  for (gene in intersect(genes_target, rownames(log_cpm))) {
    vals <- as.numeric(log_cpm[gene, ])
    young <- vals[sample_info$group == "Young"]
    aged <- vals[sample_info$group == "Aged"]
    gene_rows[[gene]] <- data.frame(
      dataset = spec$dataset,
      species = spec$species,
      boundary = boundary_name,
      clusters = paste(clusters, collapse = "+"),
      feature_type = "target_gene",
      feature = gene,
      n_samples = length(vals),
      n_young = sum(sample_info$group == "Young", na.rm = TRUE),
      n_aged = sum(sample_info$group == "Aged", na.rm = TRUE),
      mean_young = mean(young, na.rm = TRUE),
      mean_aged = mean(aged, na.rm = TRUE),
      aged_minus_young = mean(aged, na.rm = TRUE) - mean(young, na.rm = TRUE),
      wilcox_p = wilcox_safe(vals, sample_info$group),
      n_cells = length(cells),
      stringsAsFactors = FALSE
    )
  }

  program_rows <- list()
  for (program in names(gene_sets)) {
    genes_use <- intersect(gene_sets[[program]], rownames(log_cpm))
    if (length(genes_use) < 2) next
    vals <- colMeans(log_cpm[genes_use, , drop = FALSE], na.rm = TRUE)
    young <- vals[sample_info$group == "Young"]
    aged <- vals[sample_info$group == "Aged"]
    program_rows[[program]] <- data.frame(
      dataset = spec$dataset,
      species = spec$species,
      boundary = boundary_name,
      clusters = paste(clusters, collapse = "+"),
      feature_type = "program_score",
      feature = program,
      n_genes_present = length(genes_use),
      genes_present = paste(genes_use, collapse = ";"),
      n_samples = length(vals),
      n_young = sum(sample_info$group == "Young", na.rm = TRUE),
      n_aged = sum(sample_info$group == "Aged", na.rm = TRUE),
      mean_young = mean(young, na.rm = TRUE),
      mean_aged = mean(aged, na.rm = TRUE),
      aged_minus_young = mean(aged, na.rm = TRUE) - mean(young, na.rm = TRUE),
      wilcox_p = wilcox_safe(vals, sample_info$group),
      n_cells = length(cells),
      stringsAsFactors = FALSE
    )
  }

  sample_rows <- data.frame(
    dataset = spec$dataset,
    species = spec$species,
    boundary = boundary_name,
    clusters = paste(clusters, collapse = "+"),
    sample_id = sample_info$sample_id,
    group = as.character(sample_info$group),
    n_cells = as.integer(table(factor(as.character(sub_meta[[sample_col]]), levels = sample_info$sample_id))),
    stringsAsFactors = FALSE
  )

  list(
    gene_effects = bind_rows(gene_rows),
    program_effects = bind_rows(program_rows),
    sample_cells = sample_rows
  )
}

cluster_concentration <- function(spec, meta, clusters_keep) {
  sample_col <- find_sample_col(meta)
  meta$cluster_tmp <- as.character(meta$seurat_clusters)
  rows <- list()
  for (cl in clusters_keep) {
    sub <- meta[meta$cluster_tmp == cl, , drop = FALSE]
    if (nrow(sub) == 0) next
    sample_counts <- sort(table(as.character(sub[[sample_col]])), decreasing = TRUE)
    group_counts <- table(as.character(sub$group))
    p <- as.numeric(sample_counts) / sum(sample_counts)
    entropy <- -sum(p * log(p))
    rows[[cl]] <- data.frame(
      dataset = spec$dataset,
      species = spec$species,
      cluster = cl,
      n_cells = nrow(sub),
      n_samples = length(sample_counts),
      top_sample = names(sample_counts)[1],
      top_sample_fraction = as.numeric(sample_counts[1]) / sum(sample_counts),
      effective_sample_count = exp(entropy),
      young_fraction = ifelse("Young" %in% names(group_counts), as.numeric(group_counts["Young"]) / sum(group_counts), 0),
      aged_fraction = ifelse("Aged" %in% names(group_counts), as.numeric(group_counts["Aged"]) / sum(group_counts), 0),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

cluster_marker_sets <- function(spec, obj, clusters_keep) {
  set_list <- if (spec$species == "mouse") mouse_celltype_sets else human_celltype_sets
  genes_all <- unique(unlist(set_list, use.names = FALSE))
  genes_present <- genes_all[genes_all %in% rownames(obj)]
  if (length(genes_present) == 0) return(data.frame())
  vars <- unique(c(genes_present, "seurat_clusters"))
  dat <- FetchData(obj, vars = vars)
  dat$cluster <- as.character(dat$seurat_clusters)
  dat <- dat[dat$cluster %in% clusters_keep, , drop = FALSE]
  rows <- list()
  for (cl in sort(unique(dat$cluster), method = "radix")) {
    idx <- dat$cluster == cl
    for (set_name in names(set_list)) {
      genes_use <- intersect(set_list[[set_name]], colnames(dat))
      if (length(genes_use) == 0) next
      rows[[paste(cl, set_name, sep = "__")]] <- data.frame(
        dataset = spec$dataset,
        species = spec$species,
        cluster = cl,
        marker_set = set_name,
        n_genes_present = length(genes_use),
        mean_expr = mean(as.matrix(dat[idx, genes_use, drop = FALSE]), na.rm = TRUE),
        pct_detected = mean(as.matrix(dat[idx, genes_use, drop = FALSE]) > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- bind_rows(rows)
  out <- out %>%
    group_by(dataset, marker_set) %>%
    mutate(marker_set_z = z_safe(mean_expr)) %>%
    ungroup()
  out
}

cluster_centroids <- function(spec, obj, clusters_keep) {
  genes <- if (spec$species == "mouse") concept_map$mouse else concept_map$human
  concepts <- concept_map$concept
  present <- genes %in% rownames(obj)
  genes <- genes[present]
  concepts <- concepts[present]
  if (length(genes) == 0) return(data.frame())
  dat <- FetchData(obj, vars = unique(c(genes, "seurat_clusters")))
  dat$cluster <- as.character(dat$seurat_clusters)
  dat <- dat[dat$cluster %in% clusters_keep, , drop = FALSE]
  rows <- list()
  for (cl in sort(unique(dat$cluster), method = "radix")) {
    for (i in seq_along(genes)) {
      vals <- dat[[genes[i]]][dat$cluster == cl]
      rows[[paste(cl, concepts[i], sep = "__")]] <- data.frame(
        dataset = spec$dataset,
        species = spec$species,
        cluster = cl,
        concept = concepts[i],
        gene = genes[i],
        mean_expr = mean(vals, na.rm = TRUE),
        pct_detected = mean(vals > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

all_gene_effects <- list()
all_program_effects <- list()
all_sample_cells <- list()
all_concentration <- list()
all_marker_sets <- list()
all_centroids <- list()

for (spec in dataset_specs) {
  if (!file.exists(spec$file)) {
    warning("Missing file: ", spec$file)
    next
  }
  obj <- readRDS(spec$file)
  counts_pack <- get_counts(obj)
  obj <- counts_pack$obj
  counts <- counts_pack$counts
  meta <- obj@meta.data
  meta$seurat_clusters <- as.character(meta$seurat_clusters)

  variants <- make_boundaries(spec)
  existing <- sort(unique(meta$seurat_clusters), method = "radix")
  variants <- lapply(variants, function(x) intersect(as.character(x), existing))
  variants <- variants[vapply(variants, length, integer(1)) > 0]

  ds_scores <- cluster_scores %>%
    filter(dataset == spec$dataset) %>%
    arrange(marker_rank)
  clusters_keep <- unique(c(unlist(variants, use.names = FALSE), head(ds_scores$cluster, 8)))
  clusters_keep <- intersect(as.character(clusters_keep), existing)

  for (nm in names(variants)) {
    agg <- aggregate_boundary(spec, obj, counts, nm, variants[[nm]])
    if (is.null(agg)) next
    all_gene_effects[[paste(spec$dataset, nm, sep = "__")]] <- agg$gene_effects
    all_program_effects[[paste(spec$dataset, nm, sep = "__")]] <- agg$program_effects
    all_sample_cells[[paste(spec$dataset, nm, sep = "__")]] <- agg$sample_cells
  }

  all_concentration[[spec$dataset]] <- cluster_concentration(spec, meta, clusters_keep)
  all_marker_sets[[spec$dataset]] <- cluster_marker_sets(spec, obj, clusters_keep)
  all_centroids[[spec$dataset]] <- cluster_centroids(spec, obj, clusters_keep)

  rm(obj, counts, meta)
  gc()
}

gene_effects <- bind_rows(all_gene_effects)
program_effects <- bind_rows(all_program_effects)
sample_cells <- bind_rows(all_sample_cells)
concentration <- bind_rows(all_concentration)
marker_sets <- bind_rows(all_marker_sets)
centroids <- bind_rows(all_centroids)

write.csv(gene_effects, file.path(out_root, "tables", "deep_boundary_target_gene_effects.csv"), row.names = FALSE)
write.csv(program_effects, file.path(out_root, "tables", "deep_boundary_program_effects.csv"), row.names = FALSE)
write.csv(sample_cells, file.path(out_root, "tables", "deep_boundary_sample_cell_counts.csv"), row.names = FALSE)
write.csv(concentration, file.path(out_root, "tables", "deep_cluster_donor_sample_concentration.csv"), row.names = FALSE)
write.csv(marker_sets, file.path(out_root, "tables", "deep_cluster_celltype_marker_scores.csv"), row.names = FALSE)
write.csv(centroids, file.path(out_root, "tables", "deep_cross_cohort_leydig_centroid_concepts.csv"), row.names = FALSE)

direction_summary <- bind_rows(
  gene_effects %>% mutate(table_source = "target_gene") %>% select(table_source, dataset, species, boundary, clusters, feature, aged_minus_young, wilcox_p, n_young, n_aged, n_cells),
  program_effects %>% mutate(table_source = "program_score") %>% select(table_source, dataset, species, boundary, clusters, feature, aged_minus_young, wilcox_p, n_young, n_aged, n_cells)
) %>%
  group_by(table_source, dataset, feature) %>%
  mutate(
    cmb_effect = aged_minus_young[match("CMB", boundary)],
    same_direction_as_CMB = ifelse(is.na(cmb_effect) | cmb_effect == 0, NA, sign(aged_minus_young) == sign(cmb_effect)),
    attenuation_vs_CMB = ifelse(is.na(cmb_effect), NA, abs(aged_minus_young) - abs(cmb_effect))
  ) %>%
  ungroup()

write.csv(direction_summary, file.path(out_root, "tables", "deep_boundary_direction_persistence_summary.csv"), row.names = FALSE)

key_gene <- gene_effects %>%
  filter(feature %in% c("Hmgcs2", "HMGCS2", "Cyp11a1", "CYP11A1", "Star", "STAR")) %>%
  mutate(
    feature_label = recode(feature, Hmgcs2 = "HMGCS2", HMGCS2 = "HMGCS2", Cyp11a1 = "CYP11A1", CYP11A1 = "CYP11A1", Star = "STAR", STAR = "STAR"),
    boundary = factor(boundary, levels = unique(boundary))
  )

p_gene <- ggplot(key_gene, aes(x = boundary, y = aged_minus_young, group = feature_label, colour = feature_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#808080") +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  geom_point(size = 1.45, alpha = 0.9) +
  facet_wrap(~dataset, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(HMGCS2 = "#B85C50", CYP11A1 = "#4D7188", STAR = "#6F8E5A")) +
  labs(
    title = "Boundary-expansion stress test for target genes",
    x = "Boundary variant",
    y = "Aged minus Young logCPM",
    colour = NULL
  ) +
  theme_pub(6.4) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

save_plot(p_gene, file.path(out_root, "figures", "Fig_DM1_boundary_target_gene_stress_test"), 176, 156)

key_program <- program_effects %>%
  filter(feature %in% c("steroidogenesis", "ketogenesis_fao", "cholesterol_handling", "detox_antioxidant")) %>%
  mutate(boundary = factor(boundary, levels = unique(boundary)))

p_program <- ggplot(key_program, aes(x = boundary, y = aged_minus_young, group = feature, colour = feature)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#808080") +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  geom_point(size = 1.35, alpha = 0.9) +
  facet_wrap(~dataset, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = c(
    steroidogenesis = "#4D7188",
    ketogenesis_fao = "#B85C50",
    cholesterol_handling = "#6F8E5A",
    detox_antioxidant = "#8E6C9D"
  )) +
  labs(
    title = "Boundary-expansion stress test for functional programs",
    x = "Boundary variant",
    y = "Aged minus Young mean logCPM",
    colour = NULL
  ) +
  theme_pub(6.4) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

save_plot(p_program, file.path(out_root, "figures", "Fig_DM2_boundary_program_stress_test"), 176, 156)

marker_plot <- marker_sets %>%
  left_join(cluster_scores %>% select(dataset, cluster, selection, marker_rank), by = c("dataset", "cluster")) %>%
  mutate(
    cluster_label = paste0(cluster, " (", ifelse(is.na(selection), "candidate", selection), ")"),
    cluster_label = factor(cluster_label, levels = unique(cluster_label[order(dataset, marker_rank, cluster)]))
  )

p_marker <- ggplot(marker_plot, aes(x = marker_set, y = cluster_label, fill = marker_set_z)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_gradient2(low = "#4D7188", mid = "#F5F5F5", high = "#B85C50", midpoint = 0) +
  labs(
    title = "Cell-type marker audit for selected and candidate clusters",
    x = NULL,
    y = NULL,
    fill = "Set z-score"
  ) +
  theme_pub(6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p_marker, file.path(out_root, "figures", "Fig_DM3_candidate_cluster_celltype_marker_audit"), 170, 152)

centroid_plot <- centroids %>%
  left_join(cluster_scores %>% select(dataset, cluster, selection, marker_rank), by = c("dataset", "cluster")) %>%
  group_by(dataset, concept) %>%
  mutate(concept_z = z_safe(mean_expr)) %>%
  ungroup() %>%
  mutate(
    cluster_label = paste0(cluster, " (", ifelse(is.na(selection), "candidate", selection), ")"),
    cluster_label = factor(cluster_label, levels = unique(cluster_label[order(dataset, marker_rank, cluster)]))
  )

p_centroid <- ggplot(centroid_plot, aes(x = concept, y = cluster_label, fill = concept_z)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_gradient2(low = "#4D7188", mid = "#F5F5F5", high = "#B85C50", midpoint = 0) +
  labs(
    title = "Cross-cohort Leydig centroid concept map",
    x = NULL,
    y = NULL,
    fill = "Concept z-score"
  ) +
  theme_pub(6.2) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p_centroid, file.path(out_root, "figures", "Fig_DM4_cross_cohort_leydig_centroid_concepts"), 170, 152)

key <- direction_summary %>%
  filter(feature %in% c("Hmgcs2", "HMGCS2", "Cyp11a1", "CYP11A1", "steroidogenesis", "ketogenesis_fao")) %>%
  arrange(dataset, table_source, feature, boundary)
write.csv(key, file.path(out_root, "tables", "deep_key_manuscript_candidate_effects.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "logs", "4_cross_dataset_sessionInfo.txt"))
