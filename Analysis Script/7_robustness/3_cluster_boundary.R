suppressPackageStartupMessages({
  library(Seurat)
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

theme_pub <- function(base_size = 7.2) {
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

mouse_canonical <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1")
mouse_targets <- unique(c(mouse_canonical, "Hmgcs2", "Foxo3"))
human_canonical <- c("INSL3", "CYP11A1", "STAR", "HSD3B2", "HSD3B1", "CYP17A1", "LHCGR", "NR5A1")
human_targets <- unique(c(human_canonical, "HMGCS2", "FOXO3"))

dataset_specs <- list(
  list(
    dataset = "GSE270931_mouse_main",
    file = file.path(locked_root, "rds/primary_locked/GSE270931_mouse_scRNA/leydig_scRNA_phase2a__results__qc__04_filtered_clustered_seurat.rds"),
    species = "mouse",
    core = c("0", "4"),
    ext = c("0", "4", "14"),
    canonical = mouse_canonical,
    targets = mouse_targets
  ),
  list(
    dataset = "GSE270931_mouse_harmony",
    file = file.path(locked_root, "rds/primary_locked/GSE270931_harmony/leydig_scRNA_phase2a__phase2b__rds__01_harmony_clustered.rds"),
    species = "mouse",
    core = c("0", "15"),
    ext = c("0", "15"),
    canonical = mouse_canonical,
    targets = mouse_targets
  ),
  list(
    dataset = "GSE303193_mouse_external",
    file = file.path(locked_root, "rds/primary_locked/GSE303193/GSE303193_phase1__results__03_filtered_clustered_seurat.rds"),
    species = "mouse",
    core = c("16"),
    ext = c("16", "17", "19"),
    canonical = mouse_canonical,
    targets = mouse_targets
  ),
  list(
    dataset = "GSE254315_human_targeted",
    file = file.path(locked_root, "rds/primary_locked/GSE254315/GSE254315_check__results__03_filtered_clustered_human_seurat.rds"),
    species = "human",
    core = c("1", "6", "18"),
    ext = c("1", "6", "18", "7"),
    canonical = human_canonical,
    targets = human_targets
  ),
  list(
    dataset = "GSE182786_human_validation",
    file = file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__results_01_basic_clustered.rds"),
    species = "human",
    core = c("0", "17"),
    ext = c("0", "14", "17"),
    canonical = human_canonical,
    targets = human_targets
  )
)

get_fetch_data <- function(obj, genes) {
  if ("RNA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "RNA"
  }
  obj <- tryCatch(JoinLayers(obj, assay = DefaultAssay(obj)), error = function(e) obj)
  genes_present <- genes[genes %in% rownames(obj)]
  meta_vars <- intersect(c("seurat_clusters", "group", "sample_id", "orig.ident"), colnames(obj@meta.data))
  if (!"seurat_clusters" %in% meta_vars) {
    stop("Missing seurat_clusters column")
  }
  dat <- FetchData(obj, vars = unique(c(genes_present, meta_vars)))
  dat$cluster <- as.character(dat$seurat_clusters)
  if (!"group" %in% colnames(dat)) {
    dat$group <- NA_character_
  }
  dat$group <- as.character(dat$group)
  list(data = dat, genes_present = genes_present)
}

audit_one_dataset <- function(spec) {
  if (!file.exists(spec$file)) {
    warning("Missing RDS: ", spec$file)
    return(NULL)
  }

  obj <- readRDS(spec$file)
  fetched <- get_fetch_data(obj, unique(c(spec$canonical, spec$targets)))
  dat <- fetched$data
  genes_present <- fetched$genes_present
  canonical_present <- spec$canonical[spec$canonical %in% genes_present]
  target_present <- spec$targets[spec$targets %in% genes_present]
  clusters <- sort(unique(dat$cluster), method = "radix")

  marker_rows <- list()
  for (cl in clusters) {
    idx <- dat$cluster == cl
    for (gene in target_present) {
      vals <- dat[[gene]][idx]
      marker_rows[[paste(spec$dataset, cl, gene, sep = "__")]] <- data.frame(
        dataset = spec$dataset,
        species = spec$species,
        cluster = cl,
        gene = gene,
        mean_expr = mean(vals, na.rm = TRUE),
        pct_detected = mean(vals > 0, na.rm = TRUE),
        n_cells = sum(idx),
        stringsAsFactors = FALSE
      )
    }
  }
  marker_by_cluster <- bind_rows(marker_rows)

  canonical_means <- marker_by_cluster %>%
    filter(gene %in% canonical_present) %>%
    select(cluster, gene, mean_expr)

  z_rows <- list()
  for (gene in unique(canonical_means$gene)) {
    sub <- canonical_means[canonical_means$gene == gene, , drop = FALSE]
    if (nrow(sub) == 0 || stats::sd(sub$mean_expr, na.rm = TRUE) == 0) {
      sub$z_expr <- 0
    } else {
      sub$z_expr <- as.numeric(scale(sub$mean_expr))
    }
    z_rows[[gene]] <- sub
  }

  score_df <- bind_rows(z_rows) %>%
    group_by(cluster) %>%
    summarise(
      leydig_marker_score_z = mean(z_expr, na.rm = TRUE),
      n_canonical_markers_present = n_distinct(gene),
      .groups = "drop"
    )

  detect_df <- marker_by_cluster %>%
    filter(gene %in% canonical_present) %>%
    group_by(cluster) %>%
    summarise(
      canonical_detection_fraction = mean(pct_detected, na.rm = TRUE),
      n_cells = first(n_cells),
      .groups = "drop"
    )

  group_comp <- dat %>%
    count(cluster, group, name = "n_cells_group") %>%
    group_by(cluster) %>%
    mutate(
      n_cells_cluster = sum(n_cells_group),
      fraction_within_cluster = n_cells_group / n_cells_cluster
    ) %>%
    ungroup() %>%
    mutate(dataset = spec$dataset, species = spec$species) %>%
    select(dataset, species, cluster, group, n_cells_group, n_cells_cluster, fraction_within_cluster)

  selected_core <- spec$core
  selected_ext <- spec$ext
  selected_any <- unique(c(selected_core, selected_ext))
  ext_only <- setdiff(selected_ext, selected_core)

  cluster_scores <- score_df %>%
    left_join(detect_df, by = "cluster") %>%
    mutate(
      dataset = spec$dataset,
      species = spec$species,
      selection = case_when(
        cluster %in% selected_core ~ "core",
        cluster %in% ext_only ~ "extended_only",
        TRUE ~ "not_selected"
      ),
      selected_any = cluster %in% selected_any
    ) %>%
    arrange(desc(leydig_marker_score_z), desc(canonical_detection_fraction), desc(n_cells)) %>%
    mutate(marker_rank = row_number()) %>%
    select(dataset, species, cluster, selection, selected_any, marker_rank,
           leydig_marker_score_z, canonical_detection_fraction, n_cells,
           n_canonical_markers_present)

  boundary_summary <- cluster_scores %>%
    mutate(boundary = case_when(
      cluster %in% selected_core ~ "core",
      cluster %in% selected_ext ~ "extended",
      TRUE ~ "not_selected"
    )) %>%
    filter(boundary != "not_selected") %>%
    group_by(dataset, species, boundary) %>%
    summarise(
      clusters = paste(cluster[order(marker_rank)], collapse = "+"),
      n_clusters = n(),
      total_cells = sum(n_cells),
      min_marker_rank = min(marker_rank),
      max_marker_rank = max(marker_rank),
      mean_marker_score_z = mean(leydig_marker_score_z),
      mean_detection_fraction = mean(canonical_detection_fraction),
      .groups = "drop"
    )

  missed_top <- cluster_scores %>%
    filter(marker_rank <= min(6, n()), !selected_any) %>%
    mutate(
      audit_flag = "top_marker_cluster_not_in_selected_boundary"
    )

  list(
    marker_by_cluster = marker_by_cluster,
    cluster_scores = cluster_scores,
    boundary_summary = boundary_summary,
    group_composition = group_comp,
    missed_top = missed_top
  )
}

results <- lapply(dataset_specs, audit_one_dataset)
results <- results[!vapply(results, is.null, logical(1))]

marker_by_cluster <- bind_rows(lapply(results, `[[`, "marker_by_cluster"))
cluster_scores <- bind_rows(lapply(results, `[[`, "cluster_scores"))
boundary_summary <- bind_rows(lapply(results, `[[`, "boundary_summary"))
group_composition <- bind_rows(lapply(results, `[[`, "group_composition"))
missed_top <- bind_rows(lapply(results, `[[`, "missed_top"))

write.csv(marker_by_cluster, file.path(out_root, "tables", "cluster_audit_marker_by_cluster.csv"), row.names = FALSE)
write.csv(cluster_scores, file.path(out_root, "tables", "cluster_audit_cluster_scores.csv"), row.names = FALSE)
write.csv(boundary_summary, file.path(out_root, "tables", "cluster_audit_boundary_summary.csv"), row.names = FALSE)
write.csv(group_composition, file.path(out_root, "tables", "cluster_audit_group_composition.csv"), row.names = FALSE)
write.csv(missed_top, file.path(out_root, "tables", "cluster_audit_top_marker_clusters_not_selected.csv"), row.names = FALSE)

plot_scores <- cluster_scores %>%
  group_by(dataset) %>%
  mutate(cluster_order = factor(cluster, levels = cluster[order(marker_rank, decreasing = TRUE)])) %>%
  ungroup()

p_score <- ggplot(plot_scores, aes(x = reorder(cluster, leydig_marker_score_z), y = leydig_marker_score_z, fill = selection)) +
  geom_col(width = 0.76, colour = "black", linewidth = 0.12) +
  coord_flip() +
  facet_wrap(~dataset, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(core = "#2E6F8E", extended_only = "#D4A34D", not_selected = "#C8C8C8")) +
  labs(
    title = "Cluster ranking by canonical Leydig marker score",
    x = "Cluster",
    y = "Mean z-scored marker expression",
    fill = "Boundary"
  ) +
  theme_pub(6.8) +
  theme(legend.position = "top")

save_plot(p_score, file.path(out_root, "figures", "Fig_CA1_cluster_marker_score_rank"), 170, 128)

plot_marker <- marker_by_cluster %>%
  left_join(cluster_scores %>% select(dataset, cluster, selection, marker_rank), by = c("dataset", "cluster")) %>%
  filter(selection != "not_selected" | marker_rank <= 4) %>%
  mutate(
    cluster_label = paste0(cluster, " (", selection, ")"),
    cluster_label = factor(cluster_label, levels = unique(cluster_label[order(dataset, marker_rank)]))
  )

p_marker <- ggplot(plot_marker, aes(x = gene, y = cluster_label)) +
  geom_point(aes(size = pct_detected, fill = mean_expr), shape = 21, colour = "black", stroke = 0.12) +
  facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "#F2F2F2", high = "#3A536B") +
  scale_size(range = c(0.4, 3.2), limits = c(0, 1)) +
  labs(
    title = "Selected boundaries versus top marker-ranked clusters",
    x = NULL,
    y = NULL,
    fill = "Mean expr.",
    size = "Detected"
  ) +
  theme_pub(6.5) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

save_plot(p_marker, file.path(out_root, "figures", "Fig_CA2_selected_boundary_marker_dot"), 170, 150)

make_dataset_summary <- function(df, missed) {
  datasets <- unique(df$dataset)
  lines <- character()
  for (ds in datasets) {
    sub <- df[df$dataset == ds, , drop = FALSE]
    top <- sub[order(sub$marker_rank), c("cluster", "selection", "marker_rank", "leydig_marker_score_z", "canonical_detection_fraction")]
    top <- head(top, 6)
    selected <- sub[sub$selected_any, c("cluster", "selection", "marker_rank", "leydig_marker_score_z", "canonical_detection_fraction")]
    selected <- selected[order(selected$marker_rank), , drop = FALSE]
    missed_ds <- missed[missed$dataset == ds, , drop = FALSE]

    selected_text <- paste(
      sprintf("%s:%s rank %d score %.2f detect %.2f",
              selected$cluster, selected$selection, selected$marker_rank,
              selected$leydig_marker_score_z, selected$canonical_detection_fraction),
      collapse = "; "
    )
    top_text <- paste(sprintf("%s(%s,r%d)", top$cluster, top$selection, top$marker_rank), collapse = ", ")
    missed_text <- if (nrow(missed_ds) == 0) {
      "none among top marker-ranked clusters"
    } else {
      paste(sprintf("%s rank %d", missed_ds$cluster, missed_ds$marker_rank), collapse = "; ")
    }

    lines <- c(
      lines,
      paste0("### ", ds),
      "",
      paste0("- Selected clusters: ", selected_text),
      paste0("- Top marker-ranked clusters: ", top_text),
      paste0("- Top unselected marker clusters to review: ", missed_text),
      ""
    )
  }
  lines
}

writeLines(capture.output(sessionInfo()), file.path(out_root, "logs", "cluster_boundary_audit_sessionInfo.txt"))

