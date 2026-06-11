#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260602)

LOCKED_ROOT <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
PROJECT <- path.expand(Sys.getenv(
  "STATE_STRUCTURE_ROOT",
  unset = file.path(Sys.getenv("HOME"), "STATE_STRUCTURE_ROOT")
))
TABLE_DIR <- file.path(PROJECT, "tables", "2_continuum_topology")
FIG_DIR <- file.path(PROJECT, "figures", "exploratory")
DOC_DIR <- file.path(PROJECT, "docs")
LOG_DIR <- file.path(PROJECT, "logs")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

configs <- list(
  list(
    dataset = "GSE270931",
    rds = file.path(
      LOCKED_ROOT, "rds", "primary_locked", "GSE270931_mouse_scRNA",
      "leydig_scRNA_phase2a__results__qc__04_filtered_clustered_seurat.rds"
    ),
    candidate_clusters = c("0", "4", "14"),
    source_cluster = "0",
    target_cluster = "4",
    run_projection = TRUE
  ),
  list(
    dataset = "GSE303193",
    rds = file.path(
      LOCKED_ROOT, "rds", "primary_locked", "GSE303193",
      "GSE303193_phase1__results__03_filtered_clustered_seurat.rds"
    ),
    candidate_clusters = c("16", "17", "19"),
    source_cluster = NA_character_,
    target_cluster = NA_character_,
    run_projection = FALSE
  )
)

get_log_data <- function(obj) {
  assay <- obj[["RNA"]]
  data <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
    data <- tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "data"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
    )
  }
  list(obj = obj, data = data)
}

centroid_distance_table <- function(embedding, clusters, dataset, reduction) {
  cluster_levels <- unique(clusters)
  centroids <- do.call(rbind, lapply(cluster_levels, function(cluster) {
    colMeans(embedding[clusters == cluster, , drop = FALSE])
  }))
  rownames(centroids) <- cluster_levels
  d <- as.matrix(dist(centroids))
  rows <- list()
  for (i in seq_along(cluster_levels)) {
    for (j in seq_along(cluster_levels)) {
      rows[[length(rows) + 1]] <- data.frame(
        dataset = dataset,
        reduction = reduction,
        cluster_from = cluster_levels[[i]],
        cluster_to = cluster_levels[[j]],
        centroid_distance = d[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

knn_flow_table <- function(graph, clusters, candidate_clusters, dataset) {
  candidate_idx <- which(clusters %in% candidate_clusters)
  graph <- graph[candidate_idx, candidate_idx, drop = FALSE]
  candidate_labels <- clusters[candidate_idx]
  rows <- list()
  for (from_cluster in candidate_clusters) {
    from_idx <- which(candidate_labels == from_cluster)
    total_weight <- sum(graph[from_idx, , drop = FALSE])
    for (to_cluster in candidate_clusters) {
      to_idx <- which(candidate_labels == to_cluster)
      weight <- sum(graph[from_idx, to_idx, drop = FALSE])
      rows[[length(rows) + 1]] <- data.frame(
        dataset = dataset,
        cluster_from = from_cluster,
        cluster_to = to_cluster,
        n_cells_from = length(from_idx),
        n_cells_to = length(to_idx),
        directed_knn_weight = weight,
        fraction_candidate_neighbor_weight = if (total_weight > 0) weight / total_weight else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

axis_projection <- function(pca, clusters, source_cluster, target_cluster) {
  source_centroid <- colMeans(pca[clusters == source_cluster, , drop = FALSE])
  target_centroid <- colMeans(pca[clusters == target_cluster, , drop = FALSE])
  direction <- target_centroid - source_centroid
  denom <- sum(direction ^ 2)
  if (denom == 0) return(rep(NA_real_, nrow(pca)))
  as.numeric((pca - matrix(source_centroid, nrow(pca), ncol(pca), byrow = TRUE)) %*% direction / denom)
}

score_rows <- function(data, genes) {
  idx <- match(toupper(genes), toupper(rownames(data)))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(rep(NA_real_, ncol(data)))
  Matrix::colMeans(data[idx, , drop = FALSE])
}

all_distances <- list()
all_flows <- list()
all_cluster_summary <- list()
projection_cells <- NULL
projection_bins <- NULL

for (cfg in configs) {
  obj <- readRDS(cfg$rds)
  clusters <- as.character(obj$seurat_clusters)
  keep <- clusters %in% cfg$candidate_clusters
  cell_ids <- colnames(obj)
  pca <- Embeddings(obj, "pca")[, seq_len(min(20, ncol(Embeddings(obj, "pca")))), drop = FALSE]
  umap <- Embeddings(obj, "umap")
  graph <- as(obj@graphs[["RNA_nn"]], "dgCMatrix")

  all_distances[[length(all_distances) + 1]] <- centroid_distance_table(
    pca[keep, , drop = FALSE], clusters[keep], cfg$dataset, "PCA_20"
  )
  all_distances[[length(all_distances) + 1]] <- centroid_distance_table(
    umap[keep, , drop = FALSE], clusters[keep], cfg$dataset, "UMAP_2"
  )
  all_flows[[length(all_flows) + 1]] <- knn_flow_table(
    graph, clusters, cfg$candidate_clusters, cfg$dataset
  )
  counts <- as.data.frame(table(cluster = factor(clusters[keep], levels = cfg$candidate_clusters)))
  names(counts)[2] <- "n_cells"
  counts$dataset <- cfg$dataset
  all_cluster_summary[[length(all_cluster_summary) + 1]] <- counts[, c("dataset", "cluster", "n_cells")]

  if (cfg$run_projection) {
    normalized <- get_log_data(obj)
    obj <- normalized$obj
    data <- normalized$data
    axis <- axis_projection(pca, clusters, cfg$source_cluster, cfg$target_cluster)
    md <- obj@meta.data
    group <- as.character(md$group)
    compact_leaveout <- score_rows(data, c(
      "Hmgcl", "Bdh1", "Acat1", "Acadm", "Acadl", "Acadvl",
      "Cpt1a", "Cpt2", "Hadh", "Hadha", "Ppara", "Ppargc1a"
    ))
    steroidogenesis <- score_rows(data, c(
      "Insl3", "Star", "Cyp11a1", "Hsd3b1", "Cyp17a1", "Lhcgr", "Nr5a1", "Scarb1", "Tspo"
    ))
    marker_genes <- c("Hmgcs2", "Cyp11a1", "Insl3", "Star")
    marker_rows <- lapply(marker_genes, function(gene) {
      values <- score_rows(data, gene)
      data.frame(cell = cell_ids, marker = gene, expression = values, stringsAsFactors = FALSE)
    })
    marker_long <- do.call(rbind, marker_rows)
    cell_info <- data.frame(
      cell = cell_ids,
      dataset = cfg$dataset,
      cluster = clusters,
      group = group,
      axis_projection = axis,
      compact_metabolic_support_leaveout_hmgcs2 = compact_leaveout,
      steroidogenesis_comparator = steroidogenesis,
      stringsAsFactors = FALSE
    )
    cell_info <- cell_info[cell_info$cluster %in% cfg$candidate_clusters, , drop = FALSE]
    projection_cells <- cell_info
    marker_long <- merge(marker_long, cell_info[, c("cell", "cluster", "group", "axis_projection")], by = "cell")
    marker_long$axis_bin <- cut(
      marker_long$axis_projection,
      breaks = quantile(marker_long$axis_projection, probs = seq(0, 1, length.out = 11), na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )
    marker_bins <- aggregate(
      expression ~ marker + axis_bin,
      marker_long,
      mean,
      na.rm = TRUE
    )
    marker_bins$n_cells <- as.integer(table(interaction(marker_long$marker, marker_long$axis_bin))[
      interaction(marker_bins$marker, marker_bins$axis_bin)
    ])

    program_long <- rbind(
      data.frame(
        axis_projection = cell_info$axis_projection,
        program = "compact_metabolic_support_leaveout_hmgcs2",
        score = cell_info$compact_metabolic_support_leaveout_hmgcs2
      ),
      data.frame(
        axis_projection = cell_info$axis_projection,
        program = "steroidogenesis_comparator",
        score = cell_info$steroidogenesis_comparator
      )
    )
    program_long$axis_bin <- cut(
      program_long$axis_projection,
      breaks = quantile(program_long$axis_projection, probs = seq(0, 1, length.out = 11), na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )
    program_bins <- aggregate(score ~ program + axis_bin, program_long, mean, na.rm = TRUE)
    projection_bins <- list(marker_bins = marker_bins, program_bins = program_bins)
  }
  rm(obj)
  invisible(gc())
}

distance_table <- do.call(rbind, all_distances)
flow_table <- do.call(rbind, all_flows)
cluster_summary <- do.call(rbind, all_cluster_summary)

write.csv(distance_table, file.path(TABLE_DIR, "candidate_cluster_centroid_distances.csv"), row.names = FALSE)
write.csv(flow_table, file.path(TABLE_DIR, "candidate_cluster_knn_flows.csv"), row.names = FALSE)
write.csv(cluster_summary, file.path(TABLE_DIR, "candidate_cluster_cell_counts.csv"), row.names = FALSE)
write.csv(projection_cells, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_cell_values.csv"), row.names = FALSE)
write.csv(projection_bins$marker_bins, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_marker_bins.csv"), row.names = FALSE)
write.csv(projection_bins$program_bins, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_program_bins.csv"), row.names = FALSE)

projection_cluster_summary <- aggregate(
  cbind(
    axis_projection,
    compact_metabolic_support_leaveout_hmgcs2,
    steroidogenesis_comparator
  ) ~ cluster,
  projection_cells,
  median
)
projection_group_summary <- aggregate(
  cbind(
    axis_projection,
    compact_metabolic_support_leaveout_hmgcs2,
    steroidogenesis_comparator
  ) ~ group,
  projection_cells,
  median
)
projection_correlations <- data.frame(
  readout = c(
    "compact_metabolic_support_leaveout_hmgcs2",
    "steroidogenesis_comparator"
  ),
  spearman_rho_vs_axis = c(
    cor(
      projection_cells$axis_projection,
      projection_cells$compact_metabolic_support_leaveout_hmgcs2,
      method = "spearman"
    ),
    cor(
      projection_cells$axis_projection,
      projection_cells$steroidogenesis_comparator,
      method = "spearman"
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(projection_cluster_summary, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_cluster_summary.csv"), row.names = FALSE)
write.csv(projection_group_summary, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_group_summary.csv"), row.names = FALSE)
write.csv(projection_correlations, file.path(TABLE_DIR, "GSE270931_c0_to_c4_axis_correlations.csv"), row.names = FALSE)

p_marker <- ggplot(projection_bins$marker_bins, aes(x = axis_bin, y = expression, color = marker)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_wrap(~ marker, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "Hmgcs2" = "#B8573F",
    "Cyp11a1" = "#176B87",
    "Insl3" = "#4F7D5C",
    "Star" = "#B28B33"
  )) +
  labs(x = "Canonical c0 to aged-remodeled c4 projection decile", y = "Mean log-normalized expression") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

p_program <- ggplot(projection_bins$program_bins, aes(x = axis_bin, y = score, color = program)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(
    "compact_metabolic_support_leaveout_hmgcs2" = "#B8573F",
    "steroidogenesis_comparator" = "#176B87"
  )) +
  labs(x = "Canonical c0 to aged-remodeled c4 projection decile", y = "Mean log-normalized score", color = NULL) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "WP2_GSE270931_c0_to_c4_marker_trends.pdf"), p_marker, width = 10, height = 2.8)
ggsave(file.path(FIG_DIR, "WP2_GSE270931_c0_to_c4_marker_trends.png"), p_marker, width = 10, height = 2.8, dpi = 300)
ggsave(file.path(FIG_DIR, "WP2_GSE270931_c0_to_c4_program_trends.pdf"), p_program, width = 5.2, height = 3.4)
ggsave(file.path(FIG_DIR, "WP2_GSE270931_c0_to_c4_program_trends.png"), p_program, width = 5.2, height = 3.4, dpi = 300)

capture.output(sessionInfo(), file = file.path(LOG_DIR, "4_topology_sessionInfo.txt"))

