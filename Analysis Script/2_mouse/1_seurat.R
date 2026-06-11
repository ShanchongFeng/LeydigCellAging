#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

ROOT <- path.expand(Sys.getenv(
  "GSE270931_SCRNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE270931_SCRNA_ROOT")
))
OLD_BASE <- path.expand(Sys.getenv(
  "GSE270931_RAW_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE270931_RAW_ROOT")
))
DIRS <- list(
  raw_data = file.path(ROOT, "raw_data"),
  scripts = file.path(ROOT, "scripts"),
  rds = file.path(ROOT, "rds"),
  tables = file.path(ROOT, "tables"),
  figures = file.path(ROOT, "figures"),
  logs = file.path(ROOT, "logs"),
  metadata = file.path(ROOT, "metadata"),
  session = file.path(ROOT, "session"),
  checksums = file.path(ROOT, "checksums"),
  docs = file.path(ROOT, "docs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

session_file <- file.path(DIRS$session, "1_seurat_sessionInfo.txt")
on.exit({
  writeLines(capture.output(sessionInfo()), session_file)
}, add = TRUE)


raw_tar <- file.path(OLD_BASE, "GSE270931_RAW.tar")
young_file <- file.path(OLD_BASE, "GSM8365314_young.csv.gz")
aged_file <- file.path(OLD_BASE, "GSM8365315_aged.csv.gz")
source_script <- path.expand(Sys.getenv(
  "GSE270931_SOURCE_SCRIPT",
  unset = file.path(
    path.expand(Sys.getenv(
      "GSE270931_SCRIPT_ROOT",
      unset = file.path(Sys.getenv("HOME"), "GSE270931_SCRIPT_ROOT")
    )),
    "1_seurat.R"
  )
))

input_paths <- data.frame(
  role = c("raw_tar", "young_matrix", "aged_matrix", "old_source_script"),
  path = c(raw_tar, young_file, aged_file, source_script),
  exists = file.exists(c(raw_tar, young_file, aged_file, source_script)),
  stringsAsFactors = FALSE
)
write.csv(input_paths, file.path(DIRS$metadata, "GSE270931_input_paths.csv"), row.names = FALSE)
if (!all(input_paths$exists)) {
  stop("Missing GSE270931 input(s); see metadata/GSE270931_input_paths.csv")
}

params <- data.frame(
  parameter = c(
    "seed", "min.cells", "min.features", "percent.mt.pattern",
    "qc.nFeature_RNA.lower", "qc.nFeature_RNA.upper", "qc.percent.mt.upper",
    "variable.features.method", "variable.features.n",
    "pca.neighbor.dims", "cluster.resolution",
    "CMB_boundary", "ESB_boundary", "FindMarkers.ident1", "FindMarkers.ident2",
    "FindMarkers.logfc.threshold", "FindMarkers.min.pct",
    "note"
  ),
  value = c(
    "12345", "3", "200", "^mt-",
    "300", "6000", "15",
    "vst", "2000",
    "1:20", "0.4",
    "historical core clusters 0+4", "historical extended clusters 0+4+14",
    "Aged", "Young",
    "0", "0.1",
    "FindMarkers is cell-level discovery and must not be interpreted as animal-level inference"
  ),
  stringsAsFactors = FALSE
)
write.csv(params, file.path(DIRS$metadata, "GSE270931_analysis_parameters.csv"), row.names = FALSE)

young <- read.table(
  gzfile(young_file),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

aged <- read.csv(
  gzfile(aged_file),
  skip = 7,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!("Cell_Index" %in% colnames(aged))) {
  stop("Aged matrix does not contain Cell_Index column")
}

young_mat <- t(as.matrix(young))
rownames(aged) <- aged$Cell_Index
aged$Cell_Index <- NULL
aged_mat <- t(as.matrix(aged))

storage.mode(young_mat) <- "numeric"
storage.mode(aged_mat) <- "numeric"
colnames(young_mat) <- paste0("Young_", colnames(young_mat))
colnames(aged_mat) <- paste0("Aged_", colnames(aged_mat))

common_genes <- intersect(rownames(young_mat), rownames(aged_mat))
young_mat <- young_mat[common_genes, , drop = FALSE]
aged_mat <- aged_mat[common_genes, , drop = FALSE]
expr_mat <- cbind(young_mat, aged_mat)
expr_mat <- round(expr_mat)
storage.mode(expr_mat) <- "integer"

matrix_summary <- data.frame(
  item = c("young_genes", "young_cells", "aged_genes", "aged_cells", "common_genes", "merged_genes", "merged_cells"),
  value = c(nrow(young_mat), ncol(young_mat), nrow(aged_mat), ncol(aged_mat),
            length(common_genes), nrow(expr_mat), ncol(expr_mat))
)
write.csv(matrix_summary, file.path(DIRS$metadata, "GSE270931_raw_matrix_dimensions.csv"), row.names = FALSE)
saveRDS(expr_mat, file.path(DIRS$rds, "GSE270931_expr_mat_merged.rds"))

meta <- data.frame(
  cell = colnames(expr_mat),
  group = ifelse(grepl("^Young_", colnames(expr_mat)), "Young", "Aged"),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$cell

sce <- CreateSeuratObject(
  counts = expr_mat,
  meta.data = meta,
  project = "GSE270931_mouse_testis",
  min.cells = 3,
  min.features = 200
)
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
saveRDS(sce, file.path(DIRS$rds, "GSE270931_raw_seurat.rds"))

raw_qc <- data.frame(
  group = sce$group,
  nFeature_RNA = sce$nFeature_RNA,
  nCount_RNA = sce$nCount_RNA,
  percent.mt = sce$percent.mt,
  stringsAsFactors = FALSE
)
write.csv(raw_qc, file.path(DIRS$metadata, "GSE270931_raw_cell_qc_metrics.csv"), row.names = FALSE)

sce <- subset(
  sce,
  subset = nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)
saveRDS(sce, file.path(DIRS$rds, "GSE270931_filtered_preprocess_input.rds"))

filtered_qc <- data.frame(
  cell = colnames(sce),
  group = sce$group,
  nFeature_RNA = sce$nFeature_RNA,
  nCount_RNA = sce$nCount_RNA,
  percent.mt = sce$percent.mt,
  stringsAsFactors = FALSE
)
write.csv(filtered_qc, file.path(DIRS$metadata, "GSE270931_filtered_cell_qc_metrics.csv"), row.names = FALSE)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.4)
sce <- RunUMAP(sce, dims = 1:20)

saveRDS(sce, file.path(DIRS$rds, "GSE270931_filtered_clustered_seurat.rds"))

cluster_table <- as.data.frame(table(sce$seurat_clusters))
colnames(cluster_table) <- c("cluster", "n_cells")
write.csv(cluster_table, file.path(DIRS$tables, "GSE270931_cluster_cell_counts.csv"), row.names = FALSE)

group_cluster_table <- as.data.frame(table(sce$group, sce$seurat_clusters))
colnames(group_cluster_table) <- c("group", "cluster", "n_cells")
group_cluster_table <- group_cluster_table %>%
  group_by(group) %>%
  mutate(proportion_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(proportion_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(group_cluster_table, file.path(DIRS$tables, "GSE270931_cluster_group_composition.csv"), row.names = FALSE)

p_umap_cluster <- DimPlot(sce, reduction = "umap", label = TRUE) +
  ggtitle("GSE270931 mouse scRNA clusters")
p_umap_group <- DimPlot(sce, reduction = "umap", group.by = "group") +
  ggtitle("GSE270931 mouse scRNA by group")
ggsave(file.path(DIRS$figures, "GSE270931_umap_by_cluster.pdf"), p_umap_cluster, width = 7, height = 6)
ggsave(file.path(DIRS$figures, "GSE270931_umap_by_group.pdf"), p_umap_group, width = 7, height = 6)

leydig_markers <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Hmgcs2")
target_genes <- c("Hmgcs2", "Cyp11a1", "Lhcgr", "Star", "Hsd3b1", "Hsd3b", "Foxo3", "Nr5a1", "Hdac1")
marker_and_target <- unique(c(leydig_markers, target_genes))

marker_presence <- data.frame(
  gene = marker_and_target,
  present = marker_and_target %in% rownames(sce),
  stringsAsFactors = FALSE
)
write.csv(marker_presence, file.path(DIRS$tables, "GSE270931_marker_presence.csv"), row.names = FALSE)

genes_marker_use <- leydig_markers[leydig_markers %in% rownames(sce)]
genes_target_use <- target_genes[target_genes %in% rownames(sce)]

p_marker_dot <- DotPlot(sce, features = genes_marker_use) + RotatedAxis()
ggsave(file.path(DIRS$figures, "GSE270931_leydig_marker_dotplot.pdf"), p_marker_dot, width = 9, height = 4)

avg_marker <- AverageExpression(
  sce,
  features = genes_marker_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)
avg_marker_df <- as.data.frame(avg_marker$RNA)
avg_marker_df$gene <- rownames(avg_marker_df)
avg_marker_df <- avg_marker_df %>% relocate(gene)
write.csv(avg_marker_df, file.path(DIRS$tables, "GSE270931_cluster_marker_average_expression.csv"), row.names = FALSE)

cluster_marker_wide <- t(as.matrix(avg_marker$RNA))
cluster_marker_wide <- as.data.frame(cluster_marker_wide)
cluster_marker_wide$cluster <- rownames(cluster_marker_wide)

safe_col <- function(df, nm) {
  if (nm %in% colnames(df)) df[[nm]] else rep(NA_real_, nrow(df))
}

state_base <- group_cluster_table %>%
  group_by(cluster) %>%
  summarise(
    total_cells = sum(n_cells),
    young_cells = sum(n_cells[group == "Young"]),
    aged_cells = sum(n_cells[group == "Aged"]),
    young_fraction = young_cells / total_cells,
    aged_fraction = aged_cells / total_cells,
    .groups = "drop"
  ) %>%
  left_join(cluster_marker_wide, by = "cluster")

state_base$leydig_marker_score <- rowMeans(
  as.data.frame(lapply(genes_marker_use, function(g) safe_col(state_base, g))),
  na.rm = TRUE
)
names(state_base)[names(state_base) == "leydig_marker_score"] <- "leydig_marker_average_expression_score"
state_base$hmgcs2_cyp11a1_score <- rowMeans(
  data.frame(Hmgcs2 = safe_col(state_base, "Hmgcs2"), Cyp11a1 = safe_col(state_base, "Cyp11a1")),
  na.rm = TRUE
)

score_cut <- quantile(state_base$leydig_marker_average_expression_score, 0.75, na.rm = TRUE)
hc_cut <- quantile(state_base$hmgcs2_cyp11a1_score, 0.75, na.rm = TRUE)
state_base$biological_state_label <- ifelse(
  state_base$leydig_marker_average_expression_score >= score_cut &
    state_base$hmgcs2_cyp11a1_score >= hc_cut &
    state_base$young_fraction >= state_base$aged_fraction,
  "Leydig_Hmgcs2_Cyp11a1_high_Young_enriched",
  ifelse(
    state_base$leydig_marker_average_expression_score >= score_cut &
      state_base$aged_fraction > state_base$young_fraction,
    "Leydig_remodeled_Aged_enriched",
    ifelse(
      state_base$leydig_marker_average_expression_score >= score_cut,
      "Leydig_marker_high_mixed",
      "non_Leydig_or_other"
    )
  )
)

state_composition <- group_cluster_table %>%
  left_join(
    state_base %>%
      select(cluster, total_cells, young_cells, aged_cells, young_fraction, aged_fraction,
             leydig_marker_average_expression_score, hmgcs2_cyp11a1_score,
             biological_state_label),
    by = "cluster"
  ) %>%
  arrange(cluster, group)
write.csv(state_composition, file.path(DIRS$tables, "GSE270931_leydig_state_composition.csv"), row.names = FALSE)

all_markers <- FindAllMarkers(
  sce,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(all_markers, file.path(DIRS$tables, "GSE270931_all_cluster_markers.csv"), row.names = FALSE)

Idents(sce) <- "seurat_clusters"
cluster_ids <- levels(Idents(sce))
CMB_clusters <- c("0", "4")
ESB_clusters <- c("0", "4", "14")
if (!all(CMB_clusters %in% cluster_ids)) {
  stop("Historical CMB clusters 0+4 not found. Existing clusters: ", paste(cluster_ids, collapse = ", "))
}
if (!all(ESB_clusters %in% cluster_ids)) {
  warning("Historical ESB cluster 14 not found; falling back to CMB only")
  ESB_clusters <- CMB_clusters
}

leydig_CMB <- subset(sce, idents = CMB_clusters)
leydig_ESB <- subset(sce, idents = ESB_clusters)
saveRDS(leydig_CMB, file.path(DIRS$rds, "GSE270931_leydig_subset_CMB.rds"))
saveRDS(leydig_ESB, file.path(DIRS$rds, "GSE270931_leydig_subset_ESB.rds"))

subset_counts <- bind_rows(
  as.data.frame(table(leydig_CMB$group, leydig_CMB$seurat_clusters)) %>%
    setNames(c("group", "cluster", "n_cells")) %>%
    mutate(boundary = "CMB_0_4"),
  as.data.frame(table(leydig_ESB$group, leydig_ESB$seurat_clusters)) %>%
    setNames(c("group", "cluster", "n_cells")) %>%
    mutate(boundary = "ESB_0_4_14")
) %>%
  group_by(boundary, group) %>%
  mutate(proportion_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(boundary, cluster) %>%
  mutate(proportion_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(subset_counts, file.path(DIRS$tables, "GSE270931_leydig_boundary_group_composition.csv"), row.names = FALSE)

analyze_boundary <- function(obj, boundary_label) {
  genes_use <- genes_target_use[genes_target_use %in% rownames(obj)]

  avg_target <- AverageExpression(
    obj,
    features = genes_use,
    group.by = "group",
    assays = "RNA",
    layer = "data"
  )$RNA
  avg_long <- as.data.frame(avg_target)
  avg_long$gene <- rownames(avg_long)
  avg_long <- avg_long %>%
    relocate(gene) %>%
    pivot_longer(-gene, names_to = "group", values_to = "average_expression") %>%
    mutate(boundary = boundary_label)

  Idents(obj) <- "group"
  deg <- FindMarkers(
    obj,
    ident.1 = "Aged",
    ident.2 = "Young",
    logfc.threshold = 0,
    min.pct = 0.1
  )
  deg$gene <- rownames(deg)
  deg <- deg %>% relocate(gene)
  target_deg <- deg %>%
    filter(gene %in% genes_use) %>%
    mutate(
      boundary = boundary_label,
      interpretation_note = "cell-level discovery only; not animal-level inference"
    )

  list(avg_long = avg_long, deg = deg, target_deg = target_deg)
}

CMB_res <- analyze_boundary(leydig_CMB, "CMB_0_4")
ESB_res <- analyze_boundary(leydig_ESB, "ESB_0_4_14")

write.csv(CMB_res$target_deg, file.path(DIRS$tables, "GSE270931_target_gene_DEG_CMB.csv"), row.names = FALSE)
write.csv(ESB_res$target_deg, file.path(DIRS$tables, "GSE270931_target_gene_DEG_ESB.csv"), row.names = FALSE)
write.csv(CMB_res$deg, file.path(DIRS$tables, "GSE270931_all_DEG_CMB.csv"), row.names = FALSE)
write.csv(ESB_res$deg, file.path(DIRS$tables, "GSE270931_all_DEG_ESB.csv"), row.names = FALSE)
write.csv(
  bind_rows(CMB_res$avg_long, ESB_res$avg_long),
  file.path(DIRS$tables, "GSE270931_target_gene_average_expression_by_group.csv"),
  row.names = FALSE
)


writeLines(capture.output(sessionInfo()), session_file)
