#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
set.seed(12345)

ROOT <- path.expand(Sys.getenv(
  "GSE270931_SCRNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "GSE270931_SCRNA_ROOT")
))
DIRS <- list(
  rds = file.path(ROOT, "rds"),
  tables = file.path(ROOT, "tables"),
  metadata = file.path(ROOT, "metadata"),
  session = file.path(ROOT, "session"),
  docs = file.path(ROOT, "docs")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

session_file <- file.path(DIRS$session, "2_states_sessionInfo.txt")
on.exit({
  writeLines(capture.output(sessionInfo()), session_file)
}, add = TRUE)

input_rds <- file.path(DIRS$rds, "GSE270931_filtered_clustered_seurat.rds")
if (!file.exists(input_rds)) stop("Missing locked GSE270931 object: ", input_rds)

sce <- readRDS(input_rds)
DefaultAssay(sce) <- "RNA"

leydig_markers <- c("Insl3", "Cyp11a1", "Star", "Hsd3b1", "Cyp17a1", "Hmgcs2")
genes_marker_use <- leydig_markers[leydig_markers %in% rownames(sce)]

group_cluster_table <- as.data.frame(table(sce$group, sce$seurat_clusters))
colnames(group_cluster_table) <- c("group", "cluster", "n_cells")
group_cluster_table <- group_cluster_table %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(group) %>%
  mutate(proportion_within_group = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(proportion_within_cluster = n_cells / sum(n_cells)) %>%
  ungroup()
write.csv(group_cluster_table, file.path(DIRS$tables, "GSE270931_cluster_group_composition.csv"), row.names = FALSE)

avg_marker <- AverageExpression(
  sce,
  features = genes_marker_use,
  group.by = "seurat_clusters",
  assays = "RNA",
  layer = "data"
)$RNA

avg_marker_df <- as.data.frame(avg_marker)
avg_marker_df$gene <- rownames(avg_marker_df)
avg_marker_df <- avg_marker_df %>% relocate(gene)
write.csv(avg_marker_df, file.path(DIRS$tables, "GSE270931_cluster_marker_average_expression.csv"), row.names = FALSE)

marker_by_cluster <- as.data.frame(t(as.matrix(avg_marker)))
marker_by_cluster$cluster <- sub("^g", "", rownames(marker_by_cluster))

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
  left_join(marker_by_cluster, by = "cluster")

state_base$leydig_marker_average_expression_score <- rowMeans(
  as.data.frame(lapply(genes_marker_use, function(g) safe_col(state_base, g))),
  na.rm = TRUE
)
state_base$hmgcs2_cyp11a1_score <- rowMeans(
  data.frame(Hmgcs2 = safe_col(state_base, "Hmgcs2"), Cyp11a1 = safe_col(state_base, "Cyp11a1")),
  na.rm = TRUE
)

score_cut <- quantile(state_base$leydig_marker_average_expression_score, 0.75, na.rm = TRUE)
hc_cut <- quantile(state_base$hmgcs2_cyp11a1_score, 0.75, na.rm = TRUE)

state_base <- state_base %>%
  mutate(
    biological_state_label = case_when(
      leydig_marker_average_expression_score >= score_cut &
        hmgcs2_cyp11a1_score >= hc_cut &
        young_fraction >= aged_fraction ~ "Leydig_Hmgcs2_Cyp11a1_high_Young_enriched",
      leydig_marker_average_expression_score >= score_cut &
        aged_fraction > young_fraction ~ "Leydig_remodeled_Aged_enriched",
      leydig_marker_average_expression_score >= score_cut ~ "Leydig_marker_high_mixed",
      TRUE ~ "non_Leydig_or_other"
    )
  )

state_composition <- group_cluster_table %>%
  left_join(
    state_base %>%
      select(cluster, total_cells, young_cells, aged_cells, young_fraction, aged_fraction,
             all_of(genes_marker_use), leydig_marker_average_expression_score,
             hmgcs2_cyp11a1_score, biological_state_label),
    by = "cluster"
  ) %>%
  arrange(as.integer(cluster), group)

write.csv(state_composition, file.path(DIRS$tables, "GSE270931_leydig_state_composition.csv"), row.names = FALSE)
write.csv(state_base, file.path(DIRS$tables, "GSE270931_leydig_state_summary_by_cluster.csv"), row.names = FALSE)


writeLines(capture.output(sessionInfo()), session_file)
