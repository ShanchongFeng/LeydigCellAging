#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(20260607)

MANUSCRIPT_INPUT_ROOT <- path.expand(Sys.getenv(
  "MANUSCRIPT_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "MANUSCRIPT_INPUT_ROOT")
))
SCFEA_REPO_ROOT <- path.expand(Sys.getenv(
  "SCFEA_REPO_ROOT",
  unset = file.path(Sys.getenv("HOME"), "SCFEA_REPO_ROOT")
))
SCFEA_OUTPUT_ROOT <- path.expand(Sys.getenv(
  "SCFEA_OUTPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "SCFEA_OUTPUT_ROOT")
))

INPUT_DIR <- file.path(SCFEA_OUTPUT_ROOT, "input")
TABLE_DIR <- file.path(SCFEA_OUTPUT_ROOT, "tables")
dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

max_cells_per_sample <- as.integer(Sys.getenv("SCFEA_MAX_CELLS_PER_SAMPLE", unset = "400"))
if (!is.finite(max_cells_per_sample) || max_cells_per_sample < 1) {
  stop("SCFEA_MAX_CELLS_PER_SAMPLE must be a positive integer")
}

obj_path <- file.path(
  MANUSCRIPT_INPUT_ROOT,
  "rds",
  "primary_locked",
  "GSE182786_check",
  "GSE182786_check__results_01_basic_clustered.rds"
)
module_file <- file.path(SCFEA_REPO_ROOT, "data", "module_gene_m168.csv")
if (!file.exists(obj_path)) stop("Missing object: ", obj_path)
if (!file.exists(module_file)) stop("Missing scFEA module file: ", module_file)

module_tbl <- read.csv(module_file, header = FALSE, stringsAsFactors = FALSE)
module_tbl <- module_tbl[-1, , drop = FALSE]
module_genes <- unique(unlist(module_tbl[, -1, drop = FALSE], use.names = FALSE))
module_genes <- module_genes[!is.na(module_genes) & module_genes != "" & module_genes != "A"]

message("Reading ", obj_path)
obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
md$sample_id <- as.character(md$sample_id)
md$group <- as.character(md$group)
md$cluster <- as.character(md$seurat_clusters)

cmb_cells <- rownames(md)[md$cluster %in% c("0", "17")]
if (length(cmb_cells) == 0) stop("No Leydig CMB cells found.")

selected_cells <- unlist(lapply(sort(unique(md$sample_id)), function(sample_id) {
  cells <- cmb_cells[md[cmb_cells, "sample_id"] == sample_id]
  if (length(cells) > max_cells_per_sample) {
    sample(cells, max_cells_per_sample)
  } else {
    cells
  }
}), use.names = FALSE)
selected_cells <- selected_cells[!is.na(selected_cells)]

expr <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) NULL
)
if (is.null(expr)) {
  expr <- GetAssayData(obj, assay = "RNA", slot = "data")
}

counts <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  error = function(e) NULL
)
if (is.null(counts)) {
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
}

genes_present <- intersect(module_genes, rownames(expr))
if (length(genes_present) < 50) {
  stop("Too few scFEA module genes present: ", length(genes_present))
}

mat <- as.matrix(expr[genes_present, selected_cells, drop = FALSE])
meta <- md[selected_cells, c("sample_id", "group", "cluster")]
meta$cell_id <- selected_cells
meta <- meta[, c("cell_id", "sample_id", "group", "cluster")]

cell_counts <- aggregate(cell_id ~ sample_id + group + cluster, meta, length)
colnames(cell_counts)[colnames(cell_counts) == "cell_id"] <- "n_cells"

write.csv(
  mat,
  file.path(INPUT_DIR, "GSE182786_Leydig_CMB_0_17_scFEA_input_lognorm.csv"),
  quote = FALSE
)
write.csv(
  meta,
  file.path(INPUT_DIR, "GSE182786_Leydig_CMB_0_17_scFEA_cell_metadata.csv"),
  row.names = FALSE
)
write.csv(
  cell_counts,
  file.path(INPUT_DIR, "GSE182786_Leydig_CMB_0_17_scFEA_cell_counts.csv"),
  row.names = FALSE
)

all_samples <- sort(unique(md$sample_id))
sample_vec <- md[cmb_cells, "sample_id"]
agg <- matrix(
  0,
  nrow = nrow(counts),
  ncol = length(all_samples),
  dimnames = list(rownames(counts), all_samples)
)
for (sample_id in all_samples) {
  sample_cells <- cmb_cells[sample_vec == sample_id]
  if (length(sample_cells) > 0) {
    agg[, sample_id] <- Matrix::rowSums(counts[, sample_cells, drop = FALSE])
  }
}
lib_size <- colSums(agg)
denom <- matrix(lib_size + 1, nrow = nrow(agg), ncol = ncol(agg), byrow = TRUE)
logcpm <- log2(((agg + 0.5) / denom) * 1e6)

if (!"HMGCS2" %in% rownames(logcpm)) stop("HMGCS2 not found in GSE182786 object")
sample_info <- unique(md[, c("sample_id", "group")])
hmgcs2_by_sample <- data.frame(
  sample_id = all_samples,
  hmgcs2_source = "Leydig_CMB_0_17",
  hmgcs2_logcpm = as.numeric(logcpm["HMGCS2", ]),
  leydig_n_cells = as.integer(table(sample_vec)[all_samples]),
  leydig_library_size = as.numeric(lib_size),
  stringsAsFactors = FALSE
)
hmgcs2_by_sample <- merge(hmgcs2_by_sample, sample_info, by = "sample_id", all.x = TRUE)
write.csv(
  hmgcs2_by_sample,
  file.path(TABLE_DIR, "GSE182786_Leydig_CMB_0_17_hmgcs2_by_sample.csv"),
  row.names = FALSE
)

gene_qc <- data.frame(
  n_module_genes_requested = length(module_genes),
  n_module_genes_present = length(genes_present),
  n_cells_selected = length(selected_cells),
  max_cells_per_sample = max_cells_per_sample,
  stringsAsFactors = FALSE
)
write.csv(
  gene_qc,
  file.path(INPUT_DIR, "GSE182786_Leydig_CMB_0_17_scFEA_input_qc.csv"),
  row.names = FALSE
)

message("DONE: scFEA input export")
