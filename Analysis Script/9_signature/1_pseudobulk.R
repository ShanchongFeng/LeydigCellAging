suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

locked_root <- path.expand(Sys.getenv(
  "LOCKED_INPUT_ROOT",
  unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")
))
out_dir <- path.expand(Sys.getenv(
  "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT",
  unset = file.path(Sys.getenv("HOME"), "PORTABLE_SIGNATURE_PSEUDOBULK_ROOT")
))
pb_dir <- file.path(out_dir, "pseudobulk")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

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

all_design <- list()

for (cfg in datasets) {
  obj <- readRDS(cfg$rds)
  got <- get_counts(obj)
  obj <- got$obj
  counts <- got$counts
  md <- obj@meta.data
  rownames(md) <- colnames(obj)
  sample_col <- pick_col(md, c("sample_id", "donor_id", "donor", "sample"))
  group_col <- pick_col(md, c("group", "age_group", "condition", "orig.ident"))
  if (is.na(sample_col) || is.na(group_col)) {
    warning("Skipping ", cfg$dataset, " ", cfg$boundary, ": no sample/group columns")
    next
  }
  md$sample_id_export <- as.character(md[[sample_col]])
  md$group_export <- normalize_group(md[[group_col]])
  keep <- md$group_export %in% c("Young", "Aged")
  md <- md[keep, , drop = FALSE]
  counts <- counts[, rownames(md), drop = FALSE]

  samples <- unique(md$sample_id_export)
  pb_counts <- do.call(cbind, lapply(samples, function(smp) {
    cells <- rownames(md)[md$sample_id_export == smp]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  }))
  colnames(pb_counts) <- samples
  lib_size <- Matrix::colSums(pb_counts)
  logcpm <- sweep(as.matrix(pb_counts) + 0.5, 2, lib_size + 1, "/")
  logcpm <- log2(logcpm * 1e6)

  sample_md <- do.call(rbind, lapply(samples, function(smp) {
    cells <- rownames(md)[md$sample_id_export == smp]
    data.frame(
      sample = smp,
      group = unique(md$group_export[md$sample_id_export == smp])[[1]],
      n_cells = length(cells),
      library_size = sum(Matrix::colSums(counts[, cells, drop = FALSE])),
      stringsAsFactors = FALSE
    )
  }))
  sample_md$dataset <- cfg$dataset
  sample_md$boundary <- cfg$boundary
  sample_md$species <- cfg$species
  sample_md <- sample_md[, c("dataset", "boundary", "species", "sample", "group", "n_cells", "library_size")]

  prefix <- paste(cfg$dataset, cfg$boundary, sep = "__")
  write.csv(t(logcpm), file.path(pb_dir, paste0(prefix, "__sample_logcpm.csv")), row.names = TRUE)
  write.csv(sample_md, file.path(pb_dir, paste0(prefix, "__sample_metadata.csv")), row.names = FALSE)

  all_design[[length(all_design) + 1]] <- data.frame(
    dataset = cfg$dataset,
    boundary = cfg$boundary,
    species = cfg$species,
    rds = cfg$rds,
    sample_col = sample_col,
    group_col = group_col,
    n_cells = ncol(counts),
    n_genes = nrow(counts),
    n_samples = nrow(sample_md),
    n_young = sum(sample_md$group == "Young"),
    n_aged = sum(sample_md$group == "Aged"),
    min_cells_per_sample = min(sample_md$n_cells),
    stringsAsFactors = FALSE
  )
}

design <- do.call(rbind, all_design)
write.csv(design, file.path(out_dir, "decoupler_pseudobulk_design.csv"), row.names = FALSE)

