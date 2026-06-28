work_dir <- path.expand(Sys.getenv(
  "FUNCTIONAL_CATEGORY_ROOT",
  unset = file.path(Sys.getenv("HOME"), "FUNCTIONAL_CATEGORY_ROOT")
))
ortholog_root <- path.expand(Sys.getenv(
  "ORTHOLOG_CONSERVATION_ROOT",
  unset = file.path(Sys.getenv("HOME"), "ORTHOLOG_CONSERVATION_ROOT")
))
tables_dir <- file.path(work_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

effect_file <- path.expand(Sys.getenv(
  "ORTHOLOG_GENE_EFFECTS_FILE",
  unset = file.path(ortholog_root, "ortholog_gene_effects.csv")
))
stopifnot(file.exists(effect_file))
effects <- read_csv(effect_file)
effects$mouse_entrez <- as.character(effects$mouse_entrez)
effects$human_entrez <- as.character(effects$human_entrez)

go_sets <- read_csv(file.path(tables_dir, "go_external_gene_sets_20260619.csv"))
reactome_sets <- read_csv(file.path(tables_dir, "reactome_external_gene_sets_20260619.csv"))
reactome_sets$entrez_id <- as.character(reactome_sets$entrez_id)

blocks <- unique(effects[, c("dataset", "boundary", "species", "aggregation")])

calc_loo_symbol <- function(source_name, set_df, category_name) {
  out <- list()
  k <- 1L
  for (i in seq_len(nrow(blocks))) {
    block <- blocks[i, ]
    genes <- unique(set_df$gene[set_df$species == block$species & set_df$category == category_name])
    block_effects <- effects[
      effects$dataset == block$dataset &
        effects$boundary == block$boundary &
        effects$species == block$species &
        effects$aggregation == block$aggregation &
        effects$assay_symbol %in% genes,
    ]
    available <- unique(block_effects$assay_symbol)
    full_effect <- median(block_effects$aged_minus_young, na.rm = TRUE)
    for (g in available) {
      vals <- block_effects$aged_minus_young[block_effects$assay_symbol != g]
      loo_effect <- median(vals, na.rm = TRUE)
      out[[k]] <- data.frame(
        source = source_name,
        category = category_name,
        dataset = block$dataset,
        boundary = block$boundary,
        species = block$species,
        aggregation = block$aggregation,
        n_available = length(available),
        removed_gene = g,
        full_effect = full_effect,
        leave_one_effect = loo_effect,
        delta_from_full = loo_effect - full_effect,
        direction_after_removal = ifelse(loo_effect < 0, "aged_lower", "aged_higher"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

calc_loo_entrez <- function(source_name, set_df, category_name) {
  out <- list()
  k <- 1L
  for (i in seq_len(nrow(blocks))) {
    block <- blocks[i, ]
    entrez_col <- ifelse(block$species == "human", "human_entrez", "mouse_entrez")
    entrez_ids <- unique(set_df$entrez_id[set_df$species == block$species & set_df$category == category_name])
    block_effects <- effects[
      effects$dataset == block$dataset &
        effects$boundary == block$boundary &
        effects$species == block$species &
        effects$aggregation == block$aggregation &
        effects[[entrez_col]] %in% entrez_ids,
    ]
    available <- unique(block_effects[[entrez_col]])
    full_effect <- median(block_effects$aged_minus_young, na.rm = TRUE)
    for (g in available) {
      vals <- block_effects$aged_minus_young[block_effects[[entrez_col]] != g]
      loo_effect <- median(vals, na.rm = TRUE)
      removed_symbol <- paste(sort(unique(block_effects$assay_symbol[block_effects[[entrez_col]] == g])), collapse = ";")
      out[[k]] <- data.frame(
        source = source_name,
        category = category_name,
        dataset = block$dataset,
        boundary = block$boundary,
        species = block$species,
        aggregation = block$aggregation,
        n_available = length(available),
        removed_gene = removed_symbol,
        removed_entrez = g,
        full_effect = full_effect,
        leave_one_effect = loo_effect,
        delta_from_full = loo_effect - full_effect,
        direction_after_removal = ifelse(loo_effect < 0, "aged_lower", "aged_higher"),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

go_loo <- calc_loo_symbol("GO", go_sets, "go_metabolic_support_union_prespecified")
reactome_loo <- calc_loo_entrez("Reactome", reactome_sets, "reactome_metabolic_support_union_prespecified")
reactome_loo$removed_entrez <- as.character(reactome_loo$removed_entrez)
go_loo$removed_entrez <- NA_character_

loo <- rbind(
  go_loo[, c("source", "category", "dataset", "boundary", "species", "aggregation", "n_available", "removed_gene", "removed_entrez", "full_effect", "leave_one_effect", "delta_from_full", "direction_after_removal")],
  reactome_loo[, c("source", "category", "dataset", "boundary", "species", "aggregation", "n_available", "removed_gene", "removed_entrez", "full_effect", "leave_one_effect", "delta_from_full", "direction_after_removal")]
)

write.csv(loo, file.path(tables_dir, "external_metabolic_union_leave_one_gene_out_20260619.csv"), row.names = FALSE)

summary_rows <- do.call(rbind, lapply(split(loo, paste(loo$source, loo$dataset, loo$boundary, loo$aggregation, sep = "|")), function(df) {
  idx <- which.max(abs(df$delta_from_full))
  data.frame(
    source = df$source[1],
    dataset = df$dataset[1],
    boundary = df$boundary[1],
    species = df$species[1],
    aggregation = df$aggregation[1],
    n_available = df$n_available[1],
    full_effect = df$full_effect[1],
    full_direction = ifelse(df$full_effect[1] < 0, "aged_lower", "aged_higher"),
    n_direction_flips = sum(df$direction_after_removal != ifelse(df$full_effect[1] < 0, "aged_lower", "aged_higher"), na.rm = TRUE),
    max_abs_delta = abs(df$delta_from_full[idx]),
    gene_max_delta = df$removed_gene[idx],
    delta_at_gene_max = df$delta_from_full[idx],
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_rows, file.path(tables_dir, "external_metabolic_union_leave_one_gene_out_summary_20260619.csv"), row.names = FALSE)

message("Wrote external metabolic union leave-one-gene-out tables")
