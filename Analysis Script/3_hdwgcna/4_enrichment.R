#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv(
  "HDWGCNA_ROOT",
  unset = file.path(Sys.getenv("HOME"), "HDWGCNA_ROOT")
))
table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
enrich_dir <- file.path(base_dir, "enrich")
enrich_table_dir <- file.path(enrich_dir, "tables")
enrich_fig_dir <- file.path(enrich_dir, "figures")

dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_fig_dir, recursive = TRUE, showWarnings = FALSE)

mods_file <- file.path(table_dir, "hdwgcna_module_assignment.csv")

if (!file.exists(mods_file)) {
  stop("Missing module assignment file: ", mods_file)
}

mods <- read.csv(mods_file, check.names = FALSE)


modules_use <- setdiff(unique(mods$module), "grey")


# Save module gene lists
for (m in modules_use) {
  genes <- unique(mods$gene_name[mods$module == m])
  genes <- genes[!is.na(genes) & genes != ""]

  write.csv(
    data.frame(gene = genes),
    file.path(enrich_table_dir, paste0("module_", m, "_genes.csv")),
    row.names = FALSE
  )

  write.table(
    genes,
    file.path(enrich_table_dir, paste0("module_", m, "_gene_list.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

convert_symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]

  if (length(symbols) == 0) {
    return(data.frame())
  }

  gene_df <- suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
  )

  gene_df <- gene_df[!duplicated(gene_df$SYMBOL), , drop = FALSE]
  gene_df
}

plot_enrich_dot <- function(df, title, outfile, n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df2 <- df |>
    dplyr::arrange(p.adjust) |>
    dplyr::slice_head(n = n)

  df2$Description <- factor(df2$Description, levels = rev(df2$Description))

  p <- ggplot(
    df2,
    aes(
      x = GeneRatio,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL
    )

  ggsave(outfile, p, width = 8.5, height = 6)
  invisible(p)
}

run_module_enrich <- function(module_name) {

  genes <- unique(mods$gene_name[mods$module == module_name])
  genes <- genes[!is.na(genes) & genes != ""]

  gene_df <- convert_symbol_to_entrez(genes)


  write.csv(
    gene_df,
    file.path(enrich_table_dir, paste0("module_", module_name, "_SYMBOL_to_ENTREZ.csv")),
    row.names = FALSE
  )

  ego <- NULL
  ekegg <- NULL

  if (nrow(gene_df) >= 10) {
    ego <- enrichGO(
      gene = unique(gene_df$ENTREZID),
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    ego_df <- as.data.frame(ego)

    write.csv(
      ego_df,
      file.path(enrich_table_dir, paste0("module_", module_name, "_GO_BP.csv")),
      row.names = FALSE
    )

    plot_enrich_dot(
      ego_df,
      paste0(module_name, " module GO BP"),
      file.path(enrich_fig_dir, paste0("module_", module_name, "_GO_BP_dotplot.png")),
      n = 15
    )

    ekegg <- enrichKEGG(
      gene = unique(gene_df$ENTREZID),
      organism = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    }

    ekegg_df <- as.data.frame(ekegg)

    write.csv(
      ekegg_df,
      file.path(enrich_table_dir, paste0("module_", module_name, "_KEGG.csv")),
      row.names = FALSE
    )

    plot_enrich_dot(
      ekegg_df,
      paste0(module_name, " module KEGG"),
      file.path(enrich_fig_dir, paste0("module_", module_name, "_KEGG_dotplot.png")),
      n = 15
    )
  } else {
  }

  invisible(list(ego = ego, ekegg = ekegg))
}

res_list <- list()

for (m in modules_use) {
  res_list[[m]] <- run_module_enrich(m)
}

# Keyword hits
keyword_pattern <- "steroid|cholesterol|lipid|fatty|PPAR|AMPK|mitochond|oxidative|hormone|metabol|peroxisome|ketone|ketogenesis|cAMP"

enrich_files <- list.files(
  enrich_table_dir,
  pattern = "module_.*_(GO_BP|KEGG)\\.csv$",
  full.names = TRUE
)

keyword_hits <- list()

for (f in enrich_files) {
  df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  if (!("Description" %in% colnames(df))) next

  hit <- df[grep(keyword_pattern, df$Description, ignore.case = TRUE), , drop = FALSE]

  if (nrow(hit) > 0) {
    hit$source_file <- basename(f)
    keyword_hits[[basename(f)]] <- hit
  }
}

if (length(keyword_hits) > 0) {
  keyword_df <- dplyr::bind_rows(keyword_hits)

  write.csv(
    keyword_df,
    file.path(enrich_table_dir, "current_module_keyword_enrichment_hits.csv"),
    row.names = FALSE
  )

} else {
}

# Target module interpretation table
target_genes <- c(
  "Hmgcs2", "Cyp11a1", "Cyp17a1", "Star", "Hsd3b1",
  "Lhcgr", "Insl3", "Nr5a1", "Foxo3"
)

target_mods <- mods |>
  dplyr::filter(gene_name %in% target_genes) |>
  dplyr::arrange(module, gene_name)

write.csv(
  target_mods,
  file.path(enrich_table_dir, "current_target_gene_module_assignment.csv"),
  row.names = FALSE
)


