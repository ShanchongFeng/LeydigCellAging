#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv(
  "BULK_GSE287203_ROOT",
  unset = file.path(Sys.getenv("HOME"), "BULK_GSE287203_ROOT")
))
table_dir <- file.path(base_dir, "tables")
fig_dir <- file.path(base_dir, "figures")
rds_dir <- file.path(base_dir, "rds")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


deg_files <- list(
  invitro_H_vs_C = file.path(table_dir, "GSE287203_invitro_H_vs_C_DESeq2_results.csv"),
  invivo_Old_vs_Young = file.path(table_dir, "GSE287203_invivo_Old_vs_Young_DESeq2_results.csv")
)

for (f in deg_files) {
  if (!file.exists(f)) stop("Missing input file: ", f)
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

safe_dotplot <- function(enrich_obj, out_pdf, title_text, show_n = 15) {
  df <- as.data.frame(enrich_obj)

  if (is.null(enrich_obj) || nrow(df) == 0) {
    return(NULL)
  }

  pdf(out_pdf, width = 9, height = 6)
  print(
    dotplot(enrich_obj, showCategory = show_n) +
      ggtitle(title_text)
  )
  dev.off()

  invisible(TRUE)
}

prepare_gsea_list <- function(deg_df) {
  deg_clean <- deg_df %>%
    filter(!is.na(log2FoldChange)) %>%
    dplyr::select(Gene, log2FoldChange)

  gene_map <- convert_symbol_to_entrez(deg_clean$Gene)

  merged <- merge(deg_clean, gene_map, by.x = "Gene", by.y = "SYMBOL")

  merged <- merged %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene_list <- merged$log2FoldChange
  names(gene_list) <- merged$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  gene_list
}

run_one_contrast <- function(prefix, deg_file) {

  deg <- read.csv(deg_file, check.names = FALSE)

  if (!("Gene" %in% colnames(deg))) {
    stop("Gene column missing in: ", deg_file)
  }


  deg_sig <- deg %>%
    filter(!is.na(padj), padj < 0.05)

  deg_up <- deg_sig %>%
    filter(log2FoldChange >= 1)

  deg_down <- deg_sig %>%
    filter(log2FoldChange <= -1)


  write.csv(deg_sig, file.path(table_dir, paste0(prefix, "_sig_padj005.csv")), row.names = FALSE)
  write.csv(deg_up, file.path(table_dir, paste0(prefix, "_sig_up_padj005_lfc1.csv")), row.names = FALSE)
  write.csv(deg_down, file.path(table_dir, paste0(prefix, "_sig_down_padj005_lfc1.csv")), row.names = FALSE)
  gene_sets <- list(
    sig = deg_sig$Gene,
    up = deg_up$Gene,
    down = deg_down$Gene
  )

  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    map_df <- convert_symbol_to_entrez(genes)

    write.csv(
      map_df,
      file.path(table_dir, paste0(prefix, "_", set_name, "_SYMBOL_to_ENTREZ.csv")),
      row.names = FALSE
    )

    if (nrow(map_df) >= 10) {
      ego <- enrichGO(
        gene = unique(map_df$ENTREZID),
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
        file.path(table_dir, paste0(prefix, "_", set_name, "_GO_BP.csv")),
        row.names = FALSE
      )

      saveRDS(
        ego,
        file.path(rds_dir, paste0(prefix, "_", set_name, "_GO_BP.rds"))
      )

      safe_dotplot(
        ego,
        file.path(fig_dir, paste0(prefix, "_", set_name, "_GO_BP_dotplot.pdf")),
        paste0(prefix, " ", set_name, " GO BP")
      )
    } else {
    }
  }
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    map_df <- convert_symbol_to_entrez(genes)

    if (nrow(map_df) >= 10) {
      ekegg <- enrichKEGG(
        gene = unique(map_df$ENTREZID),
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
        file.path(table_dir, paste0(prefix, "_", set_name, "_KEGG.csv")),
        row.names = FALSE
      )

      saveRDS(
        ekegg,
        file.path(rds_dir, paste0(prefix, "_", set_name, "_KEGG.rds"))
      )

      safe_dotplot(
        ekegg,
        file.path(fig_dir, paste0(prefix, "_", set_name, "_KEGG_dotplot.pdf")),
        paste0(prefix, " ", set_name, " KEGG")
      )
    } else {
    }
  }
  gene_list <- prepare_gsea_list(deg)


  if (length(gene_list) >= 100) {
    gsea_kegg <- gseKEGG(
      geneList = gene_list,
      organism = "mmu",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )

    gsea_df <- as.data.frame(gsea_kegg)

    write.csv(
      gsea_df,
      file.path(table_dir, paste0(prefix, "_GSEA_KEGG_all.csv")),
      row.names = FALSE
    )

    gsea_sig <- gsea_df %>%
      filter(!is.na(p.adjust), p.adjust < 0.25)

    write.csv(
      gsea_sig,
      file.path(table_dir, paste0(prefix, "_GSEA_KEGG_padj025.csv")),
      row.names = FALSE
    )

    saveRDS(
      gsea_kegg,
      file.path(rds_dir, paste0(prefix, "_GSEA_KEGG.rds"))
    )

    safe_dotplot(
      gsea_kegg,
      file.path(fig_dir, paste0(prefix, "_GSEA_KEGG_dotplot.pdf")),
      paste0(prefix, " GSEA KEGG"),
      show_n = 20
    )
  } else {
  }
  keyword_pattern <- "steroid|cholesterol|lipid|fatty|oxidative|mitochond|cAMP|FoxO|AMPK|hormone|metabol|peroxisome|PPAR"

  keyword_tables <- list.files(
    table_dir,
    pattern = paste0("^", prefix, ".*(GO_BP|KEGG|GSEA_KEGG).*\\.csv$"),
    full.names = TRUE
  )

  keyword_hits <- list()

  for (f in keyword_tables) {
    df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next

    desc_col <- intersect(c("Description", "ID"), colnames(df))[1]
    if (is.na(desc_col)) next

    hit <- df[grep(keyword_pattern, df[[desc_col]], ignore.case = TRUE), , drop = FALSE]
    if (nrow(hit) > 0) {
      hit$source_file <- basename(f)
      keyword_hits[[basename(f)]] <- hit
    }
  }

  if (length(keyword_hits) > 0) {
    keyword_df <- bind_rows(keyword_hits)
    write.csv(
      keyword_df,
      file.path(table_dir, paste0(prefix, "_keyword_enrichment_hits.csv")),
      row.names = FALSE
    )

  } else {
  }
  summary_df <- data.frame(
    prefix = prefix,
    total_genes = nrow(deg),
    sig_padj005 = nrow(deg_sig),
    sig_up_padj005_lfc1 = nrow(deg_up),
    sig_down_padj005_lfc1 = nrow(deg_down),
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_df,
    file.path(table_dir, paste0(prefix, "_enrichment_gene_count_summary.csv")),
    row.names = FALSE
  )

  invisible(summary_df)
}

summary_all <- bind_rows(
  run_one_contrast("GSE287203_invitro_H_vs_C", deg_files$invitro_H_vs_C),
  run_one_contrast("GSE287203_invivo_Old_vs_Young", deg_files$invivo_Old_vs_Young)
)

write.csv(
  summary_all,
  file.path(table_dir, "1_bulk_enrichment_gene_count_summary_all.csv"),
  row.names = FALSE
)


