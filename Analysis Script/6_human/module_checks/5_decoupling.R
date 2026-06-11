#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv("GSE182786_SCHEME3_ROOT", unset = file.path(Sys.getenv("HOME"), "GSE182786_SCHEME3_ROOT")))
locked_root <- path.expand(Sys.getenv("LOCKED_INPUT_ROOT", unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")))
dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

module_assign_csv <- file.path(locked_root, "tables/primary_locked/hdWGCNA/leydig_scRNA_phase2a__phase2b__hdwgcna__csv__hdwgcna_module_assignment.csv")

gene_stability_csv <- file.path(base_dir, "tables", "21_threshold_gene_stability_summary.csv")
gene_sparsity_csv  <- file.path(base_dir, "tables", "18_core_gene_sparsity_metrics.csv")
effect_table_csv   <- file.path(base_dir, "tables", "15_final_core_vs_ext_effect_table.csv")
attenuation_csv    <- file.path(base_dir, "tables", "16_final_boundary_attenuation_table.csv")
loo_summary_csv    <- file.path(base_dir, "tables", "8_core_loo_stability_summary.csv")

need_files <- c(
  module_assign_csv,
  gene_stability_csv,
  gene_sparsity_csv,
  effect_table_csv,
  attenuation_csv,
  loo_summary_csv
)

missing_files <- need_files[!file.exists(path.expand(need_files))]
if (length(missing_files) > 0) {
  stop("以下文件不存在：\n", paste(missing_files, collapse = "\n"))
}


mods <- read_csv(path.expand(module_assign_csv), show_col_types = FALSE)

if (!("gene_name" %in% colnames(mods))) {
  stop("模块分配表中缺少 gene_name 列。")
}
if (!("module" %in% colnames(mods))) {
  stop("模块分配表中缺少 module 列。")
}

mods_clean <- mods %>%
  mutate(
    gene_upper = toupper(gene_name),
    mouse_module = case_when(
      str_to_lower(module) %in% c("gray", "grey") ~ "grey",
      TRUE ~ module
    )
  ) %>%
  distinct(gene_upper, .keep_all = TRUE) %>%
  select(gene_name, gene_upper, mouse_module)


target_gene_tbl <- tribble(
  ~Gene,      ~gene_role,
  "HMGCS2",   "core_candidate",
  "CYP11A1",  "secondary_candidate",
  "LHCGR",    "support_limited_candidate",
  "FOXO3",    "support_limited_candidate",
  "NR5A1",    "support_limited_candidate",
  "INSL3",    "reference_marker",
  "STAR",     "reference_marker",
  "HSD3B2",   "reference_marker",
  "CYP17A1",  "reference_marker"
)

gene_stability <- read_csv(path.expand(gene_stability_csv), show_col_types = FALSE)
gene_sparsity  <- read_csv(path.expand(gene_sparsity_csv), show_col_types = FALSE)

gene_human <- target_gene_tbl %>%
  left_join(gene_stability, by = "Gene") %>%
  left_join(
    gene_sparsity %>%
      select(Gene, mean_count, median_count, min_count, max_count, n_zero, prop_zero, n_lt5, prop_lt5, n_lt10, prop_lt10),
    by = "Gene",
    suffix = c("", "_sparsity")
  )

gene_human <- gene_human %>%
  mutate(
    present_in_matrix = coalesce(present_in_matrix, FALSE),
    n_thresholds = coalesce(n_thresholds, 0L),
    n_passed = coalesce(n_passed, 0L),
    any_passed = coalesce(any_passed, FALSE),
    n_tested_with_lfc = coalesce(n_tested_with_lfc, 0L)
  )

classify_single_gene <- function(n_passed, n_tested_with_lfc, direction_consistent,
                                 mean_log2FC, cross_zero_across_thresholds) {
  if (is.na(n_passed) || n_passed == 0 || n_tested_with_lfc == 0) {
    return("not_testable_due_to_sparsity")
  }
  if (!is.na(cross_zero_across_thresholds) && cross_zero_across_thresholds) {
    return("threshold_sensitive_direction_mixed")
  }
  if (!is.na(direction_consistent) && direction_consistent) {
    if (!is.na(mean_log2FC) && mean_log2FC < 0) return("stable_negative")
    if (!is.na(mean_log2FC) && mean_log2FC > 0) return("stable_positive")
    return("stable_zero")
  }
  return("uncertain")
}

gene_human <- gene_human %>%
  rowwise() %>%
  mutate(
    single_gene_status = classify_single_gene(
      n_passed = n_passed,
      n_tested_with_lfc = n_tested_with_lfc,
      direction_consistent = direction_consistent,
      mean_log2FC = mean_log2FC,
      cross_zero_across_thresholds = cross_zero_across_thresholds
    )
  ) %>%
  ungroup()

effect_tbl <- read_csv(path.expand(effect_table_csv), show_col_types = FALSE)
atten_tbl  <- read_csv(path.expand(attenuation_csv), show_col_types = FALSE)
loo_tbl    <- read_csv(path.expand(loo_summary_csv), show_col_types = FALSE)

core_effect <- effect_tbl %>%
  filter(boundary == "core_0_17") %>%
  transmute(
    mouse_module = module,
    core_effect = aged_minus_young,
    core_ci_low = ci_low,
    core_ci_high = ci_high
  )

ext_effect <- effect_tbl %>%
  filter(boundary == "ext_0_14_17") %>%
  transmute(
    mouse_module = module,
    ext_effect = aged_minus_young,
    ext_ci_low = ci_low,
    ext_ci_high = ci_high
  )

attenuation_tbl <- atten_tbl %>%
  transmute(
    mouse_module = module,
    attenuation_ext_minus_core = attenuation_ext_minus_core,
    attenuation_ci_low = ci_low,
    attenuation_ci_high = ci_high
  )

loo_mod_tbl <- loo_tbl %>%
  transmute(
    mouse_module = module,
    full_aged_minus_young,
    full_direction,
    direction_consistency_rate,
    n_cross_zero_vs_full,
    cross_zero_any,
    min_diff,
    max_diff,
    range_diff
  )

classify_module_status <- function(core_effect, direction_consistency_rate, cross_zero_any) {
  if (is.na(core_effect)) return("module_not_evaluable")
  if (core_effect < 0) {
    if (!is.na(direction_consistency_rate) && direction_consistency_rate >= 1 && !is.na(cross_zero_any) && !cross_zero_any) {
      return("supported_human_program_level")
    } else {
      return("negative_but_less_stable")
    }
  }
  if (core_effect > 0) return("direction_discordant_human_module")
  return("module_zero_or_unclear")
}

final_tbl <- target_gene_tbl %>%
  left_join(
    mods_clean %>% select(gene_upper, mouse_module),
    by = c("Gene" = "gene_upper")
  ) %>%
  mutate(
    mouse_module = case_when(
      is.na(mouse_module) & Gene == "FOXO3" ~ "grey",
      TRUE ~ mouse_module
    )
  ) %>%
  left_join(gene_human, by = c("Gene", "gene_role")) %>%
  left_join(core_effect, by = "mouse_module") %>%
  left_join(ext_effect, by = "mouse_module") %>%
  left_join(attenuation_tbl, by = "mouse_module") %>%
  left_join(loo_mod_tbl, by = "mouse_module")

final_tbl <- final_tbl %>%
  rowwise() %>%
  mutate(
    module_status = classify_module_status(
      core_effect = core_effect,
      direction_consistency_rate = direction_consistency_rate,
      cross_zero_any = cross_zero_any
    )
  ) %>%
  ungroup()

final_tbl <- final_tbl %>%
  mutate(
    module_interpretation = case_when(
      mouse_module == "tan" ~ "HMGCS2-linked program",
      mouse_module == "greenyellow" ~ "CYP11A1-linked program",
      mouse_module == "green" ~ "NR5A1-linked program",
      mouse_module == "pink" ~ "LHCGR-linked program",
      mouse_module == "grey" ~ "unassigned",
      TRUE ~ "other_or_unannotated"
    )
  )

final_tbl <- final_tbl %>%
  mutate(
    final_evidence_tier = case_when(
      single_gene_status == "not_testable_due_to_sparsity" &
        module_status == "supported_human_program_level" ~ "program_level_only_support",

      single_gene_status == "stable_negative" &
        module_status == "supported_human_program_level" ~ "single_gene_plus_program_support",

      single_gene_status == "stable_negative" &
        module_status == "direction_discordant_human_module" ~ "single_gene_only_human_support",

      single_gene_status == "stable_positive" &
        module_status == "direction_discordant_human_module" ~ "direction_discordant_single_gene_and_module",

      single_gene_status == "stable_positive" &
        is.na(core_effect) ~ "direction_discordant_single_gene_only",

      single_gene_status == "not_testable_due_to_sparsity" &
        module_status %in% c("module_not_evaluable", "module_zero_or_unclear") ~ "not_supported_in_current_human_stress_test",

      single_gene_status == "threshold_sensitive_direction_mixed" ~ "threshold_sensitive_support_limited",

      TRUE ~ "support_limited_or_mixed"
    )
  )

out_tbl <- final_tbl %>%
  select(
    Gene, gene_role, mouse_module, module_interpretation,
    single_gene_status, n_passed, n_tested_with_lfc,
    mean_log2FC, min_log2FC, max_log2FC, range_log2FC,
    prop_zero, prop_lt5, prop_lt10,
    core_effect, core_ci_low, core_ci_high,
    ext_effect, ext_ci_low, ext_ci_high,
    attenuation_ext_minus_core, attenuation_ci_low, attenuation_ci_high,
    direction_consistency_rate, cross_zero_any, range_diff,
    module_status, final_evidence_tier
  ) %>%
  arrange(
    factor(Gene, levels = c("HMGCS2", "CYP11A1", "LHCGR", "FOXO3", "NR5A1", "INSL3", "STAR", "HSD3B2", "CYP17A1"))
  )

write_csv(out_tbl, file.path(base_dir, "tables", "22_gene_module_decoupling_master_table.csv"))

summary_tbl <- out_tbl %>%
  transmute(
    Gene,
    gene_role,
    mouse_module,
    single_gene_status,
    core_module_effect = round(core_effect, 4),
    ext_module_effect = round(ext_effect, 4),
    attenuation_ext_minus_core = round(attenuation_ext_minus_core, 4),
    loo_direction_consistency = round(direction_consistency_rate, 3),
    module_status,
    final_evidence_tier
  )

write_csv(summary_tbl, file.path(base_dir, "tables", "23_gene_module_decoupling_summary_table.csv"))

plot_df <- out_tbl %>%
  transmute(
    Gene,
    single_gene_status,
    module_status,
    final_evidence_tier
  ) %>%
  pivot_longer(
    cols = c(single_gene_status, module_status, final_evidence_tier),
    names_to = "evidence_axis",
    values_to = "label"
  ) %>%
  mutate(
    Gene = factor(Gene, levels = rev(c("HMGCS2", "CYP11A1", "LHCGR", "FOXO3", "NR5A1", "INSL3", "STAR", "HSD3B2", "CYP17A1"))),
    evidence_axis = factor(evidence_axis, levels = c("single_gene_status", "module_status", "final_evidence_tier"))
  )

p <- ggplot(plot_df, aes(x = evidence_axis, y = Gene, fill = label)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  theme_bw(base_size = 11) +
  labs(
    title = "Gene-level and module-level evidence integration",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  filename = file.path(base_dir, "figures", "10_gene_module_decoupling_tile.png"),
  plot = p,
  width = 11,
  height = 6,
  dpi = 300
)




