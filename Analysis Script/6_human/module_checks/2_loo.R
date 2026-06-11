#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv("GSE182786_SCHEME3_ROOT", unset = file.path(Sys.getenv("HOME"), "GSE182786_SCHEME3_ROOT")))
locked_root <- path.expand(Sys.getenv("LOCKED_INPUT_ROOT", unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")))
dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

core_rds <- file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__21_human_leydig_core_0_17_with_module_scores.rds")

old_loo_csv <- file.path(locked_root, "tables/primary_locked/GSE182786_check/GSE182786_check__41_core_loo_module_diffs.csv")

check_file <- function(x) {
  if (!file.exists(path.expand(x))) {
    stop("文件不存在：", x)
  }
}

check_file(core_rds)


core_obj <- readRDS(path.expand(core_rds))


find_score_col <- function(meta_cols, prefix) {
  exact_target <- paste0("^", prefix, "_module1$")
  hit_exact <- grep(exact_target, meta_cols, value = TRUE, ignore.case = FALSE)
  if (length(hit_exact) == 1) return(hit_exact)

  loose_target1 <- paste0("^", prefix, "_module")
  hit1 <- grep(loose_target1, meta_cols, value = TRUE, ignore.case = FALSE)
  if (length(hit1) >= 1) return(hit1[1])

  hit2 <- grep(prefix, meta_cols, value = TRUE, ignore.case = TRUE)
  if (length(hit2) >= 1) return(hit2[1])

  stop("找不到模块分数字段：", prefix)
}

core_meta_cols <- colnames(core_obj@meta.data)

tan_col         <- find_score_col(core_meta_cols, "tan")
greenyellow_col <- find_score_col(core_meta_cols, "greenyellow")
green_col       <- find_score_col(core_meta_cols, "green")


vars_use <- c("sample_id", "group", tan_col, greenyellow_col, green_col)
df <- FetchData(core_obj, vars = vars_use)
colnames(df) <- c("sample_id", "group", "tan_score", "greenyellow_score", "green_score")

donor_wide <- df %>%
  group_by(sample_id, group) %>%
  summarise(
    tan_mean = mean(tan_score, na.rm = TRUE),
    greenyellow_mean = mean(greenyellow_score, na.rm = TRUE),
    green_mean = mean(green_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(group, sample_id)

donor_long <- donor_wide %>%
  pivot_longer(
    cols = c(tan_mean, greenyellow_mean, green_mean),
    names_to = "module",
    values_to = "score_mean"
  ) %>%
  mutate(
    module = recode(
      module,
      tan_mean = "tan",
      greenyellow_mean = "greenyellow",
      green_mean = "green"
    ),
    group = factor(group, levels = c("Young", "Aged")),
    module = factor(module, levels = c("tan", "greenyellow", "green"))
  ) %>%
  arrange(module, group, sample_id)

write_csv(donor_wide, file.path(base_dir, "tables", "4_core_donor_score_wide_for_loo.csv"))
write_csv(donor_long, file.path(base_dir, "tables", "5_core_donor_score_long_for_loo.csv"))

full_diff <- donor_long %>%
  group_by(module, group) %>%
  summarise(mean_score = mean(score_mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_score) %>%
  mutate(
    full_aged_minus_young = Aged - Young,
    full_direction = case_when(
      full_aged_minus_young > 0 ~ "positive",
      full_aged_minus_young < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  arrange(module)

samples <- unique(donor_wide$sample_id)

loo_long <- lapply(samples, function(s) {
  tmp <- donor_long %>% filter(sample_id != s)

  sum_df <- tmp %>%
    group_by(module, group) %>%
    summarise(mean_score = mean(score_mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_score) %>%
    mutate(
      left_out = s,
      diff_aged_minus_young = Aged - Young
    ) %>%
    select(left_out, module, Young, Aged, diff_aged_minus_young)

  sum_df
}) %>%
  bind_rows() %>%
  left_join(full_diff %>% select(module, full_aged_minus_young, full_direction), by = "module") %>%
  mutate(
    loo_direction = case_when(
      diff_aged_minus_young > 0 ~ "positive",
      diff_aged_minus_young < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    same_direction_as_full = loo_direction == full_direction,
    crosses_zero_vs_full = case_when(
      full_aged_minus_young == 0 ~ diff_aged_minus_young != 0,
      TRUE ~ sign(diff_aged_minus_young) != sign(full_aged_minus_young)
    )
  ) %>%
  arrange(module, left_out)

write_csv(loo_long, file.path(base_dir, "tables", "6_core_loo_module_diffs_long.csv"))

loo_wide <- loo_long %>%
  select(left_out, module, diff_aged_minus_young) %>%
  pivot_wider(names_from = module, values_from = diff_aged_minus_young) %>%
  arrange(left_out)

write_csv(loo_wide, file.path(base_dir, "tables", "7_core_loo_module_diffs_wide.csv"))

loo_summary <- loo_long %>%
  group_by(module, full_aged_minus_young, full_direction) %>%
  summarise(
    n_leave_one_out_runs = n(),
    n_same_direction_as_full = sum(same_direction_as_full, na.rm = TRUE),
    direction_consistency_rate = mean(same_direction_as_full, na.rm = TRUE),
    n_cross_zero_vs_full = sum(crosses_zero_vs_full, na.rm = TRUE),
    cross_zero_any = any(crosses_zero_vs_full, na.rm = TRUE),
    min_diff = min(diff_aged_minus_young, na.rm = TRUE),
    max_diff = max(diff_aged_minus_young, na.rm = TRUE),
    mean_diff = mean(diff_aged_minus_young, na.rm = TRUE),
    sd_diff = sd(diff_aged_minus_young, na.rm = TRUE),
    range_diff = max_diff - min_diff,
    .groups = "drop"
  ) %>%
  arrange(module)

write_csv(loo_summary, file.path(base_dir, "tables", "8_core_loo_stability_summary.csv"))

loo_influence <- loo_long %>%
  mutate(
    abs_shift_from_full = abs(diff_aged_minus_young - full_aged_minus_young)
  ) %>%
  group_by(module) %>%
  arrange(desc(abs_shift_from_full), .by_group = TRUE) %>%
  mutate(rank_within_module = row_number()) %>%
  ungroup()

write_csv(loo_influence, file.path(base_dir, "tables", "9_core_loo_influence_rank.csv"))

top_influential <- loo_influence %>%
  group_by(module) %>%
  slice_head(n = 5) %>%
  ungroup()

write_csv(top_influential, file.path(base_dir, "tables", "10_core_loo_top5_influential_donors.csv"))

p1 <- ggplot(loo_long, aes(x = left_out, y = diff_aged_minus_young)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_hline(aes(yintercept = full_aged_minus_young), color = "red", linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE182786 core leave-one-donor-out module robustness",
    x = "Left-out donor",
    y = "Aged minus Young"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "2_core_loo_module_diffs.png"),
  plot = p1,
  width = 9,
  height = 8,
  dpi = 300
)

p2 <- ggplot(loo_influence, aes(x = reorder(left_out, abs_shift_from_full), y = abs_shift_from_full)) +
  geom_col() +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Influence of leaving out each donor",
    x = "Left-out donor",
    y = "|LOO diff - full diff|"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "3_core_loo_influence_rank.png"),
  plot = p2,
  width = 9,
  height = 8,
  dpi = 300
)

if (file.exists(path.expand(old_loo_csv))) {
  old_loo <- suppressWarnings(read_csv(path.expand(old_loo_csv), show_col_types = FALSE))
  write_csv(old_loo, file.path(base_dir, "tables", "old_41_core_loo_module_diffs_backup.csv"))
}





