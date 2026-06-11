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
ext_rds  <- file.path(locked_root, "rds/primary_locked/GSE182786_check/GSE182786_check__22_human_leydig_ext_0_14_17_with_module_scores.rds")

check_file <- function(x) {
  if (!file.exists(path.expand(x))) {
    stop("文件不存在：", x)
  }
}

check_file(core_rds)
check_file(ext_rds)


core_obj <- readRDS(path.expand(core_rds))
ext_obj  <- readRDS(path.expand(ext_rds))


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
ext_meta_cols  <- colnames(ext_obj@meta.data)

core_tan_col         <- find_score_col(core_meta_cols, "tan")
core_greenyellow_col <- find_score_col(core_meta_cols, "greenyellow")
core_green_col       <- find_score_col(core_meta_cols, "green")

ext_tan_col         <- find_score_col(ext_meta_cols, "tan")
ext_greenyellow_col <- find_score_col(ext_meta_cols, "greenyellow")
ext_green_col       <- find_score_col(ext_meta_cols, "green")


build_donor_summary <- function(seu, boundary_label,
                                tan_col, greenyellow_col, green_col) {

  required_meta <- c("sample_id", "group")
  miss_meta <- setdiff(required_meta, colnames(seu@meta.data))
  if (length(miss_meta) > 0) {
    stop("对象缺少必要元数据列：", paste(miss_meta, collapse = ", "))
  }

  vars_use <- c("sample_id", "group", tan_col, greenyellow_col, green_col)
  df <- FetchData(seu, vars = vars_use)

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
    mutate(boundary = boundary_label) %>%
    relocate(boundary, sample_id, group, n_cells)

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
      )
    ) %>%
    relocate(boundary, sample_id, group, module, score_mean, n_cells)

  list(
    donor_wide = donor_wide,
    donor_long = donor_long
  )
}

core_res <- build_donor_summary(
  seu = core_obj,
  boundary_label = "core_0_17",
  tan_col = core_tan_col,
  greenyellow_col = core_greenyellow_col,
  green_col = core_green_col
)

ext_res <- build_donor_summary(
  seu = ext_obj,
  boundary_label = "ext_0_14_17",
  tan_col = ext_tan_col,
  greenyellow_col = ext_greenyellow_col,
  green_col = ext_green_col
)

master_long <- bind_rows(core_res$donor_long, ext_res$donor_long) %>%
  mutate(
    group = factor(group, levels = c("Young", "Aged")),
    module = factor(module, levels = c("tan", "greenyellow", "green")),
    boundary = factor(boundary, levels = c("core_0_17", "ext_0_14_17"))
  ) %>%
  arrange(boundary, module, group, sample_id)

master_wide <- bind_rows(core_res$donor_wide, ext_res$donor_wide) %>%
  mutate(
    group = factor(group, levels = c("Young", "Aged")),
    boundary = factor(boundary, levels = c("core_0_17", "ext_0_14_17"))
  ) %>%
  arrange(boundary, group, sample_id)

effect_summary <- master_long %>%
  group_by(boundary, module, group) %>%
  summarise(
    mean_score = mean(score_mean, na.rm = TRUE),
    sd_score = sd(score_mean, na.rm = TRUE),
    n_donor = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mean_score, sd_score, n_donor)
  ) %>%
  mutate(
    aged_minus_young = mean_score_Aged - mean_score_Young
  ) %>%
  arrange(boundary, module)

write_csv(master_long,  file.path(base_dir, "tables", "1_module_score_by_donor_long.csv"))
write_csv(master_wide,  file.path(base_dir, "tables", "2_module_score_by_donor_wide.csv"))
write_csv(effect_summary, file.path(base_dir, "tables", "3_boundary_module_effect_summary.csv"))

write_csv(core_res$donor_wide, file.path(base_dir, "tables", "core_score_by_sample_rebuilt.csv"))
write_csv(ext_res$donor_wide,  file.path(base_dir, "tables", "ext_score_by_sample_rebuilt.csv"))

p <- ggplot(master_long, aes(x = group, y = score_mean)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
  facet_grid(module ~ boundary, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE182786 donor-level module summary",
    x = NULL,
    y = "Mean module score per donor"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "1_donor_level_module_summary.png"),
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)





