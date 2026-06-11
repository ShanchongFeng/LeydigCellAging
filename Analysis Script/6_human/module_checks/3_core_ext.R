#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

options(stringsAsFactors = FALSE)

base_dir <- path.expand(Sys.getenv("GSE182786_SCHEME3_ROOT", unset = file.path(Sys.getenv("HOME"), "GSE182786_SCHEME3_ROOT")))
locked_root <- path.expand(Sys.getenv("LOCKED_INPUT_ROOT", unset = file.path(Sys.getenv("HOME"), "LOCKED_INPUT_ROOT")))
dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

long_csv <- file.path(base_dir, "tables", "1_module_score_by_donor_long.csv")

if (!file.exists(path.expand(long_csv))) {
  stop("找不到输入文件：", long_csv)
}


df <- read_csv(path.expand(long_csv), show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = c("Young", "Aged")),
    boundary = factor(boundary, levels = c("core_0_17", "ext_0_14_17")),
    module = factor(module, levels = c("tan", "greenyellow", "green"))
  )


effect_summary <- df %>%
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
  arrange(module, boundary)

write_csv(effect_summary, file.path(base_dir, "tables", "11_core_vs_ext_effect_summary.csv"))

attenuation_summary <- effect_summary %>%
  select(boundary, module, aged_minus_young) %>%
  pivot_wider(names_from = boundary, values_from = aged_minus_young) %>%
  mutate(
    attenuation_ext_minus_core = ext_0_14_17 - core_0_17
  ) %>%
  arrange(module)

write_csv(attenuation_summary, file.path(base_dir, "tables", "12_boundary_attenuation_summary.csv"))

set.seed(12345)

boot_effect_once <- function(sub_df) {
  young_vals <- sub_df %>% filter(group == "Young") %>% pull(score_mean)
  aged_vals  <- sub_df %>% filter(group == "Aged") %>% pull(score_mean)

  boot_young <- sample(young_vals, size = length(young_vals), replace = TRUE)
  boot_aged  <- sample(aged_vals,  size = length(aged_vals),  replace = TRUE)

  mean(boot_aged) - mean(boot_young)
}

boot_effect_ci <- function(sub_df, n_boot = 5000) {
  vals <- replicate(n_boot, boot_effect_once(sub_df))
  tibble(
    boot_mean = mean(vals),
    ci_low = unname(quantile(vals, 0.025)),
    ci_high = unname(quantile(vals, 0.975))
  )
}

boot_boundary_effects <- df %>%
  group_by(boundary, module) %>%
  group_modify(~ boot_effect_ci(.x, n_boot = 5000)) %>%
  ungroup() %>%
  arrange(module, boundary)

write_csv(boot_boundary_effects, file.path(base_dir, "tables", "13_bootstrap_boundary_effect_cis.csv"))

boot_attenuation_once <- function(sub_df_core, sub_df_ext) {
  eff_core <- boot_effect_once(sub_df_core)
  eff_ext  <- boot_effect_once(sub_df_ext)
  eff_ext - eff_core
}

boot_attenuation_ci <- function(df_core, df_ext, n_boot = 5000) {
  vals <- replicate(n_boot, boot_attenuation_once(df_core, df_ext))
  tibble(
    boot_mean = mean(vals),
    ci_low = unname(quantile(vals, 0.025)),
    ci_high = unname(quantile(vals, 0.975))
  )
}

boot_attenuation <- map_dfr(levels(df$module), function(m) {
  sub_core <- df %>% filter(module == m, boundary == "core_0_17")
  sub_ext  <- df %>% filter(module == m, boundary == "ext_0_14_17")

  out <- boot_attenuation_ci(sub_core, sub_ext, n_boot = 5000)
  out$module <- m
  out
}) %>%
  relocate(module)

write_csv(boot_attenuation, file.path(base_dir, "tables", "14_bootstrap_boundary_attenuation_cis.csv"))

final_effect_table <- effect_summary %>%
  left_join(boot_boundary_effects, by = c("boundary", "module")) %>%
  arrange(module, boundary)

write_csv(final_effect_table, file.path(base_dir, "tables", "15_final_core_vs_ext_effect_table.csv"))

final_attenuation_table <- attenuation_summary %>%
  left_join(boot_attenuation, by = "module") %>%
  arrange(module)

write_csv(final_attenuation_table, file.path(base_dir, "tables", "16_final_boundary_attenuation_table.csv"))

p1 <- ggplot(df, aes(x = group, y = score_mean)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
  facet_grid(module ~ boundary, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Core vs extended boundary: donor-level module scores",
    x = NULL,
    y = "Mean module score per donor"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "4_core_vs_ext_donor_level_boxplot.png"),
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)

plot_effect_df <- final_effect_table %>%
  mutate(boundary = factor(boundary, levels = c("core_0_17", "ext_0_14_17")))

p2 <- ggplot(plot_effect_df, aes(x = boundary, y = aged_minus_young)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.8) +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Boundary-specific module effects (Aged minus Young)",
    x = NULL,
    y = "Effect size"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "5_boundary_specific_effect_sizes.png"),
  plot = p2,
  width = 7,
  height = 8,
  dpi = 300
)

p3 <- ggplot(final_attenuation_table, aes(x = module, y = attenuation_ext_minus_core)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = "Boundary attenuation (ext effect - core effect)",
    x = NULL,
    y = "Attenuation"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(base_dir, "figures", "6_boundary_attenuation.png"),
  plot = p3,
  width = 6,
  height = 4.8,
  dpi = 300
)





