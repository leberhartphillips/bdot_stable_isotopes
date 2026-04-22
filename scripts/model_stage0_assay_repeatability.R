#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

collapse_unique_chr <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = "|")
}

input_path <- file.path(derived_dir, "live_cn_measurements_raw_qc_long.csv")

target_isotopes <- c("normalised_d13c", "normalised_d15n", "normalised_d2h")
target_tissues <- c("Breast", "Primary")

assay_raw <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(
    feather_type %in% target_tissues,
    isotope %in% target_isotopes,
    assay_qc_audit_retained %in% TRUE,
    !is.na(raw_value)
  )

repeat_rows <- assay_raw %>%
  filter(assay_qc_repeat_pair_detected %in% TRUE)

repeat_group_detail <- repeat_rows %>%
  group_by(
    feather_type,
    isotope,
    sample_id_base,
    ring,
    assay_qc_repeat_resolution_status,
    assay_qc_model_value_method
  ) %>%
  summarise(
    source_branches = collapse_unique_chr(source_branch),
    source_files = collapse_unique_chr(source_file),
    source_sheets = collapse_unique_chr(source_sheet),
    n_runs = n(),
    sample_id_raws = paste(unique(sample_id_raw), collapse = "|"),
    source_sheet_rows = paste(unique(source_sheet_row), collapse = "|"),
    raw_values = paste(formatC(raw_value, format = "f", digits = 6), collapse = "|"),
    mean_raw_value = mean(raw_value),
    sd_raw_value = safe_sd(raw_value),
    range_raw_value = safe_range(raw_value),
    min_raw_value = min(raw_value),
    max_raw_value = max(raw_value),
    mean_abs_deviation = mean(abs(raw_value - mean(raw_value))),
    adjudicated_model_value = dplyr::first(value[!is.na(value)]),
    .groups = "drop"
  ) %>%
  mutate(
    evidence_design = "technical_repeat_runs_same_sample_tissue"
  ) %>%
  arrange(feather_type, isotope, sample_id_base)

write_csv(
  repeat_group_detail,
  file.path(derived_dir, "model_stage0_assay_repeatability_repeat_groups.csv")
)

direct_summary <- tidyr::expand_grid(
  feather_type = target_tissues,
  isotope = target_isotopes
) %>%
  rowwise() %>%
  do({
    feather_type_i <- .$feather_type
    isotope_i <- .$isotope

    subset_df <- repeat_rows %>%
      filter(feather_type == feather_type_i, isotope == isotope_i)

    pooled <- pooled_within_summary(subset_df, "sample_id_base", "raw_value")
    icc <- one_way_random_effects_metrics(subset_df, "sample_id_base", "raw_value")
    source_summary <- tibble(
      source_branches = collapse_unique_chr(subset_df$source_branch),
      source_files = collapse_unique_chr(subset_df$source_file),
      source_sheets = collapse_unique_chr(subset_df$source_sheet),
      evidence_design = if (nrow(subset_df) > 0) {
        "technical_repeat_runs_same_sample_tissue"
      } else {
        NA_character_
      }
    )

    bind_cols(
      tibble(
        feather_type = feather_type_i,
        isotope = isotope_i
      ),
      source_summary,
      pooled,
      icc
    )
  }) %>%
  ungroup() %>%
  mutate(
    direct_estimate_available = n_groups > 0,
    evidence_strength = mapply(classify_evidence_strength, n_groups, total_within_df),
    direct_notes = case_when(
      feather_type == "Primary" & isotope == "normalised_d2h" ~
        "No retained primary H repeat runs were available for direct estimation.",
      direct_estimate_available & evidence_strength == "weak" ~
        "Direct estimate is based on very few repeat groups and should be treated cautiously.",
      direct_estimate_available & evidence_strength == "moderate" ~
        "Direct estimate is usable but still based on a limited number of repeat groups.",
      direct_estimate_available ~
        "Direct estimate is supported by multiple retained repeat groups.",
      TRUE ~
        "No direct repeat-based estimate is available for this tissue-isotope combination."
    )
  ) %>%
  arrange(feather_type, isotope)

write_csv(
  direct_summary,
  file.path(derived_dir, "model_stage0_assay_repeatability_summary.csv")
)

breast_h_direct <- direct_summary %>%
  filter(feather_type == "Breast", isotope == "normalised_d2h")

variance_recommendation <- direct_summary %>%
  transmute(
    feather_type,
    isotope,
    direct_estimate_available,
    direct_pooled_variance = pooled_within_variance,
    direct_pooled_sd = pooled_within_sd,
    direct_repeat_groups = n_groups,
    direct_repeat_rows = n_rows,
    evidence_strength,
    working_variance = dplyr::case_when(
      direct_estimate_available ~ pooled_within_variance,
      feather_type == "Primary" & isotope == "normalised_d2h" &
        nrow(breast_h_direct) == 1L ~ breast_h_direct$pooled_within_variance,
      TRUE ~ NA_real_
    ),
    working_sd = dplyr::case_when(
      direct_estimate_available ~ pooled_within_sd,
      feather_type == "Primary" & isotope == "normalised_d2h" &
        nrow(breast_h_direct) == 1L ~ breast_h_direct$pooled_within_sd,
      TRUE ~ NA_real_
    ),
    working_source = dplyr::case_when(
      direct_estimate_available ~ "direct_same_tissue_isotope",
      feather_type == "Primary" & isotope == "normalised_d2h" &
        nrow(breast_h_direct) == 1L ~ "proxy_breast_normalised_d2h",
      TRUE ~ "none"
    ),
    working_notes = dplyr::case_when(
      direct_estimate_available ~ direct_notes,
      feather_type == "Primary" & isotope == "normalised_d2h" &
        nrow(breast_h_direct) == 1L ~
        "No direct primary H repeat estimate exists. Breast H is the only isotope-matched technical proxy and should be used only as a sensitivity/default working value.",
      TRUE ~ direct_notes
    )
  ) %>%
  arrange(feather_type, isotope)

write_csv(
  variance_recommendation,
  file.path(derived_dir, "model_stage0_pooled_assay_variance.csv")
)

if (nrow(repeat_rows) > 0) {
  repeat_plot <- repeat_rows %>%
    mutate(sample_id_base = factor(sample_id_base, levels = unique(repeat_group_detail$sample_id_base)))

  p <- ggplot(
    repeat_plot,
    aes(x = sample_id_base, y = raw_value, group = sample_id_base)
  ) +
    geom_line(alpha = 0.5, colour = "#4C6A92") +
    geom_point(size = 2.3, alpha = 0.9, colour = "#1B2838") +
    facet_grid(isotope ~ feather_type, scales = "free_y") +
    labs(
      title = "Stage 0 retained repeat-run diagnostics",
      x = "Latent sample",
      y = "Retained raw assay value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(derived_dir, "model_stage0_assay_repeatability_repeat_runs.png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 200
  )
}
