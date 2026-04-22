#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(tibble)
})

source(file.path("R", "model_stage_repeatability_helpers.R"))

derived_dir <- file.path("data", "derived")
ensure_dir(derived_dir)

stage0_summary_path <- file.path(derived_dir, "model_stage0_assay_repeatability_summary.csv")
stage0_group_path <- file.path(derived_dir, "model_stage0_assay_repeatability_repeat_groups.csv")
museum_long_path <- file.path(derived_dir, "museum_measurements_long.csv")
museum_region_path <- file.path(derived_dir, "museum_specimen_region_isotopes_long.csv")
live_tissue_path <- file.path(derived_dir, "live_cn_tissue_summary.csv")
live_pair_path <- file.path(derived_dir, "live_cn_paired_by_ring.csv")
screening_path <- file.path(derived_dir, "live_screening_ready_paired_labelled.csv")

read_homogenization_method_comparison <- function(path) {
  read_excel(
    path,
    sheet = "Sample results_C&N",
    col_types = "text",
    skip = 9,
    col_names = TRUE
  ) %>%
    dplyr::select(1, 2, 10:16, 18:24, 28:35, 40:46, 48) %>%
    rename_with(
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all(" ", "_") %>%
        str_replace_all("__", "_") %>%
        str_replace_all("___", "_") %>%
        str_replace_all("%", "_perc_") %>%
        str_replace_all(":", "") %>%
        str_replace_all("#", "")
    ) %>%
    rename(
      normalised_d15n_____A = normalised_d15n...12,
      normalised_d15n_____B = normalised_d15n...20,
      normalised_d13c_____C = normalised_d13c...30,
      normalised_d13c_____D = normalised_d13c...42
    ) %>%
    pivot_longer(
      cols = c(
        normalised_d15n_____A,
        normalised_d15n_____B,
        normalised_d13c_____C,
        normalised_d13c_____D
      ),
      names_to = c("metric", "variable"),
      names_sep = "_____"
    ) %>%
    mutate(
      method_code = ifelse(variable %in% c("A", "C"), "B", "A"),
      isotope = ifelse(variable %in% c("A", "B"), "normalised_d15n", "normalised_d13c"),
      value = suppressWarnings(as.numeric(value)),
      specimen_id = str_extract(identifier_1, "^[^ ]+"),
      feather_id = str_trim(str_replace(identifier_1, "^[^ ]+ ", ""))
    ) %>%
    filter(
      metric %in% c("normalised_d15n", "normalised_d13c"),
      !is.na(value)
    ) %>%
    transmute(
      specimen_id,
      feather_id,
      isotope,
      method_code,
      method_label = case_when(
        method_code == "A" ~ "whole_half_vane_homogenate",
        method_code == "B" ~ "vane_subsample_homogenate",
        TRUE ~ NA_character_
      ),
      value
    )
}

museum_long <- read_csv(museum_long_path, show_col_types = FALSE)
museum_region <- read_csv(museum_region_path, show_col_types = FALSE)
live_tissue <- read_csv(live_tissue_path, show_col_types = FALSE)
live_pair <- read_csv(live_pair_path, show_col_types = FALSE)
live_screening <- read_csv(screening_path, show_col_types = FALSE)
stage0_summary <- if (file.exists(stage0_summary_path)) {
  read_csv(stage0_summary_path, show_col_types = FALSE)
} else {
  tibble()
}
stage0_groups <- if (file.exists(stage0_group_path)) {
  read_csv(stage0_group_path, show_col_types = FALSE)
} else {
  tibble()
}

triplicate_rows <- museum_long %>%
  filter(
    feather_type == "Breast",
    region == "homogenate",
    isotope %in% c("normalised_d13c", "normalised_d15n", "normalised_d2h"),
    !is.na(feather_replicate),
    feather_replicate != "",
    !is.na(value)
  )

triplicate_detail <- triplicate_rows %>%
  group_by(specimen_id, feather_replicate, isotope) %>%
  summarise(
    n_measurements = n(),
    value = mean(value),
    geo_region = dplyr::first(geo_region),
    source_measurement_ids = paste(unique(measurement_id), collapse = "|"),
    source_branches = paste(unique(source_branch), collapse = "|"),
    source_files = paste(unique(source_file), collapse = "|"),
    source_sheets = paste(unique(source_sheet), collapse = "|"),
    evidence_design = "distinct_breast_feathers_same_bird",
    .groups = "drop"
  ) %>%
  arrange(isotope, specimen_id, feather_replicate)

write_csv(
  triplicate_detail,
  file.path(derived_dir, "model_stage1_breast_feather_repeatability_detail.csv")
)

triplicate_summary <- triplicate_detail %>%
  group_by(isotope) %>%
  group_modify(
    ~ bind_cols(
      pooled_within_summary(.x, "specimen_id", "value"),
      one_way_random_effects_metrics(.x, "specimen_id", "value")
    )
  ) %>%
  ungroup() %>%
  mutate(
    feather_type = "Breast",
    n_specimens = vapply(
      isotope,
      function(iso) dplyr::n_distinct(triplicate_detail$specimen_id[triplicate_detail$isotope == iso]),
      integer(1)
    ),
    n_replicate_feathers = vapply(
      isotope,
      function(iso) sum(triplicate_detail$isotope == iso),
      integer(1)
    ),
    source_branches = vapply(
      isotope,
      function(iso) paste(unique(triplicate_detail$source_branches[triplicate_detail$isotope == iso]), collapse = "|"),
      character(1)
    ),
    source_files = vapply(
      isotope,
      function(iso) paste(unique(triplicate_detail$source_files[triplicate_detail$isotope == iso]), collapse = "|"),
      character(1)
    ),
    source_sheets = vapply(
      isotope,
      function(iso) paste(unique(triplicate_detail$source_sheets[triplicate_detail$isotope == iso]), collapse = "|"),
      character(1)
    ),
    evidence_strength = mapply(classify_evidence_strength, n_groups, total_within_df),
    evidence_stream = "within_individual_repeated_breast_feathers",
    evidence_design = "distinct_breast_feathers_same_bird"
  ) %>%
  select(
    feather_type,
    isotope,
    evidence_stream,
    evidence_design,
    source_branches,
    source_files,
    source_sheets,
    n_specimens,
    n_replicate_feathers,
    everything()
  ) %>%
  arrange(isotope)

write_csv(
  triplicate_summary,
  file.path(derived_dir, "model_stage1_breast_feather_repeatability_summary.csv")
)

homogenization_detail <- read_homogenization_method_comparison(
  file.path("data", "Master_SI_Replicate test Feathers_N_9_LEH_Jan25_v1_SB.xls")
) %>%
  arrange(isotope, specimen_id, method_code)

write_csv(
  homogenization_detail,
  file.path(derived_dir, "model_stage1_homogenization_method_detail.csv")
)

homogenization_summary <- homogenization_detail %>%
  select(specimen_id, feather_id, isotope, method_code, value) %>%
  pivot_wider(names_from = method_code, values_from = value, names_prefix = "method_") %>%
  filter(!is.na(method_A), !is.na(method_B)) %>%
  mutate(
    diff_subsample_minus_whole = method_B - method_A,
    mean_method_value = (method_A + method_B) / 2
  ) %>%
  group_by(isotope) %>%
  summarise(
    n_pairs = n(),
    mean_diff_subsample_minus_whole = mean(diff_subsample_minus_whole),
    sd_diff_subsample_minus_whole = safe_sd(diff_subsample_minus_whole),
    mean_abs_diff = mean(abs(diff_subsample_minus_whole)),
    max_abs_diff = max(abs(diff_subsample_minus_whole)),
    rmse_diff = sqrt(mean(diff_subsample_minus_whole^2)),
    pearson_r = ifelse(n() <= 1, NA_real_, cor(method_A, method_B)),
    between_specimen_sd = safe_sd(mean_method_value),
    rmse_to_between_specimen_sd = rmse_diff / between_specimen_sd,
    .groups = "drop"
  ) %>%
  mutate(
    evidence_stream = "same_feather_homogenization_method_comparison",
    method_A_label = "whole_half_vane_homogenate",
    method_B_label = "vane_subsample_homogenate"
  ) %>%
  arrange(isotope)

write_csv(
  homogenization_summary,
  file.path(derived_dir, "model_stage1_homogenization_method_summary.csv")
)

region_detail <- museum_region %>%
  filter(
    isotope %in% c("normalised_d13c", "normalised_d15n"),
    feather_type %in% c("Breast", "Primary")
  ) %>%
  group_by(specimen_id, feather_type, isotope) %>%
  filter(n() >= 2) %>%
  summarise(
    n_regions = n(),
    regions_present = paste(region, collapse = "|"),
    values_by_region = paste(paste(region, formatC(value_mean, format = "f", digits = 4), sep = ":"), collapse = "|"),
    region_range = safe_range(value_mean),
    region_sd = safe_sd(value_mean),
    homogenate_value = ifelse(any(region == "homogenate"), value_mean[region == "homogenate"][1], NA_real_),
    nonhomogenate_mean = ifelse(any(region != "homogenate"), mean(value_mean[region != "homogenate"]), NA_real_),
    homogenate_minus_nonhomogenate_mean = homogenate_value - nonhomogenate_mean,
    .groups = "drop"
  ) %>%
  arrange(feather_type, isotope, specimen_id)

write_csv(
  region_detail,
  file.path(derived_dir, "model_stage1_region_structure_detail.csv")
)

region_summary <- region_detail %>%
  group_by(feather_type, isotope) %>%
  summarise(
    n_specimens = n(),
    mean_region_range = mean(region_range, na.rm = TRUE),
    median_region_range = median(region_range, na.rm = TRUE),
    max_region_range = max(region_range, na.rm = TRUE),
    mean_region_sd = mean(region_sd, na.rm = TRUE),
    mean_abs_homogenate_minus_nonhomogenate = mean(abs(homogenate_minus_nonhomogenate_mean), na.rm = TRUE),
    max_abs_homogenate_minus_nonhomogenate = max(abs(homogenate_minus_nonhomogenate_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(evidence_stream = "within_feather_region_structure") %>%
  arrange(feather_type, isotope)

write_csv(
  region_summary,
  file.path(derived_dir, "model_stage1_region_structure_summary.csv")
)

live_tissue_coverage <- live_tissue %>%
  group_by(feather_type) %>%
  summarise(
    n_tissues = n(),
    n_with_d13c = sum(!is.na(normalised_d13c)),
    n_with_d15n = sum(!is.na(normalised_d15n)),
    n_with_d2h = sum(!is.na(normalised_d2h)),
    n_multi_assay = sum(n_assays > 1, na.rm = TRUE),
    median_n_assays = median(n_assays, na.rm = TRUE),
    max_n_assays = max(n_assays, na.rm = TRUE),
    median_n_unique_samples = median(n_unique_samples, na.rm = TRUE),
    max_n_unique_samples = max(n_unique_samples, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(feather_type)

write_csv(
  live_tissue_coverage,
  file.path(derived_dir, "model_stage1_live_tissue_coverage_summary.csv")
)

live_pair_coverage <- bind_rows(
  live_pair %>%
    count(status_known, has_complete_cnh_pair, name = "n_rings") %>%
    mutate(dataset = "live_cn_paired_by_ring"),
  live_screening %>%
    summarise(
      status_known = all(status_known),
      has_complete_cnh_pair = all(has_complete_cnh_pair),
      n_rings = n()
    ) %>%
    mutate(dataset = "live_screening_ready_paired_labelled")
) %>%
  relocate(dataset)

write_csv(
  live_pair_coverage,
  file.path(derived_dir, "model_stage1_live_paired_coverage_summary.csv")
)

stage_evidence_sources <- bind_rows(
  stage0_groups %>%
    transmute(
      stage = "Stage0",
      feather_type,
      isotope,
      evidence_stream = "technical_repeat_runs_same_sample_tissue",
      evidence_design = dplyr::coalesce(evidence_design, "technical_repeat_runs_same_sample_tissue"),
      source_branches,
      source_files,
      source_sheets,
      n_groups = 1L,
      n_rows = n_runs,
      n_specimens = NA_real_,
      notes = "Retained repeat assay runs from the live sample:tissue audit table."
    ) %>%
    group_by(stage, feather_type, isotope, evidence_stream, evidence_design, source_branches, source_files, source_sheets, notes) %>%
    summarise(
      n_groups = sum(n_groups),
      n_rows = sum(n_rows, na.rm = TRUE),
      n_specimens = NA_real_,
      .groups = "drop"
    ),
  triplicate_rows %>%
    group_by(feather_type, isotope, source_branch, source_file, source_sheet) %>%
    summarise(
      n_rows = n(),
      n_specimens = n_distinct(specimen_id),
      .groups = "drop"
    ) %>%
    transmute(
      stage = "Stage1",
      feather_type,
      isotope,
      evidence_stream = "within_individual_repeated_breast_feathers",
      evidence_design = "distinct_breast_feathers_same_bird",
      source_branches = source_branch,
      source_files = source_file,
      source_sheets = source_sheet,
      n_groups = n_specimens,
      n_rows,
      n_specimens,
      notes = "Separate breast feathers from the same museum bird; suitable for biological/tissue repeatability, not technical assay repeatability."
    ),
  homogenization_detail %>%
    group_by(isotope) %>%
    summarise(
      n_rows = n(),
      n_specimens = n_distinct(specimen_id),
      .groups = "drop"
    ) %>%
    transmute(
      stage = "Stage1",
      feather_type = "Breast",
      isotope,
      evidence_stream = "same_feather_homogenization_method_comparison",
      evidence_design = "paired_preparation_methods_same_feather",
      source_branches = "museum_homogenization_comparison_cn",
      source_files = "data/Master_SI_Replicate test Feathers_N_9_LEH_Jan25_v1_SB.xls",
      source_sheets = "Sample results_C&N",
      n_groups = n_specimens,
      n_rows,
      n_specimens,
      notes = "Paired preparation methods measured on the same feather."
    )
) %>%
  arrange(stage, feather_type, isotope, evidence_stream, source_files)

write_csv(
  stage_evidence_sources,
  file.path(derived_dir, "model_stage_evidence_sources.csv")
)

adequacy_summary <- tibble(
  feather_type = c("Breast", "Breast", "Breast", "Primary", "Primary", "Primary"),
  isotope = c("normalised_d13c", "normalised_d15n", "normalised_d2h", "normalised_d13c", "normalised_d15n", "normalised_d2h")
) %>%
  left_join(
    stage0_summary %>%
      select(
        feather_type,
        isotope,
        stage0_direct_repeat_groups = n_groups,
        stage0_direct_pooled_sd = pooled_within_sd,
        stage0_evidence_strength = evidence_strength,
        stage0_evidence_design = evidence_design,
        stage0_source_files = source_files
      ),
    by = c("feather_type", "isotope")
  ) %>%
  left_join(
    triplicate_summary %>%
      select(
        isotope,
        stage1_breast_repeat_specimens_same_isotope = n_specimens,
        stage1_breast_repeat_icc_same_isotope = repeatability_icc,
        stage1_breast_repeat_evidence_design = evidence_design,
        stage1_breast_repeat_source_files = source_files
      ),
    by = "isotope"
  ) %>%
  left_join(
    homogenization_summary %>%
      select(
        isotope,
        homogenization_n_pairs = n_pairs,
        homogenization_rmse_to_between_specimen_sd = rmse_to_between_specimen_sd
      ),
    by = "isotope"
  ) %>%
  left_join(
    region_summary %>%
      select(
        feather_type,
        isotope,
        region_structure_n_specimens = n_specimens,
        region_structure_mean_range = mean_region_range
      ),
    by = c("feather_type", "isotope")
  ) %>%
  left_join(
    live_tissue_coverage %>%
      transmute(
        feather_type,
        live_n_tissues = n_tissues,
        live_n_with_d2h = n_with_d2h,
        live_n_multi_assay = n_multi_assay
      ),
    by = "feather_type"
  ) %>%
  mutate(
    stage1_breast_repeat_specimens_same_isotope = ifelse(
      feather_type == "Breast",
      stage1_breast_repeat_specimens_same_isotope,
      NA_real_
    ),
    stage1_breast_repeat_icc_same_isotope = ifelse(
      feather_type == "Breast",
      stage1_breast_repeat_icc_same_isotope,
      NA_real_
    ),
    stage1_breast_repeat_evidence_design = ifelse(
      feather_type == "Breast",
      stage1_breast_repeat_evidence_design,
      NA_character_
    ),
    stage1_breast_repeat_source_files = ifelse(
      feather_type == "Breast",
      stage1_breast_repeat_source_files,
      NA_character_
    ),
    biological_repeatability_evidence = case_when(
      feather_type == "Breast" & isotope %in% c("normalised_d13c", "normalised_d15n", "normalised_d2h") ~
        "direct_repeated_breast_feathers",
      TRUE ~ "not_directly_available"
    ),
    homogenization_evidence = case_when(
      isotope %in% c("normalised_d13c", "normalised_d15n") ~ "same_feather_method_comparison_available",
      TRUE ~ "not_directly_available"
    ),
    adequacy_for_stage2 = case_when(
      feather_type == "Breast" & isotope %in% c("normalised_d13c", "normalised_d15n") ~
        "adequate",
      feather_type == "Breast" & isotope == "normalised_d2h" ~
        "adequate_with_technical_caution",
      feather_type == "Primary" & isotope %in% c("normalised_d13c", "normalised_d15n") ~
        "adequate_with_limited_biological_evidence",
      feather_type == "Primary" & isotope == "normalised_d2h" ~
        "provisionally_adequate_with_high_caution",
      TRUE ~ "undetermined"
    ),
    adequacy_notes = case_when(
      feather_type == "Breast" & isotope %in% c("normalised_d13c", "normalised_d15n") ~
        "Breast C/N have direct repeated-feather evidence and same-feather homogenization-method evidence.",
      feather_type == "Breast" & isotope == "normalised_d2h" ~
        "Breast H has Stage 1 repeated-breast-feather evidence from distinct feathers sampled from the same museum birds and strong live coverage, but only two live assay repeat pairs and no direct homogenization-method comparison.",
      feather_type == "Primary" & isotope %in% c("normalised_d13c", "normalised_d15n") ~
        "Primary C/N are well covered in the live paired table, but direct repeated-feather evidence is sparse and region-structure evidence is limited to two specimens.",
      feather_type == "Primary" & isotope == "normalised_d2h" ~
        "Primary H is available in the live paired table, but there are no direct primary H assay repeats and no direct biological repeatability dataset.",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(feather_type, isotope)

write_csv(
  adequacy_summary,
  file.path(derived_dir, "model_stage1_tissue_summary_adequacy.csv")
)

if (nrow(triplicate_detail) > 0) {
  p_triplicate <- ggplot(
    triplicate_detail,
    aes(x = feather_replicate, y = value, group = specimen_id, colour = specimen_id)
  ) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    facet_wrap(~ isotope, scales = "free_y") +
    labs(
      title = "Stage 1 breast repeated-feather diagnostics",
      x = "Feather replicate",
      y = "Isotope value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(derived_dir, "model_stage1_breast_feather_repeatability.png"),
    plot = p_triplicate,
    width = 10,
    height = 6,
    dpi = 200
  )
}

if (nrow(homogenization_detail) > 0) {
  p_homogenization <- homogenization_detail %>%
    mutate(method_label = factor(
      method_label,
      levels = c("whole_half_vane_homogenate", "vane_subsample_homogenate")
    )) %>%
    ggplot(
      aes(x = method_label, y = value, group = specimen_id)
    ) +
    geom_line(alpha = 0.6, colour = "#7A5C61") +
    geom_point(size = 2.2, colour = "#2C3E50") +
    facet_wrap(~ isotope, scales = "free_y") +
    labs(
      title = "Stage 1 homogenization-method comparison",
      x = "Method",
      y = "Isotope value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(derived_dir, "model_stage1_homogenization_method_comparison.png"),
    plot = p_homogenization,
    width = 10,
    height = 6,
    dpi = 200
  )
}

summary_lines <- c(
  "# Stage 0 and Stage 1 Summary",
  "",
  "## Stage 0: assay repeatability",
  ""
)

if (nrow(stage0_summary) > 0) {
  stage0_lines <- apply(stage0_summary, 1, function(row) {
    paste0(
      "- `", row[["feather_type"]], " ", row[["isotope"]], "`: ",
      row[["n_groups"]], " repeat groups, pooled SD ",
      format_num(as.numeric(row[["pooled_within_sd"]])), ", ICC ",
      format_num(as.numeric(row[["repeatability_icc"]])), ". ",
      row[["direct_notes"]]
    )
  })
  summary_lines <- c(summary_lines, stage0_lines, "", "## Stage 1: biological repeatability and homogenization", "")
}

triplicate_lines <- apply(triplicate_summary, 1, function(row) {
  paste0(
    "- Breast `", row[["isotope"]], "` repeated feathers: ",
    row[["n_specimens"]], " specimens, pooled within-specimen SD ",
    format_num(as.numeric(row[["pooled_within_sd"]])), ", ICC ",
    format_num(as.numeric(row[["repeatability_icc"]])), "."
  )
})

homogenization_lines <- apply(homogenization_summary, 1, function(row) {
  paste0(
    "- Homogenization method comparison `", row[["isotope"]], "`: ",
    row[["n_pairs"]], " paired feathers, mean subsample-minus-whole difference ",
    format_num(as.numeric(row[["mean_diff_subsample_minus_whole"]])), ", RMSE/between-specimen SD ",
    format_num(as.numeric(row[["rmse_to_between_specimen_sd"]])), "."
  )
})

region_lines <- apply(region_summary, 1, function(row) {
  paste0(
    "- Region-structure check `", row[["feather_type"]], " ", row[["isotope"]], "`: ",
    row[["n_specimens"]], " specimens with >=2 regions, mean region range ",
    format_num(as.numeric(row[["mean_region_range"]])), "."
  )
})

adequacy_lines <- apply(adequacy_summary, 1, function(row) {
  paste0(
    "- `", row[["feather_type"]], " ", row[["isotope"]], "`: ",
    row[["adequacy_for_stage2"]], ". ",
    row[["adequacy_notes"]]
  )
})

summary_lines <- c(
  summary_lines,
  triplicate_lines,
  "",
  homogenization_lines,
  "",
  region_lines,
  "",
  "## Adequacy for later paired modelling",
  "",
  adequacy_lines,
  ""
)

writeLines(
  summary_lines,
  file.path(derived_dir, "model_stage0_stage1_summary.md")
)
