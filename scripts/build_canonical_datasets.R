#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(stringr)
  library(tibble)
})

source(file.path("R", "canonical_data_helpers.R"))

derived_dir <- file.path("data", "derived")
qc_dir <- file.path(derived_dir, "qc")

dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

write_csv_base <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE, na = "")
}

read_csv_base <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    check.names = FALSE
  )
}

qc_steps <- tibble()
qc_issues <- tibble()

log_step <- function(dataset, step, rows_in, rows_out, action, detail) {
  qc_steps <<- bind_rows(
    qc_steps,
    tibble(
      dataset = dataset,
      step = step,
      rows_in = as.integer(rows_in),
      rows_out = as.integer(rows_out),
      rows_changed = as.integer(rows_out - rows_in),
      action = action,
      detail = detail
    )
  )
}

log_issue <- function(dataset, issue_type, severity, key = NA_character_, n_records = 1L, detail) {
  qc_issues <<- bind_rows(
    qc_issues,
    tibble(
      dataset = dataset,
      issue_type = issue_type,
      severity = severity,
      key = as.character(key),
      n_records = as.integer(n_records),
      detail = detail
    )
  )
}

summarise_character_or_na <- function(x) {
  out <- first_non_missing(x)
  if (length(out) == 0) {
    return(NA_character_)
  }
  as.character(out)
}

summarise_integer_or_na <- function(x) {
  out <- suppressWarnings(as.integer(first_non_missing(as.character(x))))
  if (length(out) == 0) {
    return(NA_integer_)
  }
  out
}

safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x)
}

make_harmonization_flag <- function(...) {
  flags <- c(...)
  flags <- flags[!is.na(flags) & flags != ""]
  if (length(flags) == 0) {
    return(NA_character_)
  }
  paste(unique(flags), collapse = "|")
}

run_assay_qc_extract <- function() {
  status <- system2(
    "python3",
    args = file.path("scripts", "extract_assay_qc_adjudication.py")
  )

  if (!identical(status, 0L)) {
    stop("Failed to generate assay QC adjudication tables.")
  }
}

run_assay_qc_extract()

assay_qc_adjudication <- read_csv_base(
  file.path(derived_dir, "assay_qc_adjudication.csv")
) %>%
  mutate(
    source_sheet_row = as.integer(source_sheet_row),
    raw_value = parse_numeric(raw_value),
    adjudicated_value = parse_numeric(adjudicated_value),
    repeat_pair_detected = as.logical(repeat_pair_detected),
    sensitivity_selected_source_row = suppressWarnings(as.integer(sensitivity_selected_source_row)),
    sensitivity_selected_value = parse_numeric(sensitivity_selected_value),
    audit_retained = as.logical(audit_retained),
    excluded_by_lab_instruction = as.logical(excluded_by_lab_instruction),
    contributes_to_model_value = as.logical(contributes_to_model_value)
  )

assay_qc_affected_samples <- read_csv_base(
  file.path(derived_dir, "assay_qc_affected_samples.csv")
) %>%
  mutate(adjudicated_value = parse_numeric(adjudicated_value))

assay_qc_manual_adjudication <- read_csv_base(
  file.path(derived_dir, "assay_qc_manual_adjudication.csv")
) %>%
  mutate(
    sample_base_identifier = as.character(sample_base_identifier),
    tissue_feather_type = as.character(tissue_feather_type),
    isotope = as.character(isotope),
    reason = as.character(reason),
    notes = as.character(notes)
  )

assay_qc_sheet_rules <- read_csv_base(
  file.path(derived_dir, "assay_qc_sheet_rules.csv")
)

assay_qc_repeat_pairs <- read_csv_base(
  file.path(derived_dir, "assay_qc_repeat_pairs.csv")
)

assay_qc_repeat_pair_summary <- read_csv_base(
  file.path(derived_dir, "assay_qc_repeat_pair_summary.csv")
) %>%
  mutate(n_repeat_pairs = as.integer(n_repeat_pairs))

assay_qc_dup_keys <- assay_qc_adjudication %>%
  count(raw_file, sheet, source_sheet_row, isotope, sort = TRUE) %>%
  filter(n > 1)

if (nrow(assay_qc_dup_keys) > 0) {
  stop("Assay QC adjudication rows are not unique on raw_file + sheet + source_sheet_row + isotope.")
}

log_step(
  dataset = "assay_qc",
  step = "load_adjudication_tables",
  rows_in = nrow(assay_qc_adjudication),
  rows_out = nrow(assay_qc_adjudication),
  action = "read",
  detail = "Loaded workbook-parsed assay QC adjudication, affected-sample, manual-review, repeat-pair, repeat-pair-summary, and sheet-rule tables."
)

if (nrow(assay_qc_manual_adjudication) > 0) {
  log_issue(
    dataset = "assay_qc",
    issue_type = "pending_manual_review",
    severity = "warning",
    key = collapse_unique(
      paste(
        assay_qc_manual_adjudication$sample_identifier,
        assay_qc_manual_adjudication$isotope
      )
    ),
    n_records = nrow(assay_qc_manual_adjudication),
    detail = "Some repeated assay sample-isotope pairs remain unresolved and are flagged for manual review."
  )
}

apply_assay_qc_long <- function(df) {
  qc_lookup <- assay_qc_adjudication %>%
    transmute(
      source_file = raw_file,
      source_sheet = sheet,
      source_sheet_row,
      isotope,
      assay_qc_decision_id = decision_id,
      assay_qc_action = action,
      assay_qc_reason = reason,
      assay_qc_evidence_source = evidence_source,
      assay_qc_decision_class = decision_class,
      assay_qc_rule_basis = rule_basis,
      assay_qc_notes = notes,
      assay_qc_repeat_pair_detected = repeat_pair_detected,
      assay_qc_repeat_detection_sources = repeat_detection_sources,
      assay_qc_repeat_detection_pattern = repeat_detection_pattern,
      assay_qc_repeat_resolution_status = repeat_resolution_status,
      assay_qc_audit_retained = audit_retained,
      assay_qc_excluded_by_lab_instruction = excluded_by_lab_instruction,
      assay_qc_contributes_to_model_value = contributes_to_model_value,
      assay_qc_model_value_method = model_value_method,
      assay_qc_model_source_rows = model_source_rows,
      assay_qc_strict_row_action = strict_row_action,
      assay_qc_sensitivity_row_action = sensitivity_row_action,
      assay_qc_sensitivity_selected_source_row = sensitivity_selected_source_row,
      assay_qc_sensitivity_selected_value = sensitivity_selected_value,
      qc_adjudicated_value = adjudicated_value
    )

  out <- df %>%
    left_join(
      qc_lookup,
      by = c("source_file", "source_sheet", "source_sheet_row", "isotope")
    )

  out$assay_qc_action <- dplyr::coalesce(out$assay_qc_action, "keep")
  out$assay_qc_audit_retained <- dplyr::coalesce(
    out$assay_qc_audit_retained,
    out$assay_qc_action != "exclude"
  )
  out$assay_qc_excluded_by_lab_instruction <- dplyr::coalesce(
    out$assay_qc_excluded_by_lab_instruction,
    out$assay_qc_action == "exclude"
  )
  out$assay_qc_contributes_to_model_value <- dplyr::coalesce(
    out$assay_qc_contributes_to_model_value,
    out$assay_qc_action != "exclude"
  )
  out$assay_qc_model_value_method <- dplyr::coalesce(
    out$assay_qc_model_value_method,
    ifelse(out$assay_qc_action == "exclude", NA_character_, "single_raw")
  )
  out$assay_qc_model_source_rows <- dplyr::coalesce(
    out$assay_qc_model_source_rows,
    as.character(out$source_sheet_row)
  )
  out$assay_qc_repeat_pair_detected <- dplyr::coalesce(
    out$assay_qc_repeat_pair_detected,
    FALSE
  )
  out$assay_qc_strict_row_action <- dplyr::coalesce(
    out$assay_qc_strict_row_action,
    ifelse(out$assay_qc_action == "exclude", "exclude", "include")
  )
  out$assay_qc_sensitivity_row_action <- dplyr::coalesce(
    out$assay_qc_sensitivity_row_action,
    ifelse(out$assay_qc_action == "exclude", "exclude", "include")
  )
  out$value <- ifelse(
    out$assay_qc_action == "exclude",
    NA_real_,
    ifelse(!is.na(out$qc_adjudicated_value), out$qc_adjudicated_value, out$raw_value)
  )
  out$audit_value <- ifelse(out$assay_qc_audit_retained, out$raw_value, NA_real_)
  out$strict_value <- ifelse(
    out$assay_qc_strict_row_action %in% c("exclude", "exclude_unresolved_repeat_pair"),
    NA_real_,
    ifelse(!is.na(out$qc_adjudicated_value), out$qc_adjudicated_value, out$raw_value)
  )
  out$sensitivity_value <- ifelse(
    out$assay_qc_sensitivity_row_action %in% c(
      "exclude",
      "exclude_nonselected_repeat_row",
      "manual_review_only"
    ),
    NA_real_,
    ifelse(
      !is.na(out$assay_qc_sensitivity_selected_value) &
        out$assay_qc_sensitivity_row_action == "include_non_rpt_default",
      out$assay_qc_sensitivity_selected_value,
      ifelse(!is.na(out$qc_adjudicated_value), out$qc_adjudicated_value, out$raw_value)
    )
  )
  out$included_in_canonical <- !is.na(out$audit_value)
  out$included_in_strict <- !is.na(out$strict_value)
  out$included_in_sensitivity <- !is.na(out$sensitivity_value)
  out$assay_qc_explicit <- !is.na(out$assay_qc_decision_id)
  out$assay_qc_value_changed <- !is.na(out$raw_value) &
    !is.na(out$value) &
    abs(out$raw_value - out$value) > 1e-12
  out$assay_qc_value_excluded <- !is.na(out$raw_value) &
    is.na(out$value) &
    out$assay_qc_action == "exclude"
  out
}

manual_review_lookup <- assay_qc_manual_adjudication %>%
  transmute(
    sample_base_identifier = sample_base_identifier,
    tissue_feather_type = tissue_feather_type,
    isotope = isotope,
    assay_qc_manual_review_pending = TRUE,
    assay_qc_manual_review_reason = reason,
    assay_qc_manual_review_notes = notes
  )

museum_coord_checkpoint <- readRDS(file.path("data", "CNOH_museum_samples_201125.rds")) %>%
  ungroup() %>%
  mutate(
    specimen_id = standardize_museum_specimen_id(Specimen_ID),
    feather_type = standardize_feather_type(Feather_type)
  )

coord_conflicts <- museum_coord_checkpoint %>%
  group_by(specimen_id, feather_type) %>%
  summarise(
    n_lon = distinct_non_missing_n(lon),
    n_lat = distinct_non_missing_n(lat),
    n_geo = distinct_non_missing_n(geo_region),
    .groups = "drop"
  ) %>%
  filter(n_lon > 1 | n_lat > 1 | n_geo > 1)

if (nrow(coord_conflicts) > 0) {
  log_issue(
    dataset = "museum_coordinates",
    issue_type = "conflicting_checkpoint_coordinates",
    severity = "error",
    key = collapse_unique(paste(coord_conflicts$specimen_id, coord_conflicts$feather_type)),
    n_records = nrow(coord_conflicts),
    detail = "Historical coordinate checkpoint has conflicting values for specimen_id + feather_type."
  )
  stop("Conflicting coordinates found in historical checkpoint.")
}

museum_coord_lookup <- museum_coord_checkpoint %>%
  group_by(specimen_id, feather_type) %>%
  summarise(
    coordinate_row_ids = paste0("museum_coord_", dplyr::cur_group_id()),
    lon = safe_mean(as.numeric(lon)),
    lat = safe_mean(as.numeric(lat)),
    geo_region_checkpoint = summarise_character_or_na(geo_region),
    geo_region_bin_checkpoint = summarise_integer_or_na(geo_region_bin),
    coordinate_source = "historical_checkpoint:data/CNOH_museum_samples_201125.rds",
    coordinate_collection_date = summarise_character_or_na(Collection_date),
    coordinate_month = summarise_integer_or_na(Month),
    coordinate_year = suppressWarnings(as.integer(safe_mean(as.numeric(Year)))),
    coordinate_specific_collection_location = summarise_character_or_na(Specific_collection_location),
    .groups = "drop"
  )

write_csv_base(
  museum_coord_lookup,
  file.path(derived_dir, "museum_coordinate_lookup_checkpoint.csv")
)

museum_inventory_raw <- read_excel(
  path = file.path("data", "MASTER_museum_ALL SAMPLES INVENTORY_banded_dotterels_v4.xlsx"),
  col_types = "text"
) %>%
  mutate(
    metadata_row_id = sprintf("museum_meta_%03d", row_number()),
    source_file = "data/MASTER_museum_ALL SAMPLES INVENTORY_banded_dotterels_v4.xlsx",
    source_sheet = "Sheet1",
    source_sheet_row = excel_data_row(row_number(), skip = 0),
    specimen_id_raw = normalize_space(Specimen_ID),
    feather_type_raw = normalize_space(Feather_type),
    specimen_id = standardize_museum_specimen_id(Specimen_ID),
    feather_type = dplyr::coalesce(
      standardize_feather_type(Feather_type),
      standardize_feather_type(Specimen_ID)
    ),
    collection_date_raw = normalize_space(Collection_date),
    collection_month = extract_month_bdot(collection_date_raw),
    collection_year = extract_year_bdot(collection_date_raw),
    general_collection_location = normalize_space(General_collection_location),
    specific_collection_location = normalize_space(Specific_collection_location),
    location_label = ifelse(
      is.na(specific_collection_location),
      general_collection_location,
      paste(specific_collection_location, general_collection_location, sep = ", ")
    ),
    geo_region = derive_geo_region(general_collection_location),
    geo_region_bin = case_when(
      geo_region == "Australia" ~ 1L,
      geo_region == "New Zealand" ~ 0L,
      TRUE ~ NA_integer_
    ),
    prefix_removed = case_when(
      str_detect(specimen_id_raw, regex("^AJB[[:space:]]+", ignore_case = TRUE)) ~ "AJB",
      str_detect(specimen_id_raw, regex("^WAM[[:space:]]+", ignore_case = TRUE)) ~ "WAM",
      TRUE ~ NA_character_
    ),
    feather_type_inferred = is.na(standardize_feather_type(feather_type_raw)) &
      !is.na(standardize_feather_type(specimen_id_raw)),
    harmonization_flags = case_when(
      !is.na(prefix_removed) & feather_type_inferred ~ paste0(
        "removed_prefix_", prefix_removed, "|feather_type_inferred_from_specimen_id"
      ),
      !is.na(prefix_removed) ~ paste0("removed_prefix_", prefix_removed),
      feather_type_inferred ~ "feather_type_inferred_from_specimen_id",
      TRUE ~ NA_character_
    )
  )

log_step(
  dataset = "museum_metadata",
  step = "read_inventory",
  rows_in = nrow(museum_inventory_raw),
  rows_out = nrow(museum_inventory_raw),
  action = "read",
  detail = "Loaded museum inventory workbook."
)

log_issue(
  dataset = "museum_metadata",
  issue_type = "harmonization_summary",
  severity = "info",
  n_records = sum(!is.na(museum_inventory_raw$prefix_removed)),
  detail = paste0(
    "Removed museum prefixes from inventory specimen IDs: AJB=",
    sum(museum_inventory_raw$prefix_removed == "AJB", na.rm = TRUE),
    ", WAM=",
    sum(museum_inventory_raw$prefix_removed == "WAM", na.rm = TRUE),
    "."
  )
)

log_issue(
  dataset = "museum_metadata",
  issue_type = "harmonization_summary",
  severity = "info",
  n_records = sum(museum_inventory_raw$feather_type_inferred, na.rm = TRUE),
  detail = paste0(
    "Feather type inferred from Specimen_ID for ",
    sum(museum_inventory_raw$feather_type_inferred, na.rm = TRUE),
    " inventory rows where Feather_type column was missing."
  )
)

museum_inventory_invalid <- museum_inventory_raw %>%
  filter(
    is.na(specimen_id) |
      specimen_id %in% c("", "?") |
      is.na(feather_type)
  )

if (nrow(museum_inventory_invalid) > 0) {
  log_issue(
    dataset = "museum_metadata",
    issue_type = "invalid_inventory_key",
    severity = "warning",
    key = collapse_unique(museum_inventory_invalid$specimen_id_raw),
    n_records = nrow(museum_inventory_invalid),
    detail = "Inventory rows with missing or placeholder specimen_id / feather_type were excluded from canonical metadata."
  )
}

museum_inventory_valid <- museum_inventory_raw %>%
  filter(
    !is.na(specimen_id),
    specimen_id != "",
    specimen_id != "?",
    !is.na(feather_type)
  )

log_step(
  dataset = "museum_metadata",
  step = "drop_invalid_inventory_keys",
  rows_in = nrow(museum_inventory_raw),
  rows_out = nrow(museum_inventory_valid),
  action = "filter",
  detail = "Excluded inventory rows with missing or placeholder canonical keys."
)

museum_metadata_conflicts <- museum_inventory_valid %>%
  group_by(specimen_id, feather_type) %>%
  summarise(
    n_rows = n(),
    n_general_location = distinct_non_missing_n(general_collection_location),
    n_specific_location = distinct_non_missing_n(specific_collection_location),
    n_collection_date = distinct_non_missing_n(collection_date_raw),
    n_sex = distinct_non_missing_n(Sex),
    .groups = "drop"
  ) %>%
  filter(
    n_general_location > 1 |
      n_specific_location > 1 |
      n_collection_date > 1 |
      n_sex > 1
  )

if (nrow(museum_metadata_conflicts) > 0) {
  log_issue(
    dataset = "museum_metadata",
    issue_type = "conflicting_metadata_labels",
    severity = "error",
    key = collapse_unique(paste(museum_metadata_conflicts$specimen_id, museum_metadata_conflicts$feather_type)),
    n_records = nrow(museum_metadata_conflicts),
    detail = "Canonical museum metadata has conflicting labels for at least one specimen_id + feather_type key."
  )
  stop("Conflicting labels found in museum metadata.")
}

museum_metadata_canonical <- museum_inventory_valid %>%
  group_by(specimen_id, feather_type) %>%
  summarise(
    source_metadata_row_ids = collapse_unique(metadata_row_id),
    source_metadata_row_count = n(),
    source_sheet_rows = collapse_unique(source_sheet_row),
    museum = summarise_character_or_na(Museum),
    species = summarise_character_or_na(Species),
    specimen_id_raw_examples = collapse_unique(specimen_id_raw),
    feather_type_raw_examples = collapse_unique(feather_type_raw),
    general_collection_location = summarise_character_or_na(general_collection_location),
    specific_collection_location = summarise_character_or_na(specific_collection_location),
    location_label = summarise_character_or_na(location_label),
    collection_date_raw = summarise_character_or_na(collection_date_raw),
    collection_month = summarise_integer_or_na(collection_month),
    collection_year = summarise_integer_or_na(collection_year),
    sex = summarise_character_or_na(Sex),
    notes = summarise_character_or_na(Notes),
    geo_region = summarise_character_or_na(geo_region),
    geo_region_bin = summarise_integer_or_na(geo_region_bin),
    metadata_harmonization_flags = collapse_unique(harmonization_flags),
    .groups = "drop"
  ) %>%
  left_join(museum_coord_lookup, by = c("specimen_id", "feather_type")) %>%
  mutate(
    coordinate_match_found = !is.na(lon) & !is.na(lat),
    coordinate_specific_location_matches_metadata = case_when(
      !coordinate_match_found ~ NA,
      is.na(specific_collection_location) | is.na(coordinate_specific_collection_location) ~ NA,
      TRUE ~ specific_collection_location == coordinate_specific_collection_location
    ),
    coordinate_collection_date_matches_metadata = case_when(
      !coordinate_match_found ~ NA,
      is.na(collection_date_raw) | is.na(coordinate_collection_date) ~ NA,
      TRUE ~ collection_date_raw == coordinate_collection_date
    )
  )

log_step(
  dataset = "museum_metadata",
  step = "collapse_to_canonical_metadata",
  rows_in = nrow(museum_inventory_valid),
  rows_out = nrow(museum_metadata_canonical),
  action = "deduplicate",
  detail = "Collapsed inventory rows to one canonical row per specimen_id + feather_type."
)

geo_region_mismatches <- museum_metadata_canonical %>%
  filter(
    !is.na(geo_region),
    !is.na(geo_region_checkpoint),
    geo_region != geo_region_checkpoint
  )

if (nrow(geo_region_mismatches) > 0) {
  log_issue(
    dataset = "museum_metadata",
    issue_type = "geo_region_mismatch_with_checkpoint",
    severity = "warning",
    key = collapse_unique(paste(geo_region_mismatches$specimen_id, geo_region_mismatches$feather_type)),
    n_records = nrow(geo_region_mismatches),
    detail = "Metadata-derived geo_region differs from historical coordinate checkpoint for some museum records."
  )
}

coordinate_location_mismatches <- museum_metadata_canonical %>%
  filter(coordinate_specific_location_matches_metadata %in% FALSE)

if (nrow(coordinate_location_mismatches) > 0) {
  log_issue(
    dataset = "museum_metadata",
    issue_type = "coordinate_location_mismatch_with_checkpoint",
    severity = "warning",
    key = collapse_unique(paste(coordinate_location_mismatches$specimen_id, coordinate_location_mismatches$feather_type)),
    n_records = nrow(coordinate_location_mismatches),
    detail = "Historical checkpoint coordinates disagree with current inventory specific_collection_location for some specimen_id + feather_type keys."
  )
}

coordinate_date_mismatches <- museum_metadata_canonical %>%
  filter(coordinate_collection_date_matches_metadata %in% FALSE)

if (nrow(coordinate_date_mismatches) > 0) {
  log_issue(
    dataset = "museum_metadata",
    issue_type = "coordinate_date_mismatch_with_checkpoint",
    severity = "warning",
    key = collapse_unique(paste(coordinate_date_mismatches$specimen_id, coordinate_date_mismatches$feather_type)),
    n_records = nrow(coordinate_date_mismatches),
    detail = "Historical checkpoint coordinates disagree with current inventory collection_date for some specimen_id + feather_type keys."
  )
}

write_csv_base(
  museum_metadata_canonical,
  file.path(derived_dir, "museum_metadata_canonical.csv")
)

read_triplicate_cn <- read_excel(
  path = file.path("data", "old file versions", "Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added.xls"),
  sheet = "Sample results_C&N",
  col_types = "text",
  skip = 10
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 10)
  )

log_step(
  dataset = "museum_triplicate_cn",
  step = "read_sheet",
  rows_in = nrow(read_triplicate_cn),
  rows_out = nrow(read_triplicate_cn),
  action = "read",
  detail = "Loaded triplicate C/N sheet."
)

triplicate_cn_long <- read_triplicate_cn %>%
  transmute(
    source_branch = "museum_triplicate_cn",
    source_file = "data/old file versions/Master_SI_Feathers_LEH_Sept24_v1_SB_v3_sample wts added.xls",
    source_sheet = "Sample results_C&N",
    source_sheet_row = source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`SIA Sample ID`),
    specimen_id_raw = str_extract(`Identifier 1`, "^[^ ]+"),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    feather_type_raw = str_trim(str_replace(`Identifier 1`, "^[^ ]+ ", "")) %>% str_replace("_.*", ""),
    feather_type = standardize_feather_type(feather_type_raw),
    feather_replicate = str_sub(normalize_space(`Identifier 1`), -1),
    region_raw = "homogenate",
    region = "homogenate",
    amount_mg = parse_numeric(Amount),
    normalised_d15n = parse_numeric(`normalised d15N`),
    normalised_d13c = parse_numeric(`normalised d13C`)
  ) %>%
  pivot_longer(
    cols = c(normalised_d15n, normalised_d13c),
    names_to = "isotope",
    values_to = "raw_value"
  )

log_step(
  dataset = "museum_triplicate_cn",
  step = "pivot_to_long",
  rows_in = nrow(read_triplicate_cn),
  rows_out = nrow(triplicate_cn_long),
  action = "reshape",
  detail = "Pivoted C/N triplicate rows to one row per isotope measurement."
)

triplicate_cn_long <- triplicate_cn_long %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf("museum_measurement_%04d", row_number())
  )

read_triplicate_h <- read_excel(
  path = file.path("data", "Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls"),
  sheet = "Sample results_H",
  col_types = "text",
  skip = 9
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 9)
  )

log_step(
  dataset = "museum_triplicate_h",
  step = "read_sheet",
  rows_in = nrow(read_triplicate_h),
  rows_out = nrow(read_triplicate_h),
  action = "read",
  detail = "Loaded triplicate H sheet."
)

triplicate_h_parsed <- read_triplicate_h %>%
  transmute(
    source_branch = "museum_triplicate_h",
    source_file = "data/Master_TCEA_H_SI_Batch 1_triplicate breast feathers_LEH_Nov25_SB.xls",
    source_sheet = "Sample results_H",
    source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`SIA ID`),
    specimen_id_raw = str_extract(`Identifier 1`, "^[^_]+"),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    specimen_id_from_sia = standardize_museum_specimen_id(`SIA ID`),
    feather_type = "Breast",
    feather_replicate = str_extract(`Identifier 1`, "(?<=_)."),
    region_raw = "homogenate",
    region = "homogenate",
    amount_mg = parse_numeric(`Amount (mg)`),
    isotope = "normalised_d2h",
    raw_value = parse_numeric(`Normalised d2H`)
  )

triplicate_h_mismatch <- triplicate_h_parsed %>%
  filter(
    !is.na(specimen_id),
    !is.na(specimen_id_from_sia),
    specimen_id != specimen_id_from_sia
  )

if (nrow(triplicate_h_mismatch) > 0) {
  log_issue(
    dataset = "museum_triplicate_h",
    issue_type = "identifier_vs_sia_id_mismatch",
    severity = "warning",
    key = collapse_unique(triplicate_h_mismatch$source_record_label),
    n_records = nrow(triplicate_h_mismatch),
    detail = "Triplicate H contains rows where Identifier 1 and SIA ID resolve to different specimen IDs. These rows remain visible in the raw QC table and must be excluded before metadata joins unless assay QC already excludes them."
  )
}

triplicate_h_long <- triplicate_h_parsed %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf("museum_measurement_%04d", nrow(triplicate_cn_long) + row_number())
  )

log_step(
  dataset = "museum_triplicate_h",
  step = "retain_raw_rows_pre_qc",
  rows_in = nrow(triplicate_h_parsed),
  rows_out = nrow(triplicate_h_long),
  action = "filter",
  detail = "Retained all non-missing H values before assay QC and metadata-join eligibility checks."
)

read_triplicate_o <- read_excel(
  path = file.path("data", "Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls"),
  sheet = "Sample results_O",
  col_types = "text",
  skip = 6
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 6)
  )

log_step(
  dataset = "museum_triplicate_o",
  step = "read_sheet",
  rows_in = nrow(read_triplicate_o),
  rows_out = nrow(read_triplicate_o),
  action = "read",
  detail = "Loaded O sheet used for both triplicate and batch 3/4 museum data."
)

triplicate_o_block <- read_triplicate_o %>%
  filter(read_row >= 2, read_row <= 16)

log_step(
  dataset = "museum_triplicate_o",
  step = "select_triplicate_block",
  rows_in = nrow(read_triplicate_o),
  rows_out = nrow(triplicate_o_block),
  action = "filter",
  detail = "Retained the first O-measurement block (rows 2-16 after read_excel) to mirror the current QMD triplicate section."
)

triplicate_o_filtered <- triplicate_o_block %>%
  filter(!is.na(`normalised d18O`))

log_step(
  dataset = "museum_triplicate_o",
  step = "drop_missing_d18o",
  rows_in = nrow(triplicate_o_block),
  rows_out = nrow(triplicate_o_filtered),
  action = "filter",
  detail = "Dropped triplicate O rows with missing d18O values."
)

triplicate_o_long <- triplicate_o_filtered %>%
  transmute(
    source_branch = "museum_triplicate_o",
    source_file = "data/Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls",
    source_sheet = "Sample results_O",
    source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`SIA ID`),
    specimen_id_raw = str_extract(`Identifier 1`, "^[^_]+"),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    feather_type = "Breast",
    feather_replicate = str_extract(`Identifier 1`, "(?<=_)."),
    region_raw = "homogenate",
    region = "homogenate",
    amount_mg = parse_numeric(Amount),
    isotope = "normalised_d18o",
    raw_value = parse_numeric(`normalised d18O`)
  ) %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf(
      "museum_measurement_%04d",
      nrow(triplicate_cn_long) + nrow(triplicate_h_long) + row_number()
    )
  )

log_step(
  dataset = "museum_triplicate_o",
  step = "retain_raw_rows_pre_qc",
  rows_in = nrow(triplicate_o_filtered),
  rows_out = nrow(triplicate_o_long),
  action = "filter",
  detail = "Retained all non-missing O values before assay QC adjudication."
)

read_batch_cn <- read_excel(
  path = file.path("data", "Master_EA_C_N_SI_Batch 3 vane_barb and Batch 4 Breast Feathers_LEH_Sept25_v1_SB.xls"),
  sheet = "Sample results_C&N",
  col_types = "text",
  skip = 9
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 9)
  )

log_step(
  dataset = "museum_batch_cn",
  step = "read_sheet",
  rows_in = nrow(read_batch_cn),
  rows_out = nrow(read_batch_cn),
  action = "read",
  detail = "Loaded batch 3/4 C/N sheet."
)

batch_cn_filtered <- read_batch_cn %>%
  filter(!is.na(`normalised d15N`))

log_step(
  dataset = "museum_batch_cn",
  step = "drop_missing_d15n_rows",
  rows_in = nrow(read_batch_cn),
  rows_out = nrow(batch_cn_filtered),
  action = "filter",
  detail = "Retained rows with non-missing d15N, mirroring the current QMD measurement selection."
)

log_issue(
  dataset = "museum_batch_cn",
  issue_type = "harmonization_summary",
  severity = "info",
  n_records = nrow(batch_cn_filtered),
  detail = paste0(
    "Batch C/N region harmonization counts: missing->homogenate=",
    sum(is.na(str_extract(batch_cn_filtered$`Identifier 1`, regex("Barb|Shaft", ignore_case = TRUE))), na.rm = TRUE),
    ", Shaft->Rachis=",
    sum(str_detect(str_extract(batch_cn_filtered$`Identifier 1`, regex("Barb|Shaft", ignore_case = TRUE)), regex("shaft", ignore_case = TRUE)), na.rm = TRUE),
    "."
  )
)

batch_cn_long <- batch_cn_filtered %>%
  transmute(
    source_branch = "museum_batch_cn",
    source_file = "data/Master_EA_C_N_SI_Batch 3 vane_barb and Batch 4 Breast Feathers_LEH_Sept25_v1_SB.xls",
    source_sheet = "Sample results_C&N",
    source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`Specimen_ID_`),
    specimen_id_raw = str_extract(`Identifier 1`, "^[^_]+"),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    feather_type_raw = str_extract(`Identifier 1`, regex("Primary|Breast", ignore_case = TRUE)),
    feather_type = standardize_feather_type(feather_type_raw),
    feather_replicate = NA_character_,
    region_raw = str_extract(`Identifier 1`, regex("Barb|Shaft", ignore_case = TRUE)),
    region = standardize_region(region_raw),
    amount_mg = parse_numeric(`Amount (mg)`),
    normalised_d15n = parse_numeric(`normalised d15N`),
    normalised_d13c = parse_numeric(`normalised d13C`)
  ) %>%
  pivot_longer(
    cols = c(normalised_d15n, normalised_d13c),
    names_to = "isotope",
    values_to = "raw_value"
  ) %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf(
      "museum_measurement_%04d",
      nrow(triplicate_cn_long) + nrow(triplicate_h_long) + nrow(triplicate_o_long) + row_number()
    )
  )

log_step(
  dataset = "museum_batch_cn",
  step = "pivot_to_long",
  rows_in = nrow(batch_cn_filtered),
  rows_out = nrow(batch_cn_long),
  action = "reshape",
  detail = "Pivoted batch 3/4 C/N rows to one row per isotope measurement."
)

read_batch_o <- read_triplicate_o

batch_o_block <- read_batch_o %>%
  filter(read_row >= 18, read_row <= 46)

log_step(
  dataset = "museum_batch_o",
  step = "select_batch_block",
  rows_in = nrow(read_batch_o),
  rows_out = nrow(batch_o_block),
  action = "filter",
  detail = "Retained the second O-measurement block (rows 18-46 after read_excel) to mirror the current QMD batch 3/4 section."
)

batch_o_long <- batch_o_block %>%
  filter(!is.na(`normalised d18O`)) %>%
  transmute(
    source_branch = "museum_batch_o",
    source_file = "data/Master_TCEA_O_SI_Breast Feathers_Replication_Batch 4_LEH_Sept25_SB.xls",
    source_sheet = "Sample results_O",
    source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`SIA ID`),
    specimen_id_raw = normalize_space(`SIA ID`),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    feather_type = "Breast",
    feather_replicate = NA_character_,
    region_raw = "homogenate",
    region = "homogenate",
    amount_mg = parse_numeric(Amount),
    isotope = "normalised_d18o",
    raw_value = parse_numeric(`normalised d18O`)
  ) %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf(
      "museum_measurement_%04d",
      nrow(triplicate_cn_long) + nrow(triplicate_h_long) + nrow(triplicate_o_long) + nrow(batch_cn_long) + row_number()
    )
  )

log_step(
  dataset = "museum_batch_o",
  step = "drop_missing_d18o",
  rows_in = nrow(batch_o_block),
  rows_out = nrow(batch_o_long),
  action = "filter",
  detail = "Dropped batch O rows with missing d18O values."
)

read_batch_h <- read_excel(
  path = file.path("data", "Master_TCEA_H_SI_BATCH 4 BREAST Feathers_LE_Nov25_v1.xls"),
  sheet = "Sample results_H",
  col_types = "text",
  skip = 9
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 9)
  )

log_step(
  dataset = "museum_batch_h",
  step = "read_sheet",
  rows_in = nrow(read_batch_h),
  rows_out = nrow(read_batch_h),
  action = "read",
  detail = "Loaded batch 4 H sheet."
)

batch_h_long <- read_batch_h %>%
  filter(!is.na(`Normalised d2H`)) %>%
  transmute(
    source_branch = "museum_batch_h",
    source_file = "data/Master_TCEA_H_SI_BATCH 4 BREAST Feathers_LE_Nov25_v1.xls",
    source_sheet = "Sample results_H",
    source_sheet_row,
    source_record_label = normalize_space(`Identifier 1`),
    source_sia_id = normalize_space(`SIA ID`),
    specimen_id_raw = normalize_space(`SIA ID`),
    specimen_id = standardize_museum_specimen_id(specimen_id_raw),
    feather_type = "Breast",
    feather_replicate = NA_character_,
    region_raw = "homogenate",
    region = "homogenate",
    amount_mg = NA_real_,
    isotope = "normalised_d2h",
    raw_value = parse_numeric(`Normalised d2H`)
  ) %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf(
      "museum_measurement_%04d",
      nrow(triplicate_cn_long) + nrow(triplicate_h_long) + nrow(triplicate_o_long) +
        nrow(batch_cn_long) + nrow(batch_o_long) + row_number()
    )
  )

log_step(
  dataset = "museum_batch_h",
  step = "drop_missing_d2h",
  rows_in = nrow(read_batch_h),
  rows_out = nrow(batch_h_long),
  action = "filter",
  detail = "Dropped batch H rows with missing d2H values."
)

museum_measurements_raw_qc <- bind_rows(
  triplicate_cn_long,
  triplicate_h_long,
  triplicate_o_long,
  batch_cn_long,
  batch_o_long,
  batch_h_long
) %>%
  mutate(
    measurement_harmonization_flags = case_when(
      region == "Rachis" ~ "region_shaft_harmonized_to_rachis",
      region == "homogenate" & source_branch == "museum_batch_cn" & is.na(region_raw) ~ "region_missing_harmonized_to_homogenate",
      TRUE ~ NA_character_
    )
  ) %>%
  apply_assay_qc_long()

log_step(
  dataset = "museum_measurements",
  step = "apply_assay_qc",
  rows_in = nrow(museum_measurements_raw_qc),
  rows_out = sum(museum_measurements_raw_qc$included_in_canonical, na.rm = TRUE),
  action = "adjudicate",
  detail = "Applied assay-level QC decisions to museum isotope measurements while preserving raw values."
)

write_csv_base(
  museum_measurements_raw_qc,
  file.path(derived_dir, "museum_measurements_raw_qc.csv")
)

museum_h_mismatch_usable <- museum_measurements_raw_qc %>%
  filter(
    !is.na(specimen_id_from_sia),
    specimen_id != specimen_id_from_sia,
    included_in_canonical
  )

if (nrow(museum_h_mismatch_usable) > 0) {
  log_issue(
    dataset = "museum_measurements",
    issue_type = "usable_identifier_vs_sia_id_mismatch",
    severity = "error",
    key = collapse_unique(museum_h_mismatch_usable$source_record_label),
    n_records = nrow(museum_h_mismatch_usable),
    detail = "At least one usable museum measurement still has conflicting specimen IDs between Identifier 1 and SIA ID."
  )
  stop("Usable museum measurements still contain Identifier 1 / SIA ID mismatches.")
}

museum_measurements_joinable <- museum_measurements_raw_qc %>%
  filter(included_in_canonical) %>%
  filter(is.na(specimen_id_from_sia) | specimen_id == specimen_id_from_sia) %>%
  select(-specimen_id_from_sia)

museum_measurements_enriched <- museum_measurements_joinable %>%
  left_join(
    museum_metadata_canonical,
    by = c("specimen_id", "feather_type")
  )

metadata_unmatched <- museum_measurements_enriched %>%
  filter(is.na(source_metadata_row_ids))

if (nrow(metadata_unmatched) > 0) {
  log_issue(
    dataset = "museum_measurements",
    issue_type = "empty_metadata_join",
    severity = "error",
    key = collapse_unique(metadata_unmatched$specimen_id),
    n_records = nrow(metadata_unmatched),
    detail = "Museum measurements failed to match canonical museum metadata."
  )
  stop("Museum measurements did not fully join to metadata.")
}

log_step(
  dataset = "museum_measurements",
  step = "join_metadata_exact",
  rows_in = nrow(museum_measurements_joinable),
  rows_out = nrow(museum_measurements_enriched),
  action = "join",
  detail = "Joined museum measurements to metadata exactly on standardized specimen_id + feather_type."
)

coord_unmatched <- museum_measurements_enriched %>%
  filter(is.na(lon) | is.na(lat))

if (nrow(coord_unmatched) > 0) {
  log_issue(
    dataset = "museum_measurements",
    issue_type = "missing_coordinates_after_checkpoint_join",
    severity = "error",
    key = collapse_unique(coord_unmatched$specimen_id),
    n_records = nrow(coord_unmatched),
    detail = "Museum measurements are missing coordinates after joining the historical checkpoint lookup."
  )
  stop("Museum measurements are missing coordinates after checkpoint join.")
}

write_csv_base(
  museum_measurements_enriched,
  file.path(derived_dir, "museum_measurements_long.csv")
)

museum_measurement_dup_keys <- museum_measurements_enriched %>%
  count(source_file, source_sheet, source_sheet_row, isotope, sort = TRUE) %>%
  filter(n > 1)

if (nrow(museum_measurement_dup_keys) > 0) {
  log_issue(
    dataset = "museum_measurements",
    issue_type = "duplicate_source_measurement_key",
    severity = "error",
    key = collapse_unique(paste(museum_measurement_dup_keys$source_file, museum_measurement_dup_keys$source_sheet_row, museum_measurement_dup_keys$isotope)),
    n_records = nrow(museum_measurement_dup_keys),
    detail = "At least one museum source measurement key is duplicated after wrangling."
  )
  stop("Duplicate museum source measurement keys detected.")
}

museum_specimen_region_isotopes_long <- museum_measurements_enriched %>%
  group_by(
    specimen_id,
    feather_type,
    region,
    isotope,
    museum,
    species,
    general_collection_location,
    specific_collection_location,
    location_label,
    collection_date_raw,
    collection_month,
    collection_year,
    geo_region,
    geo_region_bin,
    lon,
    lat
  ) %>%
  summarise(
    value_mean = safe_mean(value),
    value_sd = safe_sd(value),
    n_measurements = n(),
    source_measurement_ids = collapse_unique(measurement_id),
    source_measurement_files = collapse_unique(source_file),
    source_measurement_rows = collapse_unique(source_sheet_row),
    source_branches = collapse_unique(source_branch),
    feather_replicates = collapse_unique(feather_replicate),
    source_metadata_row_ids = summarise_character_or_na(source_metadata_row_ids),
    coordinate_row_ids = summarise_character_or_na(coordinate_row_ids),
    coordinate_source = summarise_character_or_na(coordinate_source),
    .groups = "drop"
  )

write_csv_base(
  museum_specimen_region_isotopes_long,
  file.path(derived_dir, "museum_specimen_region_isotopes_long.csv")
)

museum_wide_id_cols <- c(
  "specimen_id", "feather_type", "region", "museum", "species",
  "general_collection_location", "specific_collection_location", "location_label",
  "collection_date_raw", "collection_month", "collection_year",
  "geo_region", "geo_region_bin", "lon", "lat"
)

museum_specimen_region_isotopes_wide_values <- museum_specimen_region_isotopes_long %>%
  pivot_wider(
    id_cols = all_of(museum_wide_id_cols),
    names_from = isotope,
    values_from = c(value_mean, value_sd, n_measurements, source_measurement_ids),
    names_glue = "{.value}_{isotope}"
  )

museum_specimen_region_isotopes_wide_meta <- museum_specimen_region_isotopes_long %>%
  group_by(across(all_of(museum_wide_id_cols))) %>%
  summarise(
    isotopes_available = collapse_unique(isotope),
    isotopes_available_n = n(),
    source_branches_all = collapse_unique(source_branches),
    source_metadata_row_ids = summarise_character_or_na(source_metadata_row_ids),
    coordinate_row_ids = summarise_character_or_na(coordinate_row_ids),
    coordinate_source = summarise_character_or_na(coordinate_source),
    .groups = "drop"
  )

museum_specimen_region_isotopes_wide <- museum_specimen_region_isotopes_wide_meta %>%
  left_join(
    museum_specimen_region_isotopes_wide_values,
    by = museum_wide_id_cols
  ) %>%
  mutate(
    has_d13c = !is.na(value_mean_normalised_d13c),
    has_d15n = !is.na(value_mean_normalised_d15n),
    has_d18o = !is.na(value_mean_normalised_d18o),
    has_d2h = !is.na(value_mean_normalised_d2h),
    has_all_4_isotopes = has_d13c & has_d15n & has_d18o & has_d2h,
    has_core_3_isotopes = has_d13c & has_d15n & has_d18o
  )

write_csv_base(
  museum_specimen_region_isotopes_wide,
  file.path(derived_dir, "museum_specimen_region_isotopes_wide.csv")
)

museum_screening_step_1 <- museum_specimen_region_isotopes_wide %>%
  filter(feather_type == "Breast")

log_step(
  dataset = "museum_screening_ready",
  step = "filter_breast",
  rows_in = nrow(museum_specimen_region_isotopes_wide),
  rows_out = nrow(museum_screening_step_1),
  action = "filter",
  detail = "Restricted museum screening candidates to breast feathers."
)

museum_screening_step_2 <- museum_screening_step_1 %>%
  filter(region == "homogenate")

log_step(
  dataset = "museum_screening_ready",
  step = "filter_homogenate",
  rows_in = nrow(museum_screening_step_1),
  rows_out = nrow(museum_screening_step_2),
  action = "filter",
  detail = "Restricted museum screening candidates to homogenate measurements to avoid silent region averaging across barb / rachis / homogenate."
)

museum_screening_ready <- museum_screening_step_2 %>%
  filter(!is.na(collection_month) & collection_month %in% 5:8) %>%
  mutate(
    screening_subset = "museum_breast_homogenate_may_august"
  )

log_step(
  dataset = "museum_screening_ready",
  step = "filter_winter_months",
  rows_in = nrow(museum_screening_step_2),
  rows_out = nrow(museum_screening_ready),
  action = "filter",
  detail = "Restricted museum screening candidates to May-August collection months."
)

write_csv_base(
  museum_screening_ready,
  file.path(derived_dir, "museum_screening_ready_breast_homogenate_winter.csv")
)

live_phenotypes_raw <- read_excel(
  path = file.path("data", "LIMS_251105_MaxPlanck_Eberhart-Hertel.xlsx"),
  sheet = "Sample Submission Form",
  skip = 3
) %>%
  mutate(
    phenotype_row_id = sprintf("live_pheno_%03d", row_number()),
    source_file = "data/LIMS_251105_MaxPlanck_Eberhart-Hertel.xlsx",
    source_sheet = "Sample Submission Form",
    source_sheet_row = excel_data_row(row_number(), skip = 3),
    sample_identifier_raw = normalize_space(`Sample Identifier`),
    ring = standardize_live_ring(`Sample Identifier`),
    sample_type = normalize_space(`Sample Type`),
    sex = normalize_space(Sex),
    tag_id_raw = normalize_space(as.character(TagID)),
    date_sampled = as.character(as.Date(`Date Sampled`)),
    longitude_capture = parse_numeric(`Longitude Sampled`),
    latitude_capture = parse_numeric(`Latitude Sampled`),
    status = normalize_space(`Migratory Status`),
    migratory_status_comment = normalize_space(`Migratory Status Comment`),
    tag_id_adjusted = case_when(
      tag_id_raw == "266470" & !is.na(date_sampled) & as.Date(date_sampled) < as.Date("2024-12-09") ~ "266470_A",
      tag_id_raw == "266470" & !is.na(date_sampled) & as.Date(date_sampled) >= as.Date("2024-12-09") ~ "266470_B",
      TRUE ~ tag_id_raw
    ),
    status_bin = case_when(
      str_detect(str_to_lower(status), "migrant") ~ "migrant",
      TRUE ~ status
    )
  )

log_step(
  dataset = "live_phenotypes",
  step = "read_sheet",
  rows_in = nrow(live_phenotypes_raw),
  rows_out = nrow(live_phenotypes_raw),
  action = "read",
  detail = "Loaded live bird LIMS sample submission form."
)

live_phenotype_conflicts <- live_phenotypes_raw %>%
  filter(!is.na(ring)) %>%
  group_by(ring) %>%
  summarise(
    n_sex = distinct_non_missing_n(sex),
    n_tag = distinct_non_missing_n(tag_id_adjusted),
    n_status = distinct_non_missing_n(status),
    n_longitude = distinct_non_missing_n(longitude_capture),
    n_latitude = distinct_non_missing_n(latitude_capture),
    .groups = "drop"
  ) %>%
  filter(n_sex > 1 | n_tag > 1 | n_status > 1 | n_longitude > 1 | n_latitude > 1)

if (nrow(live_phenotype_conflicts) > 0) {
  log_issue(
    dataset = "live_phenotypes",
    issue_type = "conflicting_ring_labels",
    severity = "error",
    key = collapse_unique(live_phenotype_conflicts$ring),
    n_records = nrow(live_phenotype_conflicts),
    detail = "Live phenotype data has conflicting values for the same ring."
  )
  stop("Conflicting ring-level phenotype labels detected.")
}

live_phenotypes_canonical <- live_phenotypes_raw %>%
  filter(!is.na(ring)) %>%
  group_by(ring) %>%
  summarise(
    source_phenotype_row_ids = collapse_unique(phenotype_row_id),
    source_phenotype_row_count = n(),
    sample_identifier_examples = collapse_unique(sample_identifier_raw),
    sex = summarise_character_or_na(sex),
    tag_id_adjusted = summarise_character_or_na(tag_id_adjusted),
    tag_id_raw_examples = collapse_unique(tag_id_raw),
    status = summarise_character_or_na(status),
    status_bin = summarise_character_or_na(status_bin),
    migratory_status_comment = summarise_character_or_na(migratory_status_comment),
    longitude_capture = safe_mean(longitude_capture),
    latitude_capture = safe_mean(latitude_capture),
    .groups = "drop"
  ) %>%
  mutate(
    status_known = !is.na(status) & str_to_lower(status) != "unknown"
  )

log_step(
  dataset = "live_phenotypes",
  step = "collapse_to_ring_level",
  rows_in = nrow(live_phenotypes_raw),
  rows_out = nrow(live_phenotypes_canonical),
  action = "deduplicate",
  detail = "Collapsed LIMS phenotype rows to one canonical row per ring."
)

write_csv_base(
  live_phenotypes_canonical,
  file.path(derived_dir, "live_phenotypes_canonical.csv")
)

live_cn_raw <- read_excel(
  path = file.path("data", "Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls"),
  sheet = "Sample results_C&N",
  col_types = "text",
  skip = 10
) %>%
  mutate(
    read_row = row_number(),
    source_sheet_row = excel_data_row(read_row, skip = 10)
  )

log_step(
  dataset = "live_cn_measurements",
  step = "read_sheet",
  rows_in = nrow(live_cn_raw),
  rows_out = nrow(live_cn_raw),
  action = "read",
  detail = "Loaded live bird batch 5 C/N sheet."
)

live_cn_nonblank <- live_cn_raw %>%
  filter(!is.na(`Corrected Identifer 1`))

log_step(
  dataset = "live_cn_measurements",
  step = "drop_blank_rows",
  rows_in = nrow(live_cn_raw),
  rows_out = nrow(live_cn_nonblank),
  action = "filter",
  detail = "Dropped blank trailing rows from live bird C/N sheet."
)

live_pimary_count <- sum(
  str_detect(
    str_to_lower(live_cn_nonblank$`Corrected Identifer 1`),
    "pimary"
  ),
  na.rm = TRUE
)

log_issue(
  dataset = "live_cn_measurements",
  issue_type = "harmonization_summary",
  severity = "info",
  n_records = live_pimary_count,
  detail = paste0(
    "Corrected feather label typo 'pimary' -> 'Primary' for ",
    live_pimary_count,
    " live bird assay rows."
  )
)

live_cn_measurements_long_raw <- live_cn_nonblank %>%
  transmute(
    assay_row_id = sprintf("live_assay_%04d", row_number()),
    source_branch = "live_batch5_cn",
    source_file = "data/Master_EA_C_N_SI_Batch 5 Breast and Primary Feathers_LEH_Repeats added_April25_v1_SB_v2.xls",
    source_sheet = "Sample results_C&N",
    source_sheet_row,
    sample_id_raw = normalize_space(`Corrected Identifer 1`),
    sample_id_base = strip_live_repeat_suffix(`Corrected Identifer 1`),
    ring = standardize_live_ring(`Corrected Identifer 1`),
    feather_type_raw = str_extract(`Corrected Identifer 1`, regex("primary|breast|pimary", ignore_case = TRUE)),
    feather_type = standardize_feather_type(feather_type_raw),
    is_repeat_assay = str_detect(normalize_space(`Corrected Identifer 1`), "_rpt$") |
      (!is.na(`Rpt?`) & normalize_space(`Rpt?`) != ""),
    amount_mg = parse_numeric(`Amount (mg)`),
    cn_mass_ratio = parse_numeric(`C:N Mass Ratio`),
    normalised_d15n = parse_numeric(`normalised d15N`),
    normalised_d13c = parse_numeric(`normalised d13C`)
  ) %>%
  pivot_longer(
    cols = c(normalised_d15n, normalised_d13c),
    names_to = "isotope",
    values_to = "raw_value"
  ) %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    measurement_id = sprintf("live_measurement_%04d", row_number())
  ) %>%
  left_join(live_phenotypes_canonical, by = "ring") %>%
  apply_assay_qc_long() %>%
  left_join(
    manual_review_lookup,
    by = c(
      "sample_id_base" = "sample_base_identifier",
      "feather_type" = "tissue_feather_type",
      "isotope" = "isotope"
    )
  ) %>%
  mutate(
    assay_qc_manual_review_pending = coalesce(assay_qc_manual_review_pending, FALSE),
    assay_qc_manual_review_reason = ifelse(
      assay_qc_manual_review_pending,
      assay_qc_manual_review_reason,
      NA_character_
    ),
    assay_qc_manual_review_notes = ifelse(
      assay_qc_manual_review_pending,
      assay_qc_manual_review_notes,
      NA_character_
    )
  )

live_cn_unmatched <- live_cn_measurements_long_raw %>%
  filter(is.na(source_phenotype_row_ids))

if (nrow(live_cn_unmatched) > 0) {
  log_issue(
    dataset = "live_cn_measurements",
    issue_type = "empty_phenotype_join",
    severity = "error",
    key = collapse_unique(live_cn_unmatched$ring),
    n_records = nrow(live_cn_unmatched),
    detail = "Live bird isotope measurements failed to match ring-level phenotypes."
  )
  stop("Live bird measurements did not fully join to phenotypes.")
}

log_step(
  dataset = "live_cn_measurements",
  step = "join_ring_level_phenotypes",
  rows_in = nrow(live_cn_measurements_long_raw),
  rows_out = nrow(live_cn_measurements_long_raw),
  action = "join",
  detail = "Joined live bird assays to ring-level phenotypes exactly on ring."
)

log_step(
  dataset = "live_cn_measurements",
  step = "apply_assay_qc",
  rows_in = nrow(live_cn_measurements_long_raw),
  rows_out = sum(live_cn_measurements_long_raw$included_in_canonical, na.rm = TRUE),
  action = "adjudicate",
  detail = "Applied assay-level QC adjudications to live bird C/N measurements."
)

log_step(
  dataset = "live_cn_measurements_strict",
  step = "apply_assay_qc",
  rows_in = nrow(live_cn_measurements_long_raw),
  rows_out = sum(live_cn_measurements_long_raw$included_in_strict, na.rm = TRUE),
  action = "adjudicate",
  detail = "Applied the strict assay-QC value view to live bird C/N measurements."
)

log_step(
  dataset = "live_cn_measurements_sensitivity",
  step = "apply_assay_qc",
  rows_in = nrow(live_cn_measurements_long_raw),
  rows_out = sum(live_cn_measurements_long_raw$included_in_sensitivity, na.rm = TRUE),
  action = "adjudicate",
  detail = "Applied the sensitivity assay-QC value view to live bird C/N measurements."
)

write_csv_base(
  live_cn_measurements_long_raw,
  file.path(derived_dir, "live_cn_measurements_raw_qc_long.csv")
)

live_cn_measurements <- live_cn_measurements_long_raw %>%
  select(
    assay_row_id,
    source_branch,
    source_file,
    source_sheet,
    source_sheet_row,
    sample_id_raw,
    sample_id_base,
    ring,
    feather_type_raw,
    feather_type,
    is_repeat_assay,
    amount_mg,
    cn_mass_ratio,
    sex,
    tag_id_adjusted,
    tag_id_raw_examples,
    status,
    status_bin,
    migratory_status_comment,
    longitude_capture,
    latitude_capture,
    status_known,
    source_phenotype_row_ids,
    source_phenotype_row_count,
    sample_identifier_examples,
    isotope,
    measurement_id,
    raw_value,
    audit_value,
    value,
    strict_value,
    sensitivity_value,
    included_in_canonical,
    included_in_strict,
    included_in_sensitivity,
    assay_qc_action,
    assay_qc_reason,
    assay_qc_evidence_source,
    assay_qc_decision_id,
    assay_qc_explicit,
    assay_qc_value_changed,
    assay_qc_value_excluded,
    assay_qc_repeat_pair_detected,
    assay_qc_repeat_detection_sources,
    assay_qc_repeat_detection_pattern,
    assay_qc_repeat_resolution_status,
    assay_qc_audit_retained,
    assay_qc_excluded_by_lab_instruction,
    assay_qc_contributes_to_model_value,
    assay_qc_model_value_method,
    assay_qc_model_source_rows,
    assay_qc_strict_row_action,
    assay_qc_sensitivity_row_action,
    assay_qc_sensitivity_selected_source_row,
    assay_qc_sensitivity_selected_value,
    assay_qc_manual_review_pending,
    assay_qc_manual_review_reason,
    assay_qc_manual_review_notes
  ) %>%
  pivot_wider(
    id_cols = c(
      assay_row_id,
      source_branch,
      source_file,
      source_sheet,
      source_sheet_row,
      sample_id_raw,
      sample_id_base,
      ring,
      feather_type_raw,
      feather_type,
      is_repeat_assay,
      amount_mg,
      cn_mass_ratio,
      sex,
      tag_id_adjusted,
      tag_id_raw_examples,
      status,
      status_bin,
      migratory_status_comment,
      longitude_capture,
      latitude_capture,
      status_known,
      source_phenotype_row_ids,
      source_phenotype_row_count,
      sample_identifier_examples
    ),
    names_from = isotope,
    values_from = c(
      measurement_id,
      raw_value,
      audit_value,
      value,
      strict_value,
      sensitivity_value,
      included_in_canonical,
      included_in_strict,
      included_in_sensitivity,
      assay_qc_action,
      assay_qc_reason,
      assay_qc_evidence_source,
      assay_qc_decision_id,
      assay_qc_explicit,
      assay_qc_value_changed,
      assay_qc_value_excluded,
      assay_qc_repeat_pair_detected,
      assay_qc_repeat_detection_sources,
      assay_qc_repeat_detection_pattern,
      assay_qc_repeat_resolution_status,
      assay_qc_audit_retained,
      assay_qc_excluded_by_lab_instruction,
      assay_qc_contributes_to_model_value,
      assay_qc_model_value_method,
      assay_qc_model_source_rows,
      assay_qc_strict_row_action,
      assay_qc_sensitivity_row_action,
      assay_qc_sensitivity_selected_source_row,
      assay_qc_sensitivity_selected_value,
      assay_qc_manual_review_pending,
      assay_qc_manual_review_reason,
      assay_qc_manual_review_notes
    ),
    names_glue = "{.value}_{isotope}"
  ) %>%
  rename(
    raw_normalised_d15n = raw_value_normalised_d15n,
    raw_normalised_d13c = raw_value_normalised_d13c,
    audit_normalised_d15n = audit_value_normalised_d15n,
    audit_normalised_d13c = audit_value_normalised_d13c,
    adjudicated_normalised_d15n = value_normalised_d15n,
    adjudicated_normalised_d13c = value_normalised_d13c,
    strict_normalised_d15n = strict_value_normalised_d15n,
    strict_normalised_d13c = strict_value_normalised_d13c,
    sensitivity_normalised_d15n = sensitivity_value_normalised_d15n,
    sensitivity_normalised_d13c = sensitivity_value_normalised_d13c
  )

write_csv_base(
  live_cn_measurements,
  file.path(derived_dir, "live_cn_measurements_by_sample.csv")
)

build_live_cn_variant <- function(measurements_long_raw,
                                  variant_name,
                                  include_col,
                                  value_col,
                                  row_action_col,
                                  file_suffix = "") {
  variant_long <- measurements_long_raw
  variant_long$variant_name <- variant_name
  variant_long$variant_included <- variant_long[[include_col]]
  variant_long$variant_value <- variant_long[[value_col]]
  variant_long$variant_row_action <- variant_long[[row_action_col]]

  variant_dataset_name <- if (file_suffix == "") {
    "live_cn_measurements"
  } else {
    paste0("live_cn_measurements", file_suffix)
  }
  screening_dataset_name <- if (file_suffix == "") {
    "live_screening_ready"
  } else {
    paste0("live_screening_ready", file_suffix)
  }
  variant_retained <- variant_long %>%
    filter(variant_included)

  live_cn_measurements_adjudicated_long <- variant_retained %>%
    group_by(
      source_branch,
      source_file,
      source_sheet,
      sample_id_base,
      ring,
      feather_type_raw,
      feather_type,
      isotope,
      sex,
      tag_id_adjusted,
      tag_id_raw_examples,
      status,
      status_bin,
      migratory_status_comment,
      longitude_capture,
      latitude_capture,
      status_known,
      source_phenotype_row_ids,
      source_phenotype_row_count,
      sample_identifier_examples
    ) %>%
    summarise(
      source_sheet_rows = collapse_unique(source_sheet_row),
      source_sheet_rows_model = collapse_unique(source_sheet_row[assay_qc_contributes_to_model_value %in% TRUE]),
      sample_id_raw = collapse_unique(sample_id_raw),
      is_repeat_assay = any(is_repeat_assay),
      n_retained_assay_rows = n(),
      n_model_contributing_rows = sum(assay_qc_contributes_to_model_value %in% TRUE, na.rm = TRUE),
      raw_value_retained_mean = safe_mean(raw_value),
      audit_value_mean = safe_mean(audit_value),
      value = safe_mean(variant_value[assay_qc_contributes_to_model_value %in% TRUE]),
      amount_mg = safe_mean(amount_mg),
      cn_mass_ratio = safe_mean(cn_mass_ratio),
      source_measurement_ids = collapse_unique(measurement_id),
      source_measurement_ids_model = collapse_unique(measurement_id[assay_qc_contributes_to_model_value %in% TRUE]),
      adjudication_method = collapse_unique(assay_qc_model_value_method),
      assay_qc_variant_name = summarise_character_or_na(variant_name),
      assay_qc_variant_row_action = collapse_unique(variant_row_action),
      assay_qc_action = collapse_unique(assay_qc_action),
      assay_qc_reason = collapse_unique(assay_qc_reason),
      assay_qc_evidence_source = collapse_unique(assay_qc_evidence_source),
      assay_qc_decision_ids = collapse_unique(assay_qc_decision_id),
      assay_qc_repeat_pair_detected = any(assay_qc_repeat_pair_detected %in% TRUE),
      assay_qc_repeat_detection_sources = collapse_unique(assay_qc_repeat_detection_sources),
      assay_qc_repeat_detection_pattern = collapse_unique(assay_qc_repeat_detection_pattern),
      assay_qc_repeat_resolution_status = collapse_unique(assay_qc_repeat_resolution_status),
      assay_qc_contributes_to_model_value = any(assay_qc_contributes_to_model_value %in% TRUE),
      assay_qc_model_source_rows = collapse_unique(assay_qc_model_source_rows),
      assay_qc_manual_review_pending = any(assay_qc_manual_review_pending %in% TRUE),
      assay_qc_manual_review_reason = collapse_unique(assay_qc_manual_review_reason),
      assay_qc_manual_review_notes = collapse_unique(assay_qc_manual_review_notes),
      .groups = "drop"
    ) %>%
    mutate(
      adjudicated_sample_id = paste0("live_model_", make.names(sample_id_base)),
      adjudicated_measurement_id = paste(adjudicated_sample_id, isotope, sep = "::")
    ) %>%
    select(
      adjudicated_sample_id,
      adjudicated_measurement_id,
      source_branch,
      source_file,
      source_sheet,
      source_sheet_rows,
      source_sheet_rows_model,
      sample_id_raw,
      sample_id_base,
      ring,
      feather_type_raw,
      feather_type,
      is_repeat_assay,
      isotope,
      raw_value_retained_mean,
      audit_value_mean,
      value,
      amount_mg,
      cn_mass_ratio,
      sex,
      tag_id_adjusted,
      tag_id_raw_examples,
      status,
      status_bin,
      migratory_status_comment,
      longitude_capture,
      latitude_capture,
      status_known,
      source_phenotype_row_ids,
      source_phenotype_row_count,
      sample_identifier_examples,
      source_measurement_ids,
      source_measurement_ids_model,
      n_retained_assay_rows,
      n_model_contributing_rows,
      adjudication_method,
      assay_qc_variant_name,
      assay_qc_variant_row_action,
      assay_qc_action,
      assay_qc_reason,
      assay_qc_evidence_source,
      assay_qc_decision_ids,
      assay_qc_repeat_pair_detected,
      assay_qc_repeat_detection_sources,
      assay_qc_repeat_detection_pattern,
      assay_qc_repeat_resolution_status,
      assay_qc_contributes_to_model_value,
      assay_qc_model_source_rows,
      assay_qc_manual_review_pending,
      assay_qc_manual_review_reason,
      assay_qc_manual_review_notes
    )

  log_step(
    dataset = variant_dataset_name,
    step = "collapse_retained_assays_to_model_value",
    rows_in = nrow(variant_retained),
    rows_out = nrow(live_cn_measurements_adjudicated_long),
    action = "collapse",
    detail = paste0(
      "Collapsed retained assay rows to one adjudicated sample-isotope modelling value for the ",
      variant_name,
      " value view."
    )
  )

  write_csv_base(
    live_cn_measurements_adjudicated_long,
    file.path(derived_dir, paste0("live_cn_measurements_adjudicated_long", file_suffix, ".csv"))
  )

  live_cn_measurements_adjudicated_by_sample <- live_cn_measurements_adjudicated_long %>%
    pivot_wider(
      id_cols = c(
        adjudicated_sample_id,
        source_branch,
        source_file,
        source_sheet,
        source_sheet_rows,
        source_sheet_rows_model,
        sample_id_raw,
        sample_id_base,
        ring,
        feather_type_raw,
        feather_type,
        is_repeat_assay,
        amount_mg,
        cn_mass_ratio,
        sex,
        tag_id_adjusted,
        tag_id_raw_examples,
        status,
        status_bin,
        migratory_status_comment,
        longitude_capture,
        latitude_capture,
        status_known,
        source_phenotype_row_ids,
        source_phenotype_row_count,
        sample_identifier_examples,
        n_retained_assay_rows,
        n_model_contributing_rows,
        adjudication_method,
        assay_qc_variant_name,
        assay_qc_manual_review_pending,
        assay_qc_manual_review_reason,
        assay_qc_manual_review_notes
      ),
      names_from = isotope,
      values_from = c(
        raw_value_retained_mean,
        audit_value_mean,
        value,
        source_measurement_ids,
        source_measurement_ids_model,
        assay_qc_variant_row_action,
        assay_qc_action,
        assay_qc_reason,
        assay_qc_evidence_source,
        assay_qc_decision_ids,
        assay_qc_repeat_pair_detected,
        assay_qc_repeat_detection_sources,
        assay_qc_repeat_detection_pattern,
        assay_qc_repeat_resolution_status,
        assay_qc_contributes_to_model_value,
        assay_qc_model_source_rows
      ),
      names_glue = "{.value}_{isotope}"
    ) %>%
    rename(
      raw_normalised_d15n = raw_value_retained_mean_normalised_d15n,
      raw_normalised_d13c = raw_value_retained_mean_normalised_d13c,
      audit_normalised_d15n = audit_value_mean_normalised_d15n,
      audit_normalised_d13c = audit_value_mean_normalised_d13c,
      adjudicated_normalised_d15n = value_normalised_d15n,
      adjudicated_normalised_d13c = value_normalised_d13c
    ) %>%
    mutate(
      source_measurement_ids = mapply(
        function(d15n_ids, d13c_ids) collapse_unique(c(d15n_ids, d13c_ids)),
        source_measurement_ids_normalised_d15n,
        source_measurement_ids_normalised_d13c,
        USE.NAMES = FALSE
      )
    )

  write_csv_base(
    live_cn_measurements_adjudicated_by_sample,
    file.path(derived_dir, paste0("live_cn_measurements_adjudicated_by_sample", file_suffix, ".csv"))
  )

  live_cn_tissue_summary <- live_cn_measurements_adjudicated_by_sample %>%
    group_by(
      ring,
      feather_type,
      sex,
      status,
      status_bin,
      status_known,
      longitude_capture,
      latitude_capture,
      tag_id_adjusted,
      migratory_status_comment
    ) %>%
    summarise(
      n_assays = n(),
      n_unique_samples = dplyr::n_distinct(sample_id_base),
      n_retained_assay_rows = sum(n_retained_assay_rows, na.rm = TRUE),
      n_model_contributing_rows = sum(n_model_contributing_rows, na.rm = TRUE),
      sample_ids_raw = collapse_unique(sample_id_raw),
      sample_ids_base = collapse_unique(sample_id_base),
      normalised_d13c = safe_mean(adjudicated_normalised_d13c),
      normalised_d15n = safe_mean(adjudicated_normalised_d15n),
      cn_mass_ratio = safe_mean(cn_mass_ratio),
      amount_mg_mean = safe_mean(amount_mg),
      amount_mg_sd = safe_sd(amount_mg),
      source_measurement_ids = collapse_unique(source_measurement_ids),
      source_phenotype_row_ids = summarise_character_or_na(source_phenotype_row_ids),
      adjudication_methods = collapse_unique(adjudication_method),
      has_pending_manual_qc = any(assay_qc_manual_review_pending %in% TRUE),
      pending_manual_qc_reasons = collapse_unique(assay_qc_manual_review_reason),
      .groups = "drop"
    )

  write_csv_base(
    live_cn_tissue_summary,
    file.path(derived_dir, paste0("live_cn_tissue_summary", file_suffix, ".csv"))
  )

  live_cn_pair_base <- live_cn_tissue_summary %>%
    group_by(
      ring,
      sex,
      status,
      status_bin,
      status_known,
      longitude_capture,
      latitude_capture,
      tag_id_adjusted,
      migratory_status_comment,
      source_phenotype_row_ids
    ) %>%
    summarise(
      n_tissues_present = dplyr::n_distinct(feather_type),
      tissues_present = collapse_unique(feather_type),
      has_any_pending_manual_qc = any(has_pending_manual_qc %in% TRUE),
      pending_manual_qc_reasons_any = collapse_unique(pending_manual_qc_reasons),
      .groups = "drop"
    )

  live_cn_pair_values <- live_cn_tissue_summary %>%
    select(
      ring,
      feather_type,
      normalised_d13c,
      normalised_d15n,
      cn_mass_ratio,
      amount_mg_mean,
      n_assays,
      n_unique_samples,
      source_measurement_ids,
      has_pending_manual_qc,
      pending_manual_qc_reasons
    ) %>%
    pivot_wider(
      id_cols = ring,
      names_from = feather_type,
      values_from = c(
        normalised_d13c,
        normalised_d15n,
        cn_mass_ratio,
        amount_mg_mean,
        n_assays,
        n_unique_samples,
        source_measurement_ids,
        has_pending_manual_qc,
        pending_manual_qc_reasons
      ),
      names_glue = "{.value}_{feather_type}"
    )

  live_cn_paired_by_ring <- live_cn_pair_base %>%
    left_join(live_cn_pair_values, by = "ring") %>%
    mutate(
      has_breast = !is.na(normalised_d13c_Breast) | !is.na(normalised_d15n_Breast),
      has_primary = !is.na(normalised_d13c_Primary) | !is.na(normalised_d15n_Primary),
      has_complete_cn_pair = !is.na(normalised_d13c_Breast) &
        !is.na(normalised_d15n_Breast) &
        !is.na(normalised_d13c_Primary) &
        !is.na(normalised_d15n_Primary),
      has_pending_manual_qc = has_any_pending_manual_qc |
        (has_pending_manual_qc_Breast %in% TRUE) |
        (has_pending_manual_qc_Primary %in% TRUE),
      pending_manual_qc_reasons = mapply(
        function(any_reason, breast_reason, primary_reason) {
          collapse_unique(c(any_reason, breast_reason, primary_reason))
        },
        pending_manual_qc_reasons_any,
        pending_manual_qc_reasons_Breast,
        pending_manual_qc_reasons_Primary,
        USE.NAMES = FALSE
      ),
      delta_d13c_breast_minus_primary = ifelse(
        has_complete_cn_pair,
        normalised_d13c_Breast - normalised_d13c_Primary,
        NA_real_
      ),
      delta_d15n_breast_minus_primary = ifelse(
        has_complete_cn_pair,
        normalised_d15n_Breast - normalised_d15n_Primary,
        NA_real_
      )
    )

  write_csv_base(
    live_cn_paired_by_ring,
    file.path(derived_dir, paste0("live_cn_paired_by_ring", file_suffix, ".csv"))
  )

  live_screening_step_1 <- live_cn_paired_by_ring %>%
    filter(status_known)

  log_step(
    dataset = screening_dataset_name,
    step = "filter_known_status",
    rows_in = nrow(live_cn_paired_by_ring),
    rows_out = nrow(live_screening_step_1),
    action = "filter",
    detail = paste0(
      "Restricted live bird screening candidates to rings with known migratory status for the ",
      variant_name,
      " value view."
    )
  )

  live_screening_ready <- live_screening_step_1 %>%
    filter(has_complete_cn_pair)

  log_step(
    dataset = screening_dataset_name,
    step = "filter_complete_cn_pair",
    rows_in = nrow(live_screening_step_1),
    rows_out = nrow(live_screening_ready),
    action = "filter",
    detail = paste0(
      "Restricted live bird screening candidates to complete primary + breast C/N pairs for the ",
      variant_name,
      " value view."
    )
  )

  write_csv_base(
    live_screening_ready,
    file.path(derived_dir, paste0("live_screening_ready_paired_labelled", file_suffix, ".csv"))
  )

  list(
    adjudicated_long = live_cn_measurements_adjudicated_long,
    adjudicated_by_sample = live_cn_measurements_adjudicated_by_sample,
    tissue_summary = live_cn_tissue_summary,
    paired_by_ring = live_cn_paired_by_ring,
    screening_ready = live_screening_ready
  )
}

live_cn_variant_inclusive <- build_live_cn_variant(
  measurements_long_raw = live_cn_measurements_long_raw,
  variant_name = "inclusive",
  include_col = "included_in_canonical",
  value_col = "value",
  row_action_col = "assay_qc_action",
  file_suffix = ""
)

live_cn_variant_strict <- build_live_cn_variant(
  measurements_long_raw = live_cn_measurements_long_raw,
  variant_name = "strict",
  include_col = "included_in_strict",
  value_col = "strict_value",
  row_action_col = "assay_qc_strict_row_action",
  file_suffix = "_strict"
)

live_cn_variant_sensitivity <- build_live_cn_variant(
  measurements_long_raw = live_cn_measurements_long_raw,
  variant_name = "sensitivity",
  include_col = "included_in_sensitivity",
  value_col = "sensitivity_value",
  row_action_col = "assay_qc_sensitivity_row_action",
  file_suffix = "_sensitivity"
)

live_cn_measurements_adjudicated_long <- live_cn_variant_inclusive$adjudicated_long
live_cn_measurements_adjudicated_by_sample <- live_cn_variant_inclusive$adjudicated_by_sample
live_cn_tissue_summary <- live_cn_variant_inclusive$tissue_summary
live_cn_paired_by_ring <- live_cn_variant_inclusive$paired_by_ring
live_screening_ready <- live_cn_variant_inclusive$screening_ready

live_cn_measurements_adjudicated_long_strict <- live_cn_variant_strict$adjudicated_long
live_cn_measurements_adjudicated_by_sample_strict <- live_cn_variant_strict$adjudicated_by_sample
live_cn_tissue_summary_strict <- live_cn_variant_strict$tissue_summary
live_cn_paired_by_ring_strict <- live_cn_variant_strict$paired_by_ring
live_screening_ready_strict <- live_cn_variant_strict$screening_ready

live_cn_measurements_adjudicated_long_sensitivity <- live_cn_variant_sensitivity$adjudicated_long
live_cn_measurements_adjudicated_by_sample_sensitivity <- live_cn_variant_sensitivity$adjudicated_by_sample
live_cn_tissue_summary_sensitivity <- live_cn_variant_sensitivity$tissue_summary
live_cn_paired_by_ring_sensitivity <- live_cn_variant_sensitivity$paired_by_ring
live_screening_ready_sensitivity <- live_cn_variant_sensitivity$screening_ready

assay_qc_summary <- bind_rows(
  museum_measurements_raw_qc %>%
    transmute(
      dataset = "museum_measurements_raw_qc",
      isotope,
      assay_qc_action,
      assay_qc_explicit,
      assay_qc_value_changed,
      assay_qc_value_excluded,
      included_in_canonical,
      included_in_strict,
      included_in_sensitivity,
      assay_qc_repeat_pair_detected,
      assay_qc_repeat_resolution_status,
      assay_qc_audit_retained,
      assay_qc_excluded_by_lab_instruction,
      assay_qc_contributes_to_model_value,
      assay_qc_model_value_method,
      assay_qc_strict_row_action,
      assay_qc_sensitivity_row_action,
      assay_qc_manual_review_pending = FALSE
    ),
  live_cn_measurements_long_raw %>%
    transmute(
      dataset = "live_cn_measurements_raw_qc_long",
      isotope,
      assay_qc_action,
      assay_qc_explicit,
      assay_qc_value_changed,
      assay_qc_value_excluded,
      included_in_canonical,
      included_in_strict,
      included_in_sensitivity,
      assay_qc_repeat_pair_detected,
      assay_qc_repeat_resolution_status,
      assay_qc_audit_retained,
      assay_qc_excluded_by_lab_instruction,
      assay_qc_contributes_to_model_value,
      assay_qc_model_value_method,
      assay_qc_strict_row_action,
      assay_qc_sensitivity_row_action,
      assay_qc_manual_review_pending
    )
) %>%
  group_by(dataset, isotope) %>%
  summarise(
    n_rows = n(),
    n_included_rows = sum(included_in_canonical, na.rm = TRUE),
    n_included_strict_rows = sum(included_in_strict, na.rm = TRUE),
    n_included_sensitivity_rows = sum(included_in_sensitivity, na.rm = TRUE),
    n_audit_retained_rows = sum(assay_qc_audit_retained %in% TRUE, na.rm = TRUE),
    n_explicit_lab_exclusion_rows = sum(
      assay_qc_excluded_by_lab_instruction %in% TRUE,
      na.rm = TRUE
    ),
    n_model_contributing_rows = sum(
      assay_qc_contributes_to_model_value %in% TRUE,
      na.rm = TRUE
    ),
    n_explicit_qc_rows = sum(assay_qc_explicit, na.rm = TRUE),
    n_excluded_rows = sum(assay_qc_value_excluded, na.rm = TRUE),
    n_lab_green_average_rows = sum(
      assay_qc_model_value_method == "lab_green_average",
      na.rm = TRUE
    ),
    n_explicit_selected_rows = sum(
      assay_qc_model_value_method == "explicit_use_selected_row" &
        assay_qc_contributes_to_model_value %in% TRUE,
      na.rm = TRUE
    ),
    n_mean_of_valid_repeat_rows = sum(
      assay_qc_model_value_method == "mean_of_valid_repeats" &
        assay_qc_contributes_to_model_value %in% TRUE,
      na.rm = TRUE
    ),
    n_changed_nonexcluded_rows = sum(assay_qc_value_changed, na.rm = TRUE),
    n_repeat_pair_rows = sum(assay_qc_repeat_pair_detected %in% TRUE, na.rm = TRUE),
    n_repeat_rows_retained = sum(
      assay_qc_repeat_pair_detected %in% TRUE &
        assay_qc_audit_retained %in% TRUE,
      na.rm = TRUE
    ),
    n_conflicting_repeat_rows = sum(
      assay_qc_repeat_resolution_status == "conflicting_repeat_instructions",
      na.rm = TRUE
    ),
    n_strict_excluded_rows = sum(
      assay_qc_strict_row_action %in% c("exclude", "exclude_unresolved_repeat_pair"),
      na.rm = TRUE
    ),
    n_sensitivity_excluded_rows = sum(
      assay_qc_sensitivity_row_action %in% c(
        "exclude",
        "exclude_nonselected_repeat_row",
        "manual_review_only"
      ),
      na.rm = TRUE
    ),
    n_sensitivity_default_rows = sum(
      assay_qc_sensitivity_row_action == "include_non_rpt_default",
      na.rm = TRUE
    ),
    n_manual_review_pending_rows = sum(assay_qc_manual_review_pending %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_base(
  assay_qc_summary,
  file.path(derived_dir, "assay_qc_summary.csv")
)

missingness_summary <- bind_rows(
  build_missingness_summary(
    "museum_measurements_long",
    museum_measurements_enriched,
    c("value", "collection_month", "collection_year", "lon", "lat")
  ),
  build_missingness_summary(
    "museum_specimen_region_isotopes_wide",
    museum_specimen_region_isotopes_wide,
    c(
      "value_mean_normalised_d13c",
      "value_mean_normalised_d15n",
      "value_mean_normalised_d18o",
      "value_mean_normalised_d2h",
      "collection_month",
      "lon",
      "lat"
    )
  ),
  build_missingness_summary(
    "museum_screening_ready_breast_homogenate_winter",
    museum_screening_ready,
    c(
      "value_mean_normalised_d13c",
      "value_mean_normalised_d15n",
      "value_mean_normalised_d18o",
      "value_mean_normalised_d2h"
    )
  ),
  build_missingness_summary(
    "live_cn_measurements_by_sample",
    live_cn_measurements,
    c("adjudicated_normalised_d13c", "adjudicated_normalised_d15n", "status", "sex")
  ),
  build_missingness_summary(
    "live_cn_paired_by_ring",
    live_cn_paired_by_ring,
    c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary",
      "status",
      "sex"
    )
  ),
  build_missingness_summary(
    "live_cn_paired_by_ring_strict",
    live_cn_paired_by_ring_strict,
    c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary",
      "status",
      "sex"
    )
  ),
  build_missingness_summary(
    "live_cn_paired_by_ring_sensitivity",
    live_cn_paired_by_ring_sensitivity,
    c(
      "normalised_d13c_Breast",
      "normalised_d13c_Primary",
      "normalised_d15n_Breast",
      "normalised_d15n_Primary",
      "status",
      "sex"
    )
  )
)

write_csv_base(
  missingness_summary,
  file.path(derived_dir, "qc_missingness_summary.csv")
)

completeness_summary <- bind_rows(
  museum_screening_ready %>%
    count(isotopes_available_n, name = "n_rows") %>%
    mutate(
      dataset = "museum_screening_ready_breast_homogenate_winter",
      metric = "isotopes_available_n",
      value = as.character(isotopes_available_n)
    ) %>%
    select(dataset, metric, value, n_rows),
  museum_screening_ready %>%
    count(has_all_4_isotopes, name = "n_rows") %>%
    mutate(
      dataset = "museum_screening_ready_breast_homogenate_winter",
      metric = "has_all_4_isotopes",
      value = as.character(has_all_4_isotopes)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring %>%
    count(n_tissues_present, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring",
      metric = "n_tissues_present",
      value = as.character(n_tissues_present)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring %>%
    count(has_complete_cn_pair, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring",
      metric = "has_complete_cn_pair",
      value = as.character(has_complete_cn_pair)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring_strict %>%
    count(n_tissues_present, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring_strict",
      metric = "n_tissues_present",
      value = as.character(n_tissues_present)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring_strict %>%
    count(has_complete_cn_pair, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring_strict",
      metric = "has_complete_cn_pair",
      value = as.character(has_complete_cn_pair)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring_sensitivity %>%
    count(n_tissues_present, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring_sensitivity",
      metric = "n_tissues_present",
      value = as.character(n_tissues_present)
    ) %>%
    select(dataset, metric, value, n_rows),
  live_cn_paired_by_ring_sensitivity %>%
    count(has_complete_cn_pair, name = "n_rows") %>%
    mutate(
      dataset = "live_cn_paired_by_ring_sensitivity",
      metric = "has_complete_cn_pair",
      value = as.character(has_complete_cn_pair)
    ) %>%
    select(dataset, metric, value, n_rows)
)

write_csv_base(
  completeness_summary,
  file.path(derived_dir, "qc_completeness_summary.csv")
)

write_csv_base(
  qc_steps,
  file.path(derived_dir, "qc_step_log.csv")
)

write_csv_base(
  qc_issues,
  file.path(derived_dir, "qc_issues.csv")
)
