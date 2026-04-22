suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

safe_ref_sd <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  if (is.na(out) || out == 0) {
    return(1)
  }
  out
}

build_primary_site_reference_summary <- function(primary_reference_rows) {
  primary_reference_rows %>%
    filter(
      !is.na(sampling_site_label),
      stats::complete.cases(normalised_d13c, normalised_d15n)
    ) %>%
    group_by(sampling_site_label) %>%
    summarise(
      n_reference_rows = n(),
      mean_primary_d13c = mean(normalised_d13c, na.rm = TRUE),
      mean_primary_d15n = mean(normalised_d15n, na.rm = TRUE),
      sd_primary_d13c = safe_ref_sd(normalised_d13c),
      sd_primary_d15n = safe_ref_sd(normalised_d15n),
      .groups = "drop"
    )
}

add_primary_site_reference_features <- function(data, primary_reference_rows) {
  ref_tbl <- build_primary_site_reference_summary(primary_reference_rows)

  if (nrow(ref_tbl) == 0) {
    return(
      data %>%
        mutate(
          primary_ref_min_site_dist_cn = NA_real_,
          primary_ref_gap_to_second_cn = NA_real_,
          primary_ref_n_sites_available = 0L
        )
    )
  }

  dist_rows <- lapply(
    seq_len(nrow(data)),
    function(i) {
      obs <- data[i, , drop = FALSE]

      if (!stats::complete.cases(obs[, c("normalised_d13c_Primary", "normalised_d15n_Primary"), drop = FALSE])) {
        return(
          tibble(
            primary_ref_min_site_dist_cn = NA_real_,
            primary_ref_gap_to_second_cn = NA_real_,
            primary_ref_n_sites_available = nrow(ref_tbl)
          )
        )
      }

      dist_vec <- sqrt(
        ((obs$normalised_d13c_Primary - ref_tbl$mean_primary_d13c) / ref_tbl$sd_primary_d13c)^2 +
          ((obs$normalised_d15n_Primary - ref_tbl$mean_primary_d15n) / ref_tbl$sd_primary_d15n)^2
      )

      sorted_dist <- sort(dist_vec)
      gap_val <- if (length(sorted_dist) >= 2) sorted_dist[[2]] - sorted_dist[[1]] else NA_real_

      tibble(
        primary_ref_min_site_dist_cn = min(dist_vec, na.rm = TRUE),
        primary_ref_gap_to_second_cn = gap_val,
        primary_ref_n_sites_available = nrow(ref_tbl)
      )
    }
  )

  bind_cols(data, bind_rows(dist_rows))
}

normalize_winter_region <- function(x) {
  dplyr::case_when(
    x %in% c("Australia", "AU") ~ "Australia",
    x %in% c("New Zealand", "NZ") ~ "New Zealand",
    TRUE ~ as.character(x)
  )
}

build_winter_region_reference_summary <- function(museum_reference_rows) {
  museum_reference_rows %>%
    mutate(geo_region = normalize_winter_region(geo_region)) %>%
    filter(
      geo_region %in% c("Australia", "New Zealand"),
      stats::complete.cases(normalised_d13c_Breast, normalised_d15n_Breast, normalised_d2h_Breast)
    ) %>%
    group_by(geo_region) %>%
    summarise(
      n_reference_rows = n(),
      mean_breast_d13c = mean(normalised_d13c_Breast, na.rm = TRUE),
      mean_breast_d15n = mean(normalised_d15n_Breast, na.rm = TRUE),
      mean_breast_d2h = mean(normalised_d2h_Breast, na.rm = TRUE),
      sd_breast_d13c = safe_ref_sd(normalised_d13c_Breast),
      sd_breast_d15n = safe_ref_sd(normalised_d15n_Breast),
      sd_breast_d2h = safe_ref_sd(normalised_d2h_Breast),
      .groups = "drop"
    )
}

add_winter_region_reference_features <- function(data, museum_reference_rows) {
  ref_tbl <- build_winter_region_reference_summary(museum_reference_rows)

  nz_ref <- ref_tbl %>% filter(geo_region == "New Zealand")
  au_ref <- ref_tbl %>% filter(geo_region == "Australia")

  if (nrow(nz_ref) != 1 || nrow(au_ref) != 1) {
    return(
      data %>%
        mutate(
          winter_ref_dist_nz_cnh = NA_real_,
          winter_ref_dist_au_cnh = NA_real_,
          winter_ref_au_minus_nz_cnh = NA_real_
        )
    )
  }

  out_rows <- lapply(
    seq_len(nrow(data)),
    function(i) {
      obs <- data[i, , drop = FALSE]

      if (!stats::complete.cases(obs[, c("normalised_d13c_Breast", "normalised_d15n_Breast", "normalised_d2h_Breast"), drop = FALSE])) {
        return(
          tibble(
            winter_ref_dist_nz_cnh = NA_real_,
            winter_ref_dist_au_cnh = NA_real_,
            winter_ref_au_minus_nz_cnh = NA_real_
          )
        )
      }

      dist_nz <- sqrt(
        ((obs$normalised_d13c_Breast - nz_ref$mean_breast_d13c) / nz_ref$sd_breast_d13c)^2 +
          ((obs$normalised_d15n_Breast - nz_ref$mean_breast_d15n) / nz_ref$sd_breast_d15n)^2 +
          ((obs$normalised_d2h_Breast - nz_ref$mean_breast_d2h) / nz_ref$sd_breast_d2h)^2
      )

      dist_au <- sqrt(
        ((obs$normalised_d13c_Breast - au_ref$mean_breast_d13c) / au_ref$sd_breast_d13c)^2 +
          ((obs$normalised_d15n_Breast - au_ref$mean_breast_d15n) / au_ref$sd_breast_d15n)^2 +
          ((obs$normalised_d2h_Breast - au_ref$mean_breast_d2h) / au_ref$sd_breast_d2h)^2
      )

      tibble(
        winter_ref_dist_nz_cnh = dist_nz,
        winter_ref_dist_au_cnh = dist_au,
        winter_ref_au_minus_nz_cnh = dist_nz - dist_au
      )
    }
  )

  bind_cols(data, bind_rows(out_rows))
}
