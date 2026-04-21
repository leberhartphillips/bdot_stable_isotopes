suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
})

normalize_space <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[[:space:]]+", " ")
  x <- str_trim(x)
  x[x %in% c("", "NA")] <- NA_character_
  x
}

parse_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & as.character(x) != ""]
  if (length(x) == 0) {
    return(NA)
  }
  x[[1]]
}

collapse_unique <- function(x, sep = "|") {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(sort(x), collapse = sep)
}

distinct_non_missing_n <- function(x) {
  x <- as.character(x)
  x <- unique(x[!is.na(x) & x != ""])
  length(x)
}

excel_data_row <- function(read_row, skip) {
  as.integer(read_row) + as.integer(skip) + 1L
}

standardize_feather_type <- function(x) {
  x <- normalize_space(x)
  case_when(
    str_detect(str_to_lower(x), "breast") ~ "Breast",
    str_detect(str_to_lower(x), "primary|pimary") ~ "Primary",
    TRUE ~ NA_character_
  )
}

standardize_region <- function(x) {
  x <- normalize_space(x)
  case_when(
    is.na(x) ~ "homogenate",
    str_detect(str_to_lower(x), "barb") ~ "Barb",
    str_detect(str_to_lower(x), "shaft|rachis") ~ "Rachis",
    TRUE ~ str_to_title(x)
  )
}

standardize_museum_specimen_id <- function(x) {
  x <- normalize_space(x)
  x <- str_remove(x, regex("\\s+(Breast|Primary)(_[A-Z])?$", ignore_case = TRUE))
  x <- str_remove(x, regex("_[A-Z]\\s+feather$", ignore_case = TRUE))
  x <- str_remove(x, regex("\\s+feather$", ignore_case = TRUE))
  x <- str_remove(x, regex("_[A-Z]$", ignore_case = FALSE))
  x <- str_remove_all(x, regex("^(AJB|WAM)[[:space:]]+", ignore_case = TRUE))
  x <- str_remove_all(x, regex("[[:space:]]+-$"))
  normalize_space(x)
}

standardize_live_ring <- function(x) {
  x <- normalize_space(x)
  str_extract(x, "^[^_]+")
}

strip_live_repeat_suffix <- function(x) {
  x <- normalize_space(x)
  str_remove(x, "_rpt$")
}

derive_geo_region <- function(general_collection_location) {
  general_collection_location <- normalize_space(general_collection_location)
  case_when(
    str_detect(general_collection_location, "Australia") ~ "Australia",
    str_detect(general_collection_location, "New Zealand") ~ "New Zealand",
    TRUE ~ NA_character_
  )
}

extract_month_bdot <- function(date_str) {
  if (length(date_str) > 1) {
    return(vapply(date_str, extract_month_bdot, numeric(1)))
  }

  if (is.na(date_str)) {
    return(NA_real_)
  }

  date_str <- str_trim(as.character(date_str))

  if (grepl("^\\d{2}/\\d{2}/\\d{4}$", date_str)) {
    return(as.numeric(month(dmy(date_str))))
  }

  if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", date_str)) {
    return(as.numeric(month(ymd(date_str))))
  }

  if (grepl("^\\d{4,5}$", date_str)) {
    return(as.numeric(month(as.Date(as.numeric(date_str), origin = "1899-12-30"))))
  }

  if (nchar(date_str) == 4) {
    return(NA_real_)
  }

  parsed_text <- suppressWarnings(dmy(gsub(",", "", date_str)))
  if (!is.na(parsed_text)) {
    return(as.numeric(month(parsed_text)))
  }

  if (grepl("[A-Za-z]+[^0-9A-Za-z]*\\d{4}", date_str)) {
    parsed_month_year <- suppressWarnings(dmy(paste0("01 ", gsub(",", "", date_str))))
    if (!is.na(parsed_month_year)) {
      return(as.numeric(month(parsed_month_year)))
    }
  }

  NA_real_
}

extract_year_bdot <- function(date_str) {
  if (length(date_str) > 1) {
    return(vapply(date_str, extract_year_bdot, numeric(1)))
  }

  if (is.na(date_str)) {
    return(NA_real_)
  }

  date_str <- str_trim(as.character(date_str))

  if (grepl("^\\d+$", date_str)) {
    num <- as.numeric(date_str)
    if (num >= 1800 && num <= 2100) {
      return(as.numeric(num))
    }
    if (num >= 0 && num < 100000) {
      dt <- suppressWarnings(as.Date(num, origin = "1899-12-30"))
      if (!is.na(dt)) {
        return(as.numeric(year(dt)))
      }
    }
    return(NA_real_)
  }

  if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", date_str)) {
    parsed <- suppressWarnings(ymd(date_str))
    if (!is.na(parsed)) {
      return(as.numeric(year(parsed)))
    }
  }

  if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", date_str)) {
    parsed <- suppressWarnings(dmy(date_str))
    if (!is.na(parsed)) {
      return(as.numeric(year(parsed)))
    }
  }

  parsed_text <- suppressWarnings(dmy(gsub(",", "", date_str)))
  if (!is.na(parsed_text)) {
    return(as.numeric(year(parsed_text)))
  }

  if (grepl("[A-Za-z]+[^0-9A-Za-z]*\\d{4}", date_str)) {
    parsed_month_year <- suppressWarnings(dmy(paste0("01 ", gsub(",", "", date_str))))
    if (!is.na(parsed_month_year)) {
      return(as.numeric(year(parsed_month_year)))
    }
  }

  if (grepl("^pre-\\d{4}$", date_str, ignore.case = TRUE)) {
    return(as.numeric(sub("(?i)pre-", "", date_str, perl = TRUE)))
  }

  if (grepl("\\b\\d{4}\\b", date_str)) {
    yr <- regmatches(date_str, regexpr("\\b\\d{4}\\b", date_str))
    yr_num <- as.numeric(yr)
    if (!is.na(yr_num) && yr_num >= 1000 && yr_num <= 2100) {
      return(yr_num)
    }
  }

  NA_real_
}

build_missingness_summary <- function(dataset_name, df, vars) {
  tibble(
    dataset = dataset_name,
    variable = vars,
    n_total = nrow(df),
    n_missing = vapply(vars, function(var) sum(is.na(df[[var]])), integer(1)),
    prop_missing = if (nrow(df) == 0) NA_real_ else n_missing / n_total
  )
}
