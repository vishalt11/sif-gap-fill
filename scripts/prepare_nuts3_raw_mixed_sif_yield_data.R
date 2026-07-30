library(tidyverse)
library(sf)
library(lubridate)

# Build the raw, mixed-pixel OCO-2 SIF baselines for the winter-wheat yield
# comparison. The two baselines differ only in footprint selection:
#
#   1. raw_all:  all accepted OCO-2 footprints in a study NUTS3 region
#   2. raw_ww10: accepted footprints with at least 10% winter-wheat coverage
#
# Within each NUTS3 region and month, footprints are averaged by acquisition
# date first. The monthly value is then the unweighted mean of the date means,
# so a date with many adjacent OCO-2 soundings does not dominate the result.

input_sif_rds <- paste0(
  "data/extracted_modis_data/",
  "modis_1_12_bin_uncertainity_corrected_raw.rds"
)
nuts3_regions_rds <- "data/nuts3_regions_80pct_in_mgrs.rds"

crop_pure_wide_csv <- paste0(
  "data/winter_wheat_yield_model/monthly_crop_pure_sif/",
  "nuts3_monthly_strict_wheat_sif_nirv_wide.csv"
)

output_dir <- "data/winter_wheat_yield_model/raw_mixed_sif"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_years <- c(2019L, 2020L, 2021L, 2022L, 2024L)
analysis_months <- 3:7
winter_wheat_threshold <- 0.10
spatial_join_crs <- 3035

required_sif_cols <- c(
  "target_modis_sif",
  "final_check_modis_sif",
  "Delta_Date",
  "Latitude",
  "Longitude",
  "ww_pct"
)

month_lookup <- tibble(
  month = analysis_months,
  month_name = month.name[analysis_months]
)

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

sum_or_zero <- function(x) {
  sum(x, na.rm = TRUE)
}

weighted_mean_or_na <- function(x, weights) {
  valid <- is.finite(x) & is.finite(weights) & weights > 0

  if (!any(valid)) {
    return(NA_real_)
  }

  weighted.mean(x[valid], w = weights[valid])
}

write_table_pair <- function(data, path_without_extension) {
  saveRDS(data, paste0(path_without_extension, ".rds"))
  readr::write_csv(data, paste0(path_without_extension, ".csv"), na = "")
}

monthly_to_wide <- function(monthly_data, prefix, skeleton) {
  key_cols <- c("nuts_id", "nuts3", "mgrs_tile", "year")

  sif_wide <- monthly_data %>%
    select(all_of(key_cols), month_name, raw_sif) %>%
    pivot_wider(
      names_from = month_name,
      values_from = raw_sif,
      names_prefix = paste0(prefix, "_SIF_")
    )

  footprint_count_wide <- monthly_data %>%
    select(all_of(key_cols), month_name, n_raw_footprints) %>%
    pivot_wider(
      names_from = month_name,
      values_from = n_raw_footprints,
      names_prefix = paste0(prefix, "_n_footprints_")
    )

  date_count_wide <- monthly_data %>%
    select(all_of(key_cols), month_name, n_raw_dates) %>%
    pivot_wider(
      names_from = month_name,
      values_from = n_raw_dates,
      names_prefix = paste0(prefix, "_n_dates_")
    )

  wheat_coverage_wide <- monthly_data %>%
    select(all_of(key_cols), month_name, mean_footprint_ww_pct) %>%
    pivot_wider(
      names_from = month_name,
      values_from = mean_footprint_ww_pct,
      names_prefix = paste0(prefix, "_mean_footprint_ww_pct_")
    )

  skeleton %>%
    left_join(sif_wide, by = key_cols) %>%
    left_join(footprint_count_wide, by = key_cols) %>%
    left_join(date_count_wide, by = key_cols) %>%
    left_join(wheat_coverage_wide, by = key_cols)
}

summarize_raw_sif <- function(assigned_sif, baseline_name) {
  daily <- assigned_sif %>%
    group_by(
      nuts_id,
      nuts3,
      mgrs_tile,
      year,
      month,
      month_name,
      Delta_Date
    ) %>%
    summarise(
      daily_raw_sif = mean(target_modis_sif),
      n_raw_footprints = n(),
      mean_footprint_ww_pct = mean_or_na(ww_pct),
      .groups = "drop"
    )

  monthly <- daily %>%
    group_by(
      nuts_id,
      nuts3,
      mgrs_tile,
      year,
      month,
      month_name
    ) %>%
    summarise(
      raw_sif = mean(daily_raw_sif),
      mean_footprint_ww_pct = weighted_mean_or_na(
        mean_footprint_ww_pct,
        weights = n_raw_footprints
      ),
      n_raw_footprints = sum_or_zero(n_raw_footprints),
      n_raw_dates = n_distinct(Delta_Date),
      .groups = "drop"
    ) %>%
    mutate(
      baseline = baseline_name,
      .before = 1
    )

  list(daily = daily, monthly = monthly)
}

message("Reading raw OCO-2 SIF data...")
sif_data <- readRDS(input_sif_rds)

if (inherits(sif_data, "sf")) {
  sif_data <- st_drop_geometry(sif_data)
}

missing_sif_cols <- setdiff(required_sif_cols, names(sif_data))
if (length(missing_sif_cols) > 0) {
  stop(
    "Raw SIF data are missing required columns: ",
    paste(missing_sif_cols, collapse = ", ")
  )
}

message("Reading the 32 study NUTS3 regions...")
nuts3_regions <- readRDS(nuts3_regions_rds)

if (!inherits(nuts3_regions, "sf")) {
  stop("The NUTS3 region RDS must contain an sf object.")
}

required_nuts_cols <- c("nuts_id", "nuts3", "mgrs_tile")
missing_nuts_cols <- setdiff(required_nuts_cols, names(nuts3_regions))
if (length(missing_nuts_cols) > 0) {
  stop(
    "NUTS3 regions are missing required columns: ",
    paste(missing_nuts_cols, collapse = ", ")
  )
}

nuts3_regions <- nuts3_regions %>%
  select(nuts_id, nuts3, mgrs_tile, any_of("nuts3_area_in_mgrs_pct")) %>%
  st_make_valid()

if (is.na(st_crs(nuts3_regions))) {
  stop("The NUTS3 regions have no CRS.")
}

nuts3_lookup <- nuts3_regions %>%
  st_drop_geometry() %>%
  distinct(nuts_id, mgrs_tile, .keep_all = TRUE)

if (nrow(nuts3_lookup) != 32) {
  warning(
    "Expected 32 unique NUTS3/tile rows but found ",
    nrow(nuts3_lookup),
    "."
  )
}

model_skeleton <- crossing(
  nuts3_lookup,
  year = analysis_years
) %>%
  arrange(mgrs_tile, nuts_id, year)

accepted_sif <- sif_data %>%
  mutate(
    .raw_sif_row_id = row_number(),
    Delta_Date = as.Date(Delta_Date),
    target_modis_sif = as.numeric(target_modis_sif),
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude),
    ww_pct = as.numeric(ww_pct),
    final_check_modis_sif = str_to_lower(
      str_trim(as.character(final_check_modis_sif))
    ),
    year = year(Delta_Date),
    month = month(Delta_Date)
  ) %>%
  filter(
    final_check_modis_sif == "accept",
    year %in% analysis_years,
    month %in% analysis_months,
    is.finite(target_modis_sif),
    is.finite(Latitude),
    is.finite(Longitude)
  ) %>%
  left_join(month_lookup, by = "month")

finite_ww_pct <- accepted_sif$ww_pct[is.finite(accepted_sif$ww_pct)]
if (
  length(finite_ww_pct) > 0 &&
  (min(finite_ww_pct) < 0 || max(finite_ww_pct) > 1.001)
) {
  stop(
    "ww_pct is expected to be a fraction from 0 to 1. Observed range: ",
    paste(range(finite_ww_pct), collapse = " to ")
  )
}

message(
  "Accepted SIF rows before NUTS3 assignment: ",
  format(nrow(accepted_sif), big.mark = ",")
)

# Latitude and Longitude are the supplied OCO-2 footprint-centre coordinates.
# Using them avoids expensive polygon construction when only region membership
# of each footprint centroid is required.
sif_centroids <- accepted_sif %>%
  st_as_sf(
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
  st_transform(spatial_join_crs)

nuts3_for_join <- nuts3_regions %>%
  st_transform(spatial_join_crs)

assigned_sif <- st_join(
  sif_centroids,
  nuts3_for_join,
  join = st_intersects,
  left = FALSE
) %>%
  arrange(.raw_sif_row_id, nuts_id) %>%
  distinct(.raw_sif_row_id, .keep_all = TRUE) %>%
  st_drop_geometry()

message(
  "Accepted SIF rows assigned to the 32 NUTS3 regions: ",
  format(nrow(assigned_sif), big.mark = ",")
)

assigned_export <- assigned_sif %>%
  select(
    .raw_sif_row_id,
    nuts_id,
    nuts3,
    mgrs_tile,
    Delta_Date,
    year,
    month,
    month_name,
    target_modis_sif,
    ww_pct,
    Latitude,
    Longitude,
    any_of(c(
      "Quality_Flag",
      "Metadata.MeasurementMode",
      "source_file",
      "file_id"
    ))
  )

write_table_pair(
  assigned_export,
  file.path(output_dir, "accepted_raw_sif_assigned_to_nuts3")
)

raw_all <- summarize_raw_sif(
  assigned_sif = assigned_sif,
  baseline_name = "all_accepted"
)

assigned_sif_ww10 <- assigned_sif %>%
  filter(
    is.finite(ww_pct),
    ww_pct >= winter_wheat_threshold
  )

message(
  "Assigned rows with ww_pct >= ",
  winter_wheat_threshold,
  ": ",
  format(nrow(assigned_sif_ww10), big.mark = ",")
)

raw_ww10 <- summarize_raw_sif(
  assigned_sif = assigned_sif_ww10,
  baseline_name = "ww_pct_ge_0.10"
)

raw_all_wide <- monthly_to_wide(
  monthly_data = raw_all$monthly,
  prefix = "raw_all",
  skeleton = model_skeleton
)

raw_ww10_wide <- monthly_to_wide(
  monthly_data = raw_ww10$monthly,
  prefix = "raw_ww10",
  skeleton = model_skeleton
)

write_table_pair(
  raw_all$daily,
  file.path(output_dir, "nuts3_raw_mixed_sif_all_accepted_daily")
)
write_table_pair(
  raw_all$monthly,
  file.path(output_dir, "nuts3_raw_mixed_sif_all_accepted_long")
)
write_table_pair(
  raw_all_wide,
  file.path(output_dir, "nuts3_raw_mixed_sif_all_accepted_wide")
)

write_table_pair(
  raw_ww10$daily,
  file.path(output_dir, "nuts3_raw_mixed_sif_ww10_daily")
)
write_table_pair(
  raw_ww10$monthly,
  file.path(output_dir, "nuts3_raw_mixed_sif_ww10_long")
)
write_table_pair(
  raw_ww10_wide,
  file.path(output_dir, "nuts3_raw_mixed_sif_ww10_wide")
)

coverage_summary <- bind_rows(
  raw_all$monthly,
  raw_ww10$monthly
) %>%
  group_by(baseline, month, month_name) %>%
  summarise(
    n_region_years_with_sif = n(),
    total_footprints = sum(n_raw_footprints),
    total_dates = sum(n_raw_dates),
    median_footprints_per_region_year = median(n_raw_footprints),
    median_dates_per_region_year = median(n_raw_dates),
    .groups = "drop"
  )

write_table_pair(
  coverage_summary,
  file.path(output_dir, "raw_mixed_sif_coverage_summary")
)

if (file.exists(crop_pure_wide_csv)) {
  crop_pure_wide <- readr::read_csv(
    crop_pure_wide_csv,
    show_col_types = FALSE
  ) %>%
    mutate(year = as.integer(year))

  duplicate_crop_pure_keys <- crop_pure_wide %>%
    count(nuts_id, mgrs_tile, year) %>%
    filter(n > 1)

  if (nrow(duplicate_crop_pure_keys) > 0) {
    stop(
      "The crop-pure predictor table contains duplicate ",
      "nuts_id/mgrs_tile/year keys."
    )
  }

  comparison_data <- model_skeleton %>%
    left_join(
      crop_pure_wide %>% select(-any_of("nuts3")),
      by = c("nuts_id", "mgrs_tile", "year")
    ) %>%
    left_join(
      raw_all_wide %>%
        select(-any_of(c("nuts3", "nuts3_area_in_mgrs_pct"))),
      by = c("nuts_id", "mgrs_tile", "year")
    ) %>%
    left_join(
      raw_ww10_wide %>%
        select(-any_of(c("nuts3", "nuts3_area_in_mgrs_pct"))),
      by = c("nuts_id", "mgrs_tile", "year")
    )

  write_table_pair(
    comparison_data,
    file.path(
      output_dir,
      "nuts3_crop_pure_and_raw_mixed_sif_predictors"
    )
  )

  message(
    "Combined crop-pure/raw comparison rows: ",
    nrow(comparison_data)
  )
} else {
  warning(
    "Crop-pure predictor CSV was not found, so the combined comparison ",
    "table was not written: ",
    crop_pure_wide_csv
  )
}

message("")
message("Raw mixed-pixel SIF preparation complete.")
message("Output directory: ", output_dir)
print(coverage_summary, n = Inf)
