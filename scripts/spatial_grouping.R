library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(leaflet)
library(htmlwidgets)

# Fixed-grid spatial aggregation of OCO-2 SIF using the original MGRS tile
# assignment already stored in the input CSV. No geometry from mgrs_de.rds is
# used. Each row belongs to exactly one input tile, one date and one mode.
#
# A Sentinel product is 5,490 x 5,490 pixels at 20 m (109,800 m per side).
# The modeling domain starts at its northwest corner and uses 27 x 27 cells of
# 4,000 m (200 x 200 pixels), covering 108,000 m. This excludes the final
# 1,800 m (90 pixels) along the east and south raster edges.

sf::sf_use_s2(FALSE)

input_csv <- "data/9tiles_2_7_inrange_M01_cnn.csv"
mgrs_tif_dir <- "data/temp_data/mgrs_tifs"
output_dir <- "data/sentinel2_spatial_aggregation_4000m_original_tiles"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

aggregation_size_m <- 4000L
cells_per_side <- 27L
domain_size_m <- aggregation_size_m * cells_per_side
sentinel_pixel_size_m <- 20L
expected_raster_size_m <- 109800L
cell_pixels <- aggregation_size_m / sentinel_pixel_size_m
minimum_primary_footprints <- 3L
minimum_strict_footprints <- 5L

# Leaflet defaults. Leave NA to choose the tile with the most n >= 3 groups.
leaflet_tile_t <- NA_character_
leaflet_sample_groups <- 10L
leaflet_seed <- 42L

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

normalize_mgrs_tile <- function(x) {
  tile <- stringr::str_to_upper(stringr::str_trim(as.character(x)))
  tile <- stringr::str_remove(tile, "^T")
  if_else(is.na(tile) | tile == "", NA_character_, paste0("T", tile))
}

extract_mgrs_tile_from_path <- function(path) {
  stringr::str_extract(as.character(path), "T[0-9]{2}[A-Z]{3}")
}

make_sif_polygon <- function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
  st_polygon(list(rbind(
    c(lon1, lat1),
    c(lon2, lat2),
    c(lon3, lat3),
    c(lon4, lat4),
    c(lon1, lat1)
  )))
}

make_square_polygon <- function(xmin, ymin, xmax, ymax) {
  st_polygon(list(rbind(
    c(xmin, ymin),
    c(xmax, ymin),
    c(xmax, ymax),
    c(xmin, ymax),
    c(xmin, ymin)
  )))
}

collapse_values <- function(x) {
  values <- sort(unique(na.omit(as.character(x))))
  paste(values, collapse = ";")
}

build_tile_table <- function(tif_dir) {
  tif_paths <- list.files(tif_dir, pattern = "[.]tif$", full.names = TRUE)

  if (length(tif_paths) == 0) {
    stop("No Sentinel reference tifs found in: ", tif_dir)
  }

  tile_table <- map_dfr(tif_paths, function(tif_path) {
    raster <- terra::rast(tif_path)
    raster_extent <- terra::ext(raster)
    raster_resolution <- terra::res(raster)

    tibble(
      mgrs_tile_t = extract_mgrs_tile_from_path(basename(tif_path)),
      sentinel_tif_path = tif_path,
      raster_crs = terra::crs(raster),
      raster_xmin = raster_extent$xmin,
      raster_xmax = raster_extent$xmax,
      raster_ymin = raster_extent$ymin,
      raster_ymax = raster_extent$ymax,
      xres = raster_resolution[[1]],
      yres = raster_resolution[[2]],
      raster_ncol = terra::ncol(raster),
      raster_nrow = terra::nrow(raster)
    )
  }) %>%
    filter(!is.na(mgrs_tile_t)) %>%
    arrange(mgrs_tile_t, sentinel_tif_path) %>%
    distinct(mgrs_tile_t, .keep_all = TRUE)

  if (nrow(tile_table) == 0) {
    stop("Could not extract any T-prefixed MGRS tile IDs from tif names.")
  }

  reference_crs <- st_crs(tile_table$raster_crs[[1]])
  same_crs <- map_lgl(tile_table$raster_crs, ~ st_crs(.x) == reference_crs)

  if (!all(same_crs)) {
    stop("Sentinel reference tifs do not all use the same CRS.")
  }
  if (any(abs(tile_table$xres - sentinel_pixel_size_m) > 1e-8) ||
      any(abs(tile_table$yres - sentinel_pixel_size_m) > 1e-8)) {
    stop("Expected every Sentinel reference tif to use 20 m pixels.")
  }
  if (any(tile_table$raster_xmax - tile_table$raster_xmin < domain_size_m) ||
      any(tile_table$raster_ymax - tile_table$raster_ymin < domain_size_m)) {
    stop("At least one Sentinel raster is smaller than 108,000 m.")
  }
  if (any(abs(
        (tile_table$raster_xmax - tile_table$raster_xmin) -
          expected_raster_size_m
      ) > 1e-8) ||
      any(abs(
        (tile_table$raster_ymax - tile_table$raster_ymin) -
          expected_raster_size_m
      ) > 1e-8)) {
    stop("Expected every Sentinel reference raster to span 109,800 m.")
  }

  list(table = tile_table, crs = reference_crs)
}

if (aggregation_size_m %% sentinel_pixel_size_m != 0) {
  stop("The aggregation size must be divisible by the 20 m pixel size.")
}
if (cell_pixels != as.integer(cell_pixels)) {
  stop("The aggregation cell does not contain an integer number of pixels.")
}

message("Reading Sentinel reference tifs...")
tile_result <- build_tile_table(mgrs_tif_dir)
tile_info <- tile_result$table
sentinel_crs <- tile_result$crs

message("Sentinel tiles: ", paste(tile_info$mgrs_tile_t, collapse = ", "))

# Anchor every 108,000 m domain at the northwest corner of its 109,800 m
# Sentinel raster. Columns progress west-to-east; rows progress north-to-south.
tile_grid_info <- tile_info %>%
  mutate(
    domain_xmin = raster_xmin,
    domain_xmax = raster_xmin + domain_size_m,
    domain_ymax = raster_ymax,
    domain_ymin = raster_ymax - domain_size_m,
    domain_size_m = domain_size_m,
    aggregation_size_m = aggregation_size_m,
    aggregation_km = aggregation_size_m / 1000,
    cells_per_side = cells_per_side,
    cell_pixels = as.integer(cell_pixels),
    excluded_east_m = raster_xmax - domain_xmax,
    excluded_south_m = domain_ymin - raster_ymin,
    excluded_east_pixels = excluded_east_m / xres,
    excluded_south_pixels = excluded_south_m / yres
  )

if (any(tile_grid_info$excluded_east_m < 0) ||
    any(tile_grid_info$excluded_south_m < 0)) {
  stop("The northwest-anchored domain extends outside a Sentinel raster.")
}
if (any(abs(tile_grid_info$excluded_east_m - 1800) > 1e-8) ||
    any(abs(tile_grid_info$excluded_south_m - 1800) > 1e-8)) {
  stop("Expected exactly 1,800 m to be excluded at the east and south edges.")
}
if (any(abs(tile_grid_info$excluded_east_pixels -
            round(tile_grid_info$excluded_east_pixels)) > 1e-8) ||
    any(abs(tile_grid_info$excluded_south_pixels -
            round(tile_grid_info$excluded_south_pixels)) > 1e-8)) {
  stop("The excluded raster strips are not aligned to the 20 m pixel grid.")
}

write_csv(
  tile_grid_info,
  file.path(output_dir, "fixed_grid_4000m_tile_domains.csv")
)

message("Reading SIF rows: ", input_csv)
df <- read_csv(
  input_csv,
  col_types = cols(
    Delta_Date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

required_cols <- c(
  "Delta_Date", "Metadata.MeasurementMode", "target_modis_sif",
  "mgrs_tile", "product_path", corner_cols
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

message("Rows read: ", nrow(df))

df_prepared <- df %>%
  mutate(
    sif_row_id = row_number(),
    source_csv_row = row_number(),
    Delta_Date = as.Date(Delta_Date),
    sif_year = year(Delta_Date),
    sif_month = month(Delta_Date),
    sif_doy = yday(Delta_Date),
    measurement_mode = as.integer(.data[["Metadata.MeasurementMode"]]),
    target_modis_sif = as.numeric(target_modis_sif),
    mgrs_tile_original = stringr::str_remove(
      normalize_mgrs_tile(mgrs_tile),
      "^T"
    ),
    mgrs_tile_t = normalize_mgrs_tile(mgrs_tile),
    sentinel_source_tile_t = extract_mgrs_tile_from_path(product_path),
    sentinel_source_tile = stringr::str_remove(
      sentinel_source_tile_t,
      "^T"
    ),
    product_path_exists = file.exists(as.character(product_path)),
    source_matches_input_tile = case_when(
      is.na(sentinel_source_tile_t) | is.na(mgrs_tile_t) ~ NA,
      TRUE ~ sentinel_source_tile_t == mgrs_tile_t
    ),
    across(all_of(corner_cols), as.numeric)
  ) %>%
  filter(
    !is.na(Delta_Date),
    !is.na(measurement_mode),
    is.finite(target_modis_sif),
    if_all(all_of(corner_cols), ~ is.finite(.x))
  )

if ("final_check_modis_sif" %in% names(df_prepared)) {
  df_prepared <- df_prepared %>%
    filter(final_check_modis_sif == "accept")
}

# Keep downstream summaries stable if a nonessential descriptive column is
# absent from a future version of the input table.
if (!"sif_area_km2_evi" %in% names(df_prepared)) {
  df_prepared$sif_area_km2_evi <- NA_real_
}
if (!"state" %in% names(df_prepared)) {
  df_prepared$state <- NA_character_
}
if (!"hzs" %in% names(df_prepared)) {
  df_prepared$hzs <- NA_character_
}

sif_sf <- df_prepared %>%
  mutate(
    geometry = pmap(
      list(
        Lon_corner1, Lat_corner1,
        Lon_corner2, Lat_corner2,
        Lon_corner3, Lat_corner3,
        Lon_corner4, Lat_corner4
      ),
      make_sif_polygon
    )
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid() %>%
  st_transform(sentinel_crs)

sif_centroids <- suppressWarnings(st_centroid(sif_sf))
centroid_xy <- st_coordinates(sif_centroids)

sif_sf <- sif_sf %>%
  mutate(
    centroid_x = centroid_xy[, "X"],
    centroid_y = centroid_xy[, "Y"],
    input_tile_available = mgrs_tile_t %in% tile_info$mgrs_tile_t,
    keep_for_grouping = input_tile_available & product_path_exists
  )

row_audit <- sif_sf %>%
  st_drop_geometry() %>%
  count(
    mgrs_tile_t,
    sentinel_source_tile_t,
    input_tile_available,
    product_path_exists,
    source_matches_input_tile,
    keep_for_grouping,
    name = "n_sif_rows"
  ) %>%
  arrange(mgrs_tile_t, sentinel_source_tile_t)

write_csv(
  row_audit,
  file.path(output_dir, "original_tile_row_audit.csv")
)

n_unavailable_tile <- sum(!sif_sf$input_tile_available)
n_missing_product <- sum(!sif_sf$product_path_exists)
n_source_mismatch <- sum(
  sif_sf$source_matches_input_tile %in% FALSE,
  na.rm = TRUE
)

message("Rows with unavailable input tiles: ", n_unavailable_tile)
message("Rows with missing product paths: ", n_missing_product)
message("Rows where product tile differs from input mgrs_tile: ", n_source_mismatch)

study_sif_sf <- sif_sf %>%
  filter(keep_for_grouping)

if (nrow(study_sif_sf) == 0) {
  stop("No rows remain after input-tile and product-path checks.")
}

sif_for_grid <- study_sif_sf %>%
  st_drop_geometry() %>%
  inner_join(tile_grid_info, by = "mgrs_tile_t")

message("Assigning centroids to northwest-anchored 4,000 m cells...")
all_assignments <- sif_for_grid %>%
  mutate(
    # Row zero is the northernmost row; column zero is the westernmost column.
    cell_col = floor((centroid_x - domain_xmin) / aggregation_size_m),
    cell_row = floor((domain_ymax - centroid_y) / aggregation_size_m),
    centroid_inside_domain =
      cell_col >= 0 & cell_col < cells_per_side &
      cell_row >= 0 & cell_row < cells_per_side,
    cell_xmin = domain_xmin + cell_col * aggregation_size_m,
    cell_xmax = cell_xmin + aggregation_size_m,
    cell_ymax = domain_ymax - cell_row * aggregation_size_m,
    cell_ymin = cell_ymax - aggregation_size_m,
    cell_id = paste0(
      "s2_fixed_4000m_", mgrs_tile_t,
      "_r", str_pad(cell_row, width = 2, pad = "0"),
      "_c", str_pad(cell_col, width = 2, pad = "0")
    ),
    aggregation_id = paste0(
      cell_id,
      "_", format(Delta_Date, "%Y%m%d"),
      "_m", measurement_mode
    )
  )

domain_assignment_summary <- all_assignments %>%
  count(
    mgrs_tile_t,
    centroid_inside_domain,
    name = "n_sif_rows"
  ) %>%
  arrange(mgrs_tile_t, desc(centroid_inside_domain))

write_csv(
  domain_assignment_summary,
  file.path(output_dir, "fixed_grid_4000m_domain_assignment_summary.csv")
)

optional_assignment_cols <- c(
  "Daily_SIF_757nm", "Daily_SIF_771nm", "final_check_modis_sif",
  "sif_area_km2_evi", "state", "hzs", "Latitude", "Longitude"
)

fixed_grid_assignments <- all_assignments %>%
  filter(centroid_inside_domain) %>%
  select(
    aggregation_id,
    cell_id,
    aggregation_size_m,
    aggregation_km,
    cell_pixels,
    cells_per_side,
    cell_row,
    cell_col,
    cell_xmin,
    cell_ymin,
    cell_xmax,
    cell_ymax,
    domain_xmin,
    domain_ymin,
    domain_xmax,
    domain_ymax,
    excluded_east_m,
    excluded_south_m,
    sif_row_id,
    source_csv_row,
    Delta_Date,
    sif_year,
    sif_month,
    sif_doy,
    measurement_mode,
    mgrs_tile,
    mgrs_tile_original,
    mgrs_tile_t,
    sentinel_source_tile,
    sentinel_source_tile_t,
    source_matches_input_tile,
    sentinel_tif_path,
    product_path,
    target_modis_sif,
    any_of(optional_assignment_cols),
    centroid_x,
    centroid_y,
    all_of(corner_cols)
  ) %>%
  arrange(
    mgrs_tile_t,
    Delta_Date,
    measurement_mode,
    cell_row,
    cell_col,
    sif_row_id
  )

# A source row must occur once at most. Raster overlap cannot duplicate a row
# because its original mgrs_tile determines the only grid it can enter.
duplicate_assignments <- fixed_grid_assignments %>%
  count(sif_row_id) %>%
  filter(n != 1)

if (nrow(duplicate_assignments) > 0) {
  stop("At least one SIF row was assigned to more than one 4,000 m cell.")
}

group_integrity <- fixed_grid_assignments %>%
  group_by(aggregation_id) %>%
  summarise(
    n_dates = n_distinct(Delta_Date),
    n_modes = n_distinct(measurement_mode),
    n_tiles = n_distinct(mgrs_tile_t),
    .groups = "drop"
  )

if (any(group_integrity$n_dates != 1) ||
    any(group_integrity$n_modes != 1) ||
    any(group_integrity$n_tiles != 1)) {
  stop("An aggregate group contains multiple dates, modes or input tiles.")
}

message("Summarising aggregated SIF targets...")
fixed_grid_manifest <- fixed_grid_assignments %>%
  group_by(
    aggregation_id,
    cell_id,
    aggregation_size_m,
    aggregation_km,
    cell_pixels,
    cells_per_side,
    cell_row,
    cell_col,
    cell_xmin,
    cell_ymin,
    cell_xmax,
    cell_ymax,
    mgrs_tile_original,
    mgrs_tile_t,
    Delta_Date,
    sif_year,
    sif_month,
    sif_doy,
    measurement_mode,
    sentinel_tif_path
  ) %>%
  summarise(
    n_footprints = n(),
    aggregated_target_modis_sif = mean(target_modis_sif),
    median_target_modis_sif = median(target_modis_sif),
    min_target_modis_sif = min(target_modis_sif),
    max_target_modis_sif = max(target_modis_sif),
    sd_target_modis_sif = if (n() > 1) sd(target_modis_sif) else NA_real_,
    se_target_modis_sif = if (n() > 1) {
      sd(target_modis_sif) / sqrt(n())
    } else {
      NA_real_
    },
    mean_sif_area_km2 = mean(sif_area_km2_evi, na.rm = TRUE),
    total_sif_area_km2 = sum(sif_area_km2_evi, na.rm = TRUE),
    states = collapse_values(state),
    hzs_values = collapse_values(hzs),
    n_sentinel_source_tiles = n_distinct(
      sentinel_source_tile_t,
      na.rm = TRUE
    ),
    sentinel_source_tiles = collapse_values(sentinel_source_tile_t),
    n_product_paths = n_distinct(product_path, na.rm = TRUE),
    product_paths = collapse_values(product_path),
    all_source_tiles_match_input = all(
      coalesce(source_matches_input_tile, FALSE)
    ),
    sif_row_ids = paste(sif_row_id, collapse = ","),
    source_csv_rows = paste(source_csv_row, collapse = ","),
    .groups = "drop"
  ) %>%
  mutate(
    eligible_n3 = n_footprints >= minimum_primary_footprints,
    eligible_n5 = n_footprints >= minimum_strict_footprints,
    cell_area_km2 = (aggregation_size_m / 1000)^2
  ) %>%
  arrange(
    mgrs_tile_t,
    Delta_Date,
    measurement_mode,
    cell_row,
    cell_col
  )

fixed_grid_assignments <- fixed_grid_assignments %>%
  left_join(
    fixed_grid_manifest %>%
      select(
        aggregation_id,
        n_footprints,
        aggregated_target_modis_sif,
        sd_target_modis_sif,
        se_target_modis_sif,
        eligible_n3,
        eligible_n5
      ),
    by = "aggregation_id"
  )

write_csv(
  fixed_grid_assignments,
  file.path(output_dir, "fixed_grid_4000m_sif_assignments.csv")
)

write_csv(
  fixed_grid_manifest,
  file.path(output_dir, "fixed_grid_4000m_aggregate_manifest.csv")
)

fixed_grid_manifest_sf <- fixed_grid_manifest %>%
  mutate(
    geometry = pmap(
      list(cell_xmin, cell_ymin, cell_xmax, cell_ymax),
      make_square_polygon
    )
  ) %>%
  st_as_sf(crs = sentinel_crs)

saveRDS(
  fixed_grid_manifest_sf,
  file.path(output_dir, "fixed_grid_4000m_aggregate_polygons.rds")
)

count_distribution <- fixed_grid_manifest %>%
  count(n_footprints, name = "n_aggregate_groups") %>%
  arrange(n_footprints) %>%
  mutate(
    pct_aggregate_groups = n_aggregate_groups / sum(n_aggregate_groups),
    groups_ge_n = rev(cumsum(rev(n_aggregate_groups))),
    pct_groups_ge_n = groups_ge_n / sum(n_aggregate_groups)
  )

write_csv(
  count_distribution,
  file.path(output_dir, "fixed_grid_4000m_footprint_count_distribution.csv")
)

threshold_summary <- fixed_grid_manifest %>%
  summarise(
    aggregation_size_m = first(aggregation_size_m),
    aggregation_km = first(aggregation_km),
    cell_pixels = first(cell_pixels),
    cells_per_side = first(cells_per_side),
    n_aggregate_groups = n(),
    n_footprints_assigned = sum(n_footprints),
    mean_footprints_per_group = mean(n_footprints),
    median_footprints_per_group = median(n_footprints),
    max_footprints_per_group = max(n_footprints),
    groups_ge_2 = sum(n_footprints >= 2),
    groups_ge_3 = sum(n_footprints >= 3),
    groups_ge_5 = sum(n_footprints >= 5),
    groups_ge_10 = sum(n_footprints >= 10),
    footprints_in_groups_ge_3 = sum(n_footprints[n_footprints >= 3]),
    footprints_in_groups_ge_5 = sum(n_footprints[n_footprints >= 5]),
    pct_groups_ge_3 = groups_ge_3 / n_aggregate_groups,
    pct_groups_ge_5 = groups_ge_5 / n_aggregate_groups,
    pct_footprints_in_groups_ge_3 =
      footprints_in_groups_ge_3 / n_footprints_assigned,
    pct_footprints_in_groups_ge_5 =
      footprints_in_groups_ge_5 / n_footprints_assigned
  )

write_csv(
  threshold_summary,
  file.path(output_dir, "fixed_grid_4000m_threshold_summary.csv")
)

balance_summary <- bind_rows(
  fixed_grid_manifest %>%
    filter(eligible_n5) %>%
    count(value = as.character(sif_year), name = "n_groups") %>%
    mutate(dimension = "year"),
  fixed_grid_manifest %>%
    filter(eligible_n5) %>%
    count(value = str_pad(sif_month, 2, pad = "0"), name = "n_groups") %>%
    mutate(dimension = "month"),
  fixed_grid_manifest %>%
    filter(eligible_n5) %>%
    count(value = as.character(measurement_mode), name = "n_groups") %>%
    mutate(dimension = "measurement_mode"),
  fixed_grid_manifest %>%
    filter(eligible_n5) %>%
    count(value = mgrs_tile_t, name = "n_groups") %>%
    mutate(dimension = "mgrs_tile")
) %>%
  select(dimension, value, n_groups) %>%
  arrange(dimension, value)

write_csv(
  balance_summary,
  file.path(output_dir, "fixed_grid_4000m_n5_balance_summary.csv")
)

message("\nTraining-sample availability at 4 km:")
print(threshold_summary, n = Inf, width = Inf)

message("\nWrote corrected fixed-grid diagnostics to: ", output_dir)
message("  - fixed_grid_4000m_sif_assignments.csv")
message("  - fixed_grid_4000m_aggregate_manifest.csv")
message("  - fixed_grid_4000m_aggregate_polygons.rds")
message("  - fixed_grid_4000m_footprint_count_distribution.csv")
message("  - fixed_grid_4000m_threshold_summary.csv")
message("  - fixed_grid_4000m_n5_balance_summary.csv")

# ---------------------------------------------------------------------------
# Leaflet preview: ten sampled n >= 3 groups from one input tile.

preview_candidates <- fixed_grid_manifest %>%
  filter(eligible_n3)

if (nrow(preview_candidates) == 0) {
  warning("No n >= 3 aggregate groups are available for the Leaflet preview.")
} else {
  if (is.na(leaflet_tile_t)) {
    leaflet_tile_t <- preview_candidates %>%
      count(mgrs_tile_t, name = "n_groups") %>%
      arrange(desc(n_groups), mgrs_tile_t) %>%
      slice(1) %>%
      pull(mgrs_tile_t)
  }

  if (!leaflet_tile_t %in% preview_candidates$mgrs_tile_t) {
    warning("No n >= 3 groups found for Leaflet tile: ", leaflet_tile_t)
  } else {
    set.seed(leaflet_seed)

    tile_preview_candidates <- preview_candidates %>%
      filter(mgrs_tile_t == leaflet_tile_t)

    n_preview_groups <- min(
      leaflet_sample_groups,
      nrow(tile_preview_candidates)
    )

    sampled_groups <- tile_preview_candidates %>%
      slice_sample(n = n_preview_groups)

    sampled_cells <- fixed_grid_manifest_sf %>%
      filter(aggregation_id %in% sampled_groups$aggregation_id) %>%
      mutate(
        cell_label = paste0(
          "<b>4 km aggregate</b>",
          "<br>ID: ", aggregation_id,
          "<br>Date: ", Delta_Date,
          "<br>Mode: ", measurement_mode,
          "<br>Footprints: ", n_footprints,
          "<br>Mean SIF: ", round(aggregated_target_modis_sif, 4),
          "<br>SD: ", round(sd_target_modis_sif, 4),
          "<br>SE: ", round(se_target_modis_sif, 4)
        )
      ) %>%
      st_transform(4326)

    sampled_assignment_rows <- fixed_grid_assignments %>%
      filter(aggregation_id %in% sampled_groups$aggregation_id) %>%
      select(
        aggregation_id,
        n_footprints,
        aggregated_target_modis_sif,
        sif_row_id,
        target_modis_sif,
        Delta_Date,
        measurement_mode,
        mgrs_tile_t,
        any_of(c("state", "hzs"))
      )

    sampled_footprints <- study_sif_sf %>%
      select(sif_row_id) %>%
      inner_join(sampled_assignment_rows, by = "sif_row_id") %>%
      mutate(
        footprint_label = paste0(
          "<b>SIF row: ", sif_row_id, "</b>",
          "<br>Date: ", Delta_Date,
          "<br>Mode: ", measurement_mode,
          "<br>Observed SIF: ", round(target_modis_sif, 4),
          "<br>Group mean: ", round(aggregated_target_modis_sif, 4)
        )
      ) %>%
      st_transform(4326)

    tile_domain_sf <- tile_grid_info %>%
      filter(mgrs_tile_t == leaflet_tile_t) %>%
      mutate(
        geometry = pmap(
          list(domain_xmin, domain_ymin, domain_xmax, domain_ymax),
          make_square_polygon
        )
      ) %>%
      st_as_sf(crs = sentinel_crs) %>%
      st_transform(4326)

    target_palette <- colorNumeric(
      palette = "viridis",
      domain = sampled_footprints$target_modis_sif,
      na.color = "#cccccc"
    )

    leaflet_map <- leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Basemap") %>%
      addPolygons(
        data = tile_domain_sf,
        fill = FALSE,
        color = "#111111",
        weight = 3,
        opacity = 1,
        label = ~paste0(mgrs_tile_t, " 108 km domain"),
        group = "108 km tile domain"
      ) %>%
      addPolygons(
        data = sampled_cells,
        fillColor = "#d95f0e",
        fillOpacity = 0.14,
        color = "#d95f0e",
        weight = 1.6,
        opacity = 0.95,
        label = ~lapply(cell_label, htmltools::HTML),
        group = "4 km cells"
      ) %>%
      addPolygons(
        data = sampled_footprints,
        fillColor = ~target_palette(target_modis_sif),
        fillOpacity = 0.7,
        color = "#444444",
        weight = 0.7,
        opacity = 0.85,
        label = ~lapply(footprint_label, htmltools::HTML),
        group = "SIF footprints"
      ) %>%
      addLegend(
        position = "bottomright",
        pal = target_palette,
        values = sampled_footprints$target_modis_sif,
        title = "target_modis_sif",
        opacity = 0.85
      ) %>%
      addLayersControl(
        overlayGroups = c(
          "108 km tile domain",
          "4 km cells",
          "SIF footprints"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )

    leaflet_out <- file.path(
      output_dir,
      paste0("fixed_grid_4000m_leaflet_preview_", leaflet_tile_t, ".html")
    )

    htmlwidgets::saveWidget(leaflet_map, leaflet_out, selfcontained = TRUE)
    message("Wrote Leaflet preview: ", leaflet_out)
  }
}
