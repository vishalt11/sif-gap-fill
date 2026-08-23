library(tidyverse)
library(terra)
library(sf)
library(lubridate)
library(dbscan)
library(leaflet)
library(htmlwidgets)

# Density-centred spatial aggregation of OCO-2 SIF.
#
# DBSCAN is run once per tile/date/mode stratum to identify dense track
# neighbourhoods. Each neighbourhood is ordered along its principal spatial
# direction and split with a fast forward chunk scan. Every accepted sample is
# represented by an explicit 4,000 m x 4,000 m square aligned to the 20 m
# Sentinel-2 grid. Only footprints belonging to the most frequent land-cover
# class inside that square contribute to the aggregate target. All footprints
# encountered by an accepted square are consumed, which prevents one source
# sounding from being used in multiple aggregate targets.

sf::sf_use_s2(FALSE)

input_csv <- paste0(
  "data/main_sif_data/",
  "9tiles_2_7_M01_QF01_inoutrange_PARrm_BKR_landcover.csv"
)
mgrs_tif_dir <- "data/temp_data/mgrs_tifs"
output_dir <- "data/sentinel2_spatial_aggregation_density_4000m_landcover"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

window_size_m <- 4000L
half_window_m <- window_size_m / 2
sentinel_pixel_size_m <- 20L
window_pixels <- window_size_m / sentinel_pixel_size_m

# eps connects nearby soundings into candidate track neighbourhoods. DBSCAN
# clusters can be longer than 4 km, so their extents are never used directly as
# modelling windows. The explicit square-window check below controls support.
dbscan_eps_m <- 1800
minimum_footprints <- 4L

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

make_sif_polygon <- function(
    lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
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

safe_mean <- function(x) {
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) NA_real_ else mean(values)
}

safe_min <- function(x) {
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) NA_real_ else min(values)
}

safe_max <- function(x) {
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) NA_real_ else max(values)
}

modal_value <- function(x) {
  values <- as.character(x)
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(NA_character_)
  }

  counts <- sort(table(values), decreasing = TRUE)
  tied <- names(counts)[counts == max(counts)]
  sort(tied)[[1]]
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
    stop("Could not extract MGRS tile IDs from the reference tif names.")
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
  if (any(tile_table$raster_xmax - tile_table$raster_xmin < window_size_m) ||
      any(tile_table$raster_ymax - tile_table$raster_ymin < window_size_m)) {
    stop("At least one Sentinel reference tif is smaller than 4 km.")
  }

  list(table = tile_table, crs = reference_crs)
}

# Snap a candidate square to the reference tile's 20 m grid and, where needed,
# shift it inward so that the entire 4 km window remains inside the raster.
snap_window_to_tile <- function(center_x, center_y, tile_row) {
  xmin_raw <- center_x - half_window_m
  ymax_raw <- center_y + half_window_m

  xmin <- tile_row$raster_xmin +
    round((xmin_raw - tile_row$raster_xmin) / sentinel_pixel_size_m) *
      sentinel_pixel_size_m
  ymax <- tile_row$raster_ymax -
    round((tile_row$raster_ymax - ymax_raw) / sentinel_pixel_size_m) *
      sentinel_pixel_size_m

  xmin <- min(
    max(xmin, tile_row$raster_xmin),
    tile_row$raster_xmax - window_size_m
  )
  ymax <- min(
    max(ymax, tile_row$raster_ymin + window_size_m),
    tile_row$raster_ymax
  )

  tibble(
    cell_xmin = xmin,
    cell_ymin = ymax - window_size_m,
    cell_xmax = xmin + window_size_m,
    cell_ymax = ymax,
    center_x = xmin + half_window_m,
    center_y = ymax - half_window_m
  )
}

inside_window <- function(tbl, window_row) {
  tbl$centroid_x >= window_row$cell_xmin &
    tbl$centroid_x < window_row$cell_xmax &
    tbl$centroid_y > window_row$cell_ymin &
    tbl$centroid_y <= window_row$cell_ymax
}

land_cover_count_string <- function(tbl) {
  tbl %>%
    count(land_cover, land_cover_class, name = "n") %>%
    arrange(desc(n), land_cover, land_cover_class) %>%
    transmute(value = paste0(land_cover, ":", land_cover_class, "=", n)) %>%
    pull(value) %>%
    paste(collapse = ";")
}

# Order a DBSCAN track neighbourhood along its dominant spatial direction.
# This follows the efficient PCA/chunking strategy used by the earlier
# multi-footprint diagnostic instead of testing a window around every point.
add_pca_track_score <- function(cluster_tbl) {
  coordinates <- as.matrix(cluster_tbl[, c("centroid_x", "centroid_y")])

  if (nrow(coordinates) < 2) {
    return(cluster_tbl %>% mutate(track_score = centroid_x))
  }

  track_score <- tryCatch(
    stats::prcomp(
      coordinates,
      center = TRUE,
      scale. = FALSE
    )$x[, 1],
    error = function(e) coordinates[, 1]
  )

  cluster_tbl %>%
    mutate(track_score = track_score) %>%
    arrange(track_score, sif_row_id)
}

# Construct one grid-aligned 4 km window from centroid bounds. A zero-row
# result means that the points cannot all fit inside one valid window.
window_from_bounds <- function(
    xmin, xmax, ymin, ymax, tile_row) {
  if ((xmax - xmin) >= window_size_m ||
      (ymax - ymin) >= window_size_m) {
    return(tibble())
  }

  window <- snap_window_to_tile(
    center_x = (xmin + xmax) / 2,
    center_y = (ymin + ymax) / 2,
    tile_row = tile_row
  )

  bounds_inside <-
    xmin >= window$cell_xmin[[1]] &&
    xmax < window$cell_xmax[[1]] &&
    ymin > window$cell_ymin[[1]] &&
    ymax <= window$cell_ymax[[1]]

  if (!bounds_inside) tibble() else window
}

# Find the largest leading section of a PCA-ordered track that fits in one
# 4 km window. Bounds are updated incrementally, so this is linear in the
# number of remaining points rather than an all-points-by-all-candidates scan.
largest_fitting_prefix <- function(ordered_tbl, tile_row) {
  if (nrow(ordered_tbl) < minimum_footprints) {
    return(list(n = 0L, window = tibble()))
  }

  xmin <- Inf
  xmax <- -Inf
  ymin <- Inf
  ymax <- -Inf
  best_n <- 0L
  best_window <- tibble()

  for (row_i in seq_len(nrow(ordered_tbl))) {
    xmin <- min(xmin, ordered_tbl$centroid_x[[row_i]])
    xmax <- max(xmax, ordered_tbl$centroid_x[[row_i]])
    ymin <- min(ymin, ordered_tbl$centroid_y[[row_i]])
    ymax <- max(ymax, ordered_tbl$centroid_y[[row_i]])

    if (row_i < minimum_footprints) {
      next
    }

    candidate_window <- window_from_bounds(
      xmin, xmax, ymin, ymax, tile_row
    )
    if (nrow(candidate_window) == 0) {
      break
    }

    best_n <- as.integer(row_i)
    best_window <- candidate_window
  }

  list(n = best_n, window = best_window)
}

empty_assignment_table <- function() {
  tibble()
}

empty_exclusion_table <- function() {
  tibble(
    sif_row_id = integer(),
    exclusion_reason = character(),
    related_aggregation_id = character()
  )
}

cluster_one_dbscan_group <- function(
    cluster_tbl, tile_row, stratum_number) {
  remaining <- add_pca_track_score(cluster_tbl)
  accepted <- list()
  excluded <- list()
  accepted_number <- 0L
  dbscan_cluster <- cluster_tbl$dbscan_cluster[[1]]

  while (nrow(remaining) >= minimum_footprints) {
    fit <- largest_fitting_prefix(remaining, tile_row)

    # If the first four PCA-ordered points cannot share one 4 km support,
    # remove the leading point and retry without rerunning DBSCAN.
    if (fit$n < minimum_footprints) {
      excluded[[length(excluded) + 1L]] <- tibble(
        sif_row_id = remaining$sif_row_id[[1]],
        exclusion_reason = "cannot_form_4km_track_chunk",
        related_aggregation_id = NA_character_
      )
      remaining <- remaining[-1, , drop = FALSE]
      next
    }

    # Avoid leaving a final one-to-three-point remainder where a smaller first
    # chunk can leave at least four points for the following window.
    trailing_n <- nrow(remaining) - fit$n
    if (trailing_n > 0L && trailing_n < minimum_footprints) {
      balanced_n <- fit$n - (minimum_footprints - trailing_n)
      if (balanced_n >= minimum_footprints) {
        balanced_rows <- remaining[seq_len(balanced_n), , drop = FALSE]
        balanced_window <- window_from_bounds(
          min(balanced_rows$centroid_x),
          max(balanced_rows$centroid_x),
          min(balanced_rows$centroid_y),
          max(balanced_rows$centroid_y),
          tile_row
        )
        if (nrow(balanced_window) == 1) {
          fit$n <- as.integer(balanced_n)
          fit$window <- balanced_window
        }
      }
    }

    # Include every still-unassigned point from this DBSCAN neighbourhood that
    # lies inside the selected support, then retain only its majority class.
    in_window <- inside_window(remaining, fit$window)
    window_tbl <- remaining[in_window, , drop = FALSE]
    class_counts <- window_tbl %>%
      count(land_cover, land_cover_class, name = "n") %>%
      arrange(desc(n), land_cover, land_cover_class)
    majority <- class_counts %>% slice(1)
    majority_tbl <- window_tbl %>%
      filter(
        land_cover == majority$land_cover[[1]],
        land_cover_class == majority$land_cover_class[[1]]
      )

    if (nrow(majority_tbl) < minimum_footprints) {
      excluded[[length(excluded) + 1L]] <- tibble(
        sif_row_id = remaining$sif_row_id[[1]],
        exclusion_reason = "fewer_than_four_majority_land_cover_footprints",
        related_aggregation_id = NA_character_
      )
      remaining <- remaining[-1, , drop = FALSE]
      next
    }

    accepted_number <- accepted_number + 1L

    date_string <- format(remaining$Delta_Date[[1]], "%Y%m%d")
    tile_string <- remaining$mgrs_tile_t[[1]]
    mode_value <- remaining$measurement_mode[[1]]
    land_cover_value <- majority$land_cover[[1]]
    aggregation_id <- paste0(
      "s2_density_4000m_", tile_string,
      "_", date_string,
      "_m", mode_value,
      "_s", str_pad(stratum_number, width = 5, pad = "0"),
      "_db", str_pad(dbscan_cluster, width = 3, pad = "0"),
      "_w", str_pad(accepted_number, width = 3, pad = "0"),
      "_lc", land_cover_value
    )

    assigned_ids <- majority_tbl$sif_row_id
    consumed_ids <- window_tbl$sif_row_id

    accepted[[length(accepted) + 1L]] <- remaining %>%
      filter(sif_row_id %in% assigned_ids) %>%
      mutate(
        aggregation_id = aggregation_id,
        cell_id = aggregation_id,
        aggregation_size_m = window_size_m,
        aggregation_km = window_size_m / 1000,
        cell_pixels = as.integer(window_pixels),
        cell_xmin = fit$window$cell_xmin[[1]],
        cell_ymin = fit$window$cell_ymin[[1]],
        cell_xmax = fit$window$cell_xmax[[1]],
        cell_ymax = fit$window$cell_ymax[[1]],
        center_x = fit$window$center_x[[1]],
        center_y = fit$window$center_y[[1]],
        window_n_all_land_covers = nrow(window_tbl),
        majority_land_cover = majority$land_cover[[1]],
        majority_land_cover_class = majority$land_cover_class[[1]],
        excluded_minority_footprints =
          nrow(window_tbl) - nrow(majority_tbl),
        land_cover_majority_fraction =
          nrow(majority_tbl) / nrow(window_tbl),
        n_land_cover_classes_in_window = nrow(class_counts),
        land_cover_counts = land_cover_count_string(window_tbl)
      )

    minority_ids <- setdiff(consumed_ids, assigned_ids)
    if (length(minority_ids) > 0) {
      excluded[[length(excluded) + 1L]] <- tibble(
        sif_row_id = minority_ids,
        exclusion_reason = "minority_land_cover_in_accepted_window",
        related_aggregation_id = aggregation_id
      )
    }

    # Majority rows are assigned and minority rows are excluded. Removing every
    # row encountered by the window guarantees one-use-only supervision.
    remaining <- remaining %>%
      filter(!sif_row_id %in% consumed_ids)
  }

  if (nrow(remaining) > 0) {
    excluded[[length(excluded) + 1L]] <- tibble(
      sif_row_id = remaining$sif_row_id,
      exclusion_reason = "density_noise_or_fewer_than_four_majority_footprints",
      related_aggregation_id = NA_character_
    )
  }

  list(
    assignments = if (length(accepted) == 0) {
      empty_assignment_table()
    } else {
      bind_rows(accepted)
    },
    exclusions = if (length(excluded) == 0) {
      empty_exclusion_table()
    } else {
      bind_rows(excluded)
    }
  )
}

cluster_one_stratum <- function(stratum_tbl, stratum_number) {
  if (nrow(stratum_tbl) < minimum_footprints) {
    return(list(
      assignments = empty_assignment_table(),
      exclusions = tibble(
        sif_row_id = stratum_tbl$sif_row_id,
        exclusion_reason = "stratum_fewer_than_four_footprints",
        related_aggregation_id = NA_character_
      )
    ))
  }

  tile_row <- stratum_tbl %>%
    slice(1) %>%
    select(
      raster_xmin, raster_xmax, raster_ymin, raster_ymax,
      sentinel_tif_path
    )

  coordinates <- as.matrix(
    stratum_tbl[, c("centroid_x", "centroid_y")]
  )
  dbscan_result <- dbscan::dbscan(
    coordinates,
    eps = dbscan_eps_m,
    minPts = minimum_footprints
  )
  labelled <- stratum_tbl %>%
    mutate(dbscan_cluster = dbscan_result$cluster)

  noise <- labelled %>%
    filter(dbscan_cluster == 0L)
  clustered <- labelled %>%
    filter(dbscan_cluster > 0L)

  noise_exclusions <- if (nrow(noise) == 0) {
    empty_exclusion_table()
  } else {
    tibble(
      sif_row_id = noise$sif_row_id,
      exclusion_reason = "dbscan_noise",
      related_aggregation_id = NA_character_
    )
  }

  if (nrow(clustered) == 0) {
    return(list(
      assignments = empty_assignment_table(),
      exclusions = noise_exclusions
    ))
  }

  cluster_results <- clustered %>%
    group_by(dbscan_cluster) %>%
    group_split(.keep = TRUE) %>%
    map(
      ~ cluster_one_dbscan_group(
        .x,
        tile_row = tile_row,
        stratum_number = stratum_number
      )
    )

  assignments <- map_dfr(cluster_results, "assignments")
  exclusions <- bind_rows(
    noise_exclusions,
    map_dfr(cluster_results, "exclusions")
  )

  list(
    assignments = assignments,
    exclusions = exclusions
  )
}

if (window_size_m %% sentinel_pixel_size_m != 0) {
  stop("The 4 km window must contain an integer number of 20 m pixels.")
}

message("Reading Sentinel reference tifs...")
tile_result <- build_tile_table(mgrs_tif_dir)
tile_info <- tile_result$table
sentinel_crs <- tile_result$crs

write_csv(
  tile_info,
  file.path(output_dir, "density_4000m_tile_reference.csv")
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
  "final_check_modis_sif", "Quality_Flag", "date_align", "mgrs_tile",
  "product_path", "BKR10_ID", "BKR_NAME", "land_cover",
  "land_cover_class", corner_cols
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
    Quality_Flag = as.integer(Quality_Flag),
    date_align = stringr::str_to_lower(
      stringr::str_trim(as.character(date_align))
    ),
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
    BKR10_ID = as.character(BKR10_ID),
    BKR_NAME = stringr::str_squish(as.character(BKR_NAME)),
    land_cover = as.integer(land_cover),
    land_cover_class = stringr::str_to_lower(
      stringr::str_squish(as.character(land_cover_class))
    ),
    across(all_of(corner_cols), as.numeric)
  ) %>%
  filter(
    !is.na(Delta_Date),
    !is.na(measurement_mode),
    !is.na(Quality_Flag),
    date_align %in% c("inrange", "outrange"),
    is.finite(target_modis_sif),
    stringr::str_to_lower(final_check_modis_sif) == "accept",
    !is.na(land_cover),
    !is.na(land_cover_class),
    land_cover_class != "",
    if_all(all_of(corner_cols), ~ is.finite(.x))
  )

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
    keep_for_grouping =
      input_tile_available &
      product_path_exists &
      coalesce(source_matches_input_tile, FALSE)
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
  file.path(output_dir, "density_4000m_input_row_audit.csv")
)

study_tbl <- sif_sf %>%
  filter(keep_for_grouping) %>%
  st_drop_geometry() %>%
  inner_join(
    tile_info %>%
      select(
        mgrs_tile_t, sentinel_tif_path,
        raster_xmin, raster_xmax, raster_ymin, raster_ymax
      ),
    by = "mgrs_tile_t"
  )

if (nrow(study_tbl) == 0) {
  stop("No rows remain after tile, product-path and source-tile checks.")
}

# A stratum represents one OCO-2 track subset within one Sentinel-2 tile and
# day. Product path and date alignment are retained as metadata, but they do
# not split otherwise compatible SIF observations into separate strata.
strata <- study_tbl %>%
  arrange(
    mgrs_tile_t, Delta_Date, measurement_mode, sif_row_id
  ) %>%
  group_by(
    mgrs_tile_t, Delta_Date, measurement_mode
  ) %>%
  group_split(.keep = TRUE)

message(
  "Running density-centred 4 km aggregation for ",
  length(strata),
  " tile/date/mode strata..."
)

cluster_results <- imap(
  strata,
  function(stratum_tbl, stratum_i) {
    if (stratum_i == 1L ||
        stratum_i %% 25L == 0L ||
        stratum_i == length(strata)) {
      message(
        "  Stratum ",
        format(stratum_i, big.mark = ","),
        " / ",
        format(length(strata), big.mark = ",")
      )
    }
    cluster_one_stratum(stratum_tbl, as.integer(stratum_i))
  }
)

density_assignments_full <- map_dfr(cluster_results, "assignments")
density_exclusions <- map_dfr(cluster_results, "exclusions")

if (nrow(density_assignments_full) == 0) {
  stop(
    "No density-centred groups were accepted. Consider inspecting the input ",
    "land-cover classes or increasing dbscan_eps_m."
  )
}

duplicate_assignments <- density_assignments_full %>%
  count(sif_row_id) %>%
  filter(n != 1)
if (nrow(duplicate_assignments) > 0) {
  stop("At least one SIF row was assigned to more than one density window.")
}

group_integrity <- density_assignments_full %>%
  group_by(aggregation_id) %>%
  summarise(
    n_dates = n_distinct(Delta_Date),
    n_modes = n_distinct(measurement_mode),
    n_tiles = n_distinct(mgrs_tile_t),
    n_products = n_distinct(product_path),
    n_date_align = n_distinct(date_align),
    n_land_covers = n_distinct(land_cover),
    centroids_inside_window = all(
      centroid_x >= first(cell_xmin) &
        centroid_x < first(cell_xmax) &
        centroid_y > first(cell_ymin) &
        centroid_y <= first(cell_ymax)
    ),
    .groups = "drop"
  )

if (any(group_integrity$n_dates != 1) ||
    any(group_integrity$n_modes != 1) ||
    any(group_integrity$n_tiles != 1) ||
    any(group_integrity$n_land_covers != 1) ||
    any(!group_integrity$centroids_inside_window)) {
  stop("At least one accepted density group failed its integrity checks.")
}

density_manifest <- density_assignments_full %>%
  group_by(
    aggregation_id,
    cell_id,
    aggregation_size_m,
    aggregation_km,
    cell_pixels,
    cell_xmin,
    cell_ymin,
    cell_xmax,
    cell_ymax,
    center_x,
    center_y,
    mgrs_tile_original,
    mgrs_tile_t,
    Delta_Date,
    sif_year,
    sif_month,
    sif_doy,
    measurement_mode,
    sentinel_tif_path,
    majority_land_cover,
    majority_land_cover_class,
    window_n_all_land_covers,
    excluded_minority_footprints,
    land_cover_majority_fraction,
    n_land_cover_classes_in_window,
    land_cover_counts
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
    mean_sif_area_km2 = safe_mean(sif_area_km2_evi),
    min_sif_area_km2 = safe_min(sif_area_km2_evi),
    max_sif_area_km2 = safe_max(sif_area_km2_evi),
    total_sif_area_km2 = sum(sif_area_km2_evi, na.rm = TRUE),
    states = collapse_values(state),
    hzs_values = collapse_values(hzs),
    BKR10_ID_values = collapse_values(BKR10_ID),
    BKR_NAME_values = collapse_values(BKR_NAME),
    majority_BKR10_ID = modal_value(BKR10_ID),
    majority_BKR_NAME = modal_value(BKR_NAME),
    n_BKR10 = n_distinct(BKR10_ID, na.rm = TRUE),
    quality_flag_values = collapse_values(Quality_Flag),
    n_quality_flag_0 = sum(Quality_Flag == 0L),
    n_quality_flag_1 = sum(Quality_Flag == 1L),
    quality_flag_1_fraction = mean(Quality_Flag == 1L),
    date_align_values = collapse_values(date_align),
    n_date_inrange = sum(date_align == "inrange"),
    n_date_outrange = sum(date_align == "outrange"),
    date_outrange_fraction = mean(date_align == "outrange"),
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
    eligible_n4 = n_footprints >= 4L,
    eligible_n5 = n_footprints >= 5L,
    eligible_n6 = n_footprints >= 6L,
    eligible_n8 = n_footprints >= 8L,
    eligible_n10 = n_footprints >= 10L,
    cell_area_km2 = (aggregation_size_m / 1000)^2
  ) %>%
  arrange(
    mgrs_tile_t,
    Delta_Date,
    measurement_mode,
    desc(n_footprints),
    aggregation_id
  )

assignment_optional_cols <- c(
  "Daily_SIF_757nm", "Daily_SIF_771nm", "sif_area_km2_evi",
  "state", "hzs", "Latitude", "Longitude"
)

density_assignments <- density_assignments_full %>%
  select(
    aggregation_id,
    cell_id,
    aggregation_size_m,
    aggregation_km,
    cell_pixels,
    cell_xmin,
    cell_ymin,
    cell_xmax,
    cell_ymax,
    center_x,
    center_y,
    window_n_all_land_covers,
    majority_land_cover,
    majority_land_cover_class,
    excluded_minority_footprints,
    land_cover_majority_fraction,
    n_land_cover_classes_in_window,
    land_cover_counts,
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
    final_check_modis_sif,
    Quality_Flag,
    date_align,
    BKR10_ID,
    BKR_NAME,
    land_cover,
    land_cover_class,
    any_of(assignment_optional_cols),
    centroid_x,
    centroid_y,
    all_of(corner_cols)
  ) %>%
  arrange(aggregation_id, sif_row_id)

manifest_path <- file.path(
  output_dir,
  "density_cluster_4000m_aggregate_manifest.csv"
)
assignments_path <- file.path(
  output_dir,
  "density_cluster_4000m_sif_assignments.csv"
)
exclusions_path <- file.path(
  output_dir,
  "density_cluster_4000m_excluded_sif_rows.csv"
)

write_csv(density_manifest, manifest_path)
write_csv(density_assignments, assignments_path)
write_csv(density_exclusions, exclusions_path)

count_distribution <- density_manifest %>%
  count(n_footprints, name = "n_windows") %>%
  arrange(n_footprints) %>%
  mutate(
    pct_windows = n_windows / sum(n_windows),
    cumulative_windows_ge_n = rev(cumsum(rev(n_windows))),
    pct_windows_ge_n = cumulative_windows_ge_n / sum(n_windows)
  )

threshold_summary <- tibble(
  minimum_n = c(4L, 5L, 6L, 8L, 10L),
  n_windows = map_int(
    minimum_n,
    ~ sum(density_manifest$n_footprints >= .x)
  ),
  n_assigned_footprints = map_int(
    minimum_n,
    ~ sum(
      density_manifest$n_footprints[
        density_manifest$n_footprints >= .x
      ]
    )
  )
) %>%
  mutate(
    pct_windows = n_windows / nrow(density_manifest),
    pct_assigned_footprints =
      n_assigned_footprints / sum(density_manifest$n_footprints)
  )

land_cover_summary <- density_manifest %>%
  group_by(majority_land_cover, majority_land_cover_class) %>%
  summarise(
    n_windows = n(),
    median_footprints = median(n_footprints),
    total_footprints = sum(n_footprints),
    median_land_cover_majority_fraction =
      median(land_cover_majority_fraction),
    mean_aggregate_sif = mean(aggregated_target_modis_sif),
    min_aggregate_sif = min(aggregated_target_modis_sif),
    median_aggregate_sif = median(aggregated_target_modis_sif),
    max_aggregate_sif = max(aggregated_target_modis_sif),
    .groups = "drop"
  ) %>%
  arrange(majority_land_cover)

write_csv(
  count_distribution,
  file.path(output_dir, "density_cluster_4000m_count_distribution.csv")
)
write_csv(
  threshold_summary,
  file.path(output_dir, "density_cluster_4000m_threshold_summary.csv")
)
write_csv(
  land_cover_summary,
  file.path(output_dir, "density_cluster_4000m_land_cover_summary.csv")
)

chip_polygons <- density_manifest %>%
  mutate(
    geometry = pmap(
      list(cell_xmin, cell_ymin, cell_xmax, cell_ymax),
      make_square_polygon
    )
  ) %>%
  st_as_sf(crs = sentinel_crs) %>%
  st_make_valid()

saveRDS(
  chip_polygons,
  file.path(output_dir, "density_cluster_4000m_chip_polygons.rds")
)

# A compact Leaflet preview is written for visual inspection. The plotted SIF
# polygons are only the majority-class footprints that contribute to each target.
set.seed(leaflet_seed)
preview_groups <- density_manifest %>%
  slice_sample(n = min(leaflet_sample_groups, nrow(density_manifest)))

if (nrow(preview_groups) > 0) {
  preview_ids <- preview_groups$aggregation_id

  preview_chips <- chip_polygons %>%
    filter(aggregation_id %in% preview_ids) %>%
    st_transform(4326) %>%
    mutate(
      popup = paste0(
        "<b>", aggregation_id, "</b>",
        "<br>Date: ", Delta_Date,
        "<br>Tile: ", mgrs_tile_t,
        "<br>Mode: ", measurement_mode,
        "<br>Majority class: ", majority_land_cover_class,
        " (", majority_land_cover, ")",
        "<br>Assigned footprints: ", n_footprints,
        "<br>All classes in window: ", window_n_all_land_covers,
        "<br>Land-cover fraction: ",
        round(land_cover_majority_fraction, 3),
        "<br>Mean SIF: ", round(aggregated_target_modis_sif, 4),
        "<br>Median SIF: ", round(median_target_modis_sif, 4),
        "<br>Min SIF: ", round(min_target_modis_sif, 4),
        "<br>Max SIF: ", round(max_target_modis_sif, 4)
      )
    )

  preview_footprints <- sif_sf %>%
    filter(
      sif_row_id %in% density_assignments$sif_row_id[
        density_assignments$aggregation_id %in% preview_ids
      ]
    ) %>%
    inner_join(
      density_assignments %>%
        filter(aggregation_id %in% preview_ids) %>%
        select(sif_row_id, aggregation_id),
      by = "sif_row_id"
    ) %>%
    st_transform(4326) %>%
    mutate(
      popup = paste0(
        "<b>SIF row ", sif_row_id, "</b>",
        "<br>Window: ", aggregation_id,
        "<br>Date: ", Delta_Date,
        "<br>SIF: ", round(target_modis_sif, 4),
        "<br>Land cover: ", land_cover_class,
        " (", land_cover, ")",
        "<br>BKR: ", BKR10_ID, " - ", BKR_NAME
      )
    )

  footprint_palette <- colorNumeric(
    palette = "viridis",
    domain = preview_footprints$target_modis_sif,
    na.color = "#808080"
  )

  preview_map <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(
      data = preview_chips,
      color = "#2166AC",
      weight = 2,
      fillColor = "#F4A582",
      fillOpacity = 0.18,
      popup = ~popup,
      group = "4 km density windows"
    ) %>%
    addPolygons(
      data = preview_footprints,
      color = "#303030",
      weight = 1,
      fillColor = ~footprint_palette(target_modis_sif),
      fillOpacity = 0.65,
      popup = ~popup,
      group = "Contributing SIF footprints"
    ) %>%
    addLayersControl(
      overlayGroups = c(
        "4 km density windows",
        "Contributing SIF footprints"
      ),
      options = layersControlOptions(collapsed = FALSE)
    ) %>%
    addLegend(
      pal = footprint_palette,
      values = preview_footprints$target_modis_sif,
      title = "target_modis_sif",
      opacity = 0.8
    )

  preview_bounds <- st_bbox(preview_chips)
  preview_map <- preview_map %>%
    fitBounds(
      preview_bounds[["xmin"]],
      preview_bounds[["ymin"]],
      preview_bounds[["xmax"]],
      preview_bounds[["ymax"]]
    )

  saveWidget(
    preview_map,
    file.path(output_dir, "density_cluster_4000m_leaflet_sample.html"),
    selfcontained = TRUE
  )
}

message("Done.")
message("Accepted 4 km windows: ", nrow(density_manifest))
message("Assigned majority-class footprints: ", nrow(density_assignments))
message("Excluded or unassigned footprints: ", nrow(density_exclusions))
message("Manifest: ", manifest_path)
message("Assignments: ", assignments_path)
message(
  "Each manifest row contains mean, median, minimum and maximum SIF, ",
  "footprint count, land-cover composition and BKR composition."
)



