library(tidyverse)
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
# represented by an explicit 4,000 m x 4,000 m square aligned to a 20 m MGRS
# grid. Only footprints belonging to the most frequent land-cover
# class inside that square contribute to the aggregate target. All footprints
# encountered by an accepted square are consumed, which prevents one source
# sounding from being used in multiple aggregate targets.

sf::sf_use_s2(FALSE)

input_csv <- paste0(
  "data/main_sif_data/",
  "redtiles_sif_BKR_landcover.csv"
)
mgrs_geometry_path <- "data/mgrs_de.rds"
output_dir <- "data/density_aggregation/sentinel2_spatial_aggregation_density_4000m_landcover_redtiles"

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

leaflet_sample_groups <- 100L
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

build_tile_table <- function(mgrs_path, requested_tiles) {
  if (!file.exists(mgrs_path)) {
    stop("MGRS geometry file not found: ", mgrs_path)
  }

  mgrs_all <- readRDS(mgrs_path)
  if (!inherits(mgrs_all, "sf")) {
    stop("The MGRS geometry RDS must contain an sf object.")
  }
  if (!"mgrs_tile" %in% names(mgrs_all)) {
    stop("The MGRS geometry object has no mgrs_tile column.")
  }

  requested_tiles <- sort(unique(na.omit(normalize_mgrs_tile(requested_tiles))))
  if (length(requested_tiles) == 0) {
    stop("No valid assigned MGRS tile IDs were found in the SIF data.")
  }

  mgrs_selected <- mgrs_all %>%
    mutate(
      mgrs_tile_original = stringr::str_remove(
        normalize_mgrs_tile(mgrs_tile),
        "^T"
      ),
      mgrs_tile_t = normalize_mgrs_tile(mgrs_tile)
    ) %>%
    filter(mgrs_tile_t %in% requested_tiles)

  missing_tiles <- setdiff(requested_tiles, mgrs_selected$mgrs_tile_t)
  if (length(missing_tiles) > 0) {
    stop(
      "Assigned tiles missing from the MGRS geometry object: ",
      paste(missing_tiles, collapse = ", ")
    )
  }

  duplicate_tiles <- mgrs_selected %>%
    st_drop_geometry() %>%
    count(mgrs_tile_t) %>%
    filter(n != 1L)
  if (nrow(duplicate_tiles) > 0) {
    stop("At least one requested MGRS tile has duplicate geometries.")
  }

  if (!"crs_utm" %in% names(mgrs_selected)) {
    stop("The MGRS geometry object has no crs_utm column.")
  }

  tile_crs_values <- mgrs_selected %>%
    st_drop_geometry() %>%
    transmute(crs_utm = stringr::str_trim(as.character(crs_utm))) %>%
    filter(!is.na(crs_utm), crs_utm != "") %>%
    distinct() %>%
    pull(crs_utm)

  if (length(tile_crs_values) != 1L) {
    stop(
      "All assigned tiles must use one projected CRS for joint clustering; ",
      "found: ", paste(tile_crs_values, collapse = ", ")
    )
  }

  analysis_crs <- st_crs(tile_crs_values[[1]])
  if (is.na(analysis_crs)) {
    stop("Could not interpret the assigned tiles' crs_utm value.")
  }

  tile_boundaries <- mgrs_selected %>%
    st_make_valid() %>%
    st_transform(analysis_crs) %>%
    select(
      mgrs_tile_original,
      mgrs_tile_t,
      any_of(c("grid_zone", "square_100km", "zone", "crs_utm")),
      geometry
    ) %>%
    arrange(mgrs_tile_t)

  tile_table <- map_dfr(seq_len(nrow(tile_boundaries)), function(tile_i) {
    tile_bbox <- st_bbox(tile_boundaries[tile_i, ])

    tibble(
      mgrs_tile_original = tile_boundaries$mgrs_tile_original[[tile_i]],
      mgrs_tile_t = tile_boundaries$mgrs_tile_t[[tile_i]],
      window_crs = tile_crs_values[[1]],
      tile_xmin = unname(tile_bbox[["xmin"]]),
      tile_xmax = unname(tile_bbox[["xmax"]]),
      tile_ymin = unname(tile_bbox[["ymin"]]),
      tile_ymax = unname(tile_bbox[["ymax"]]),
      grid_xmin = ceiling((unname(tile_bbox[["xmin"]]) - 1e-6) /
        sentinel_pixel_size_m) * sentinel_pixel_size_m,
      grid_xmax = floor((unname(tile_bbox[["xmax"]]) + 1e-6) /
        sentinel_pixel_size_m) * sentinel_pixel_size_m,
      grid_ymin = ceiling((unname(tile_bbox[["ymin"]]) - 1e-6) /
        sentinel_pixel_size_m) * sentinel_pixel_size_m,
      grid_ymax = floor((unname(tile_bbox[["ymax"]]) + 1e-6) /
        sentinel_pixel_size_m) * sentinel_pixel_size_m
    )
  })

  if (any(tile_table$grid_xmax - tile_table$grid_xmin < window_size_m) ||
      any(tile_table$grid_ymax - tile_table$grid_ymin < window_size_m)) {
    stop("At least one assigned MGRS tile is smaller than the 4 km window.")
  }

  list(
    table = tile_table,
    boundaries = tile_boundaries,
    crs = analysis_crs
  )
}

# Snap a candidate square to the assigned MGRS tile's 20 m grid and, where
# needed, shift it inward so the complete 4 km window stays inside the tile.
snap_window_to_tile <- function(center_x, center_y, tile_row) {
  xmin_raw <- center_x - half_window_m
  ymax_raw <- center_y + half_window_m

  xmin <- tile_row$grid_xmin +
    round((xmin_raw - tile_row$grid_xmin) / sentinel_pixel_size_m) *
      sentinel_pixel_size_m
  ymax <- tile_row$grid_ymax -
    round((tile_row$grid_ymax - ymax_raw) / sentinel_pixel_size_m) *
      sentinel_pixel_size_m

  xmin <- min(
    max(xmin, tile_row$grid_xmin),
    tile_row$grid_xmax - window_size_m
  )
  ymax <- min(
    max(ymax, tile_row$grid_ymin + window_size_m),
    tile_row$grid_ymax
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

    # If the first minimum-sized PCA-ordered group cannot share one 4 km
    # support, remove the leading point and retry without rerunning DBSCAN.
    if (fit$n < minimum_footprints) {
      excluded[[length(excluded) + 1L]] <- tibble(
        sif_row_id = remaining$sif_row_id[[1]],
        exclusion_reason = "cannot_form_2km_track_chunk",
        related_aggregation_id = NA_character_
      )
      remaining <- remaining[-1, , drop = FALSE]
      next
    }

    # Avoid leaving a final undersized remainder where a smaller first chunk
    # can leave enough points for the following window.
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
        exclusion_reason = "fewer_than_min_n_majority_land_cover_footprints",
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
      exclusion_reason = "density_noise_or_fewer_than_minimum_footprints",
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
        exclusion_reason = "stratum_fewer_than_minimum_footprints",
        related_aggregation_id = NA_character_
      )
    ))
  }

  tile_row <- stratum_tbl %>%
    slice(1) %>%
    select(
      tile_xmin, tile_xmax, tile_ymin, tile_ymax,
      grid_xmin, grid_xmax, grid_ymin, grid_ymax
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
  "final_check_modis_sif", "Quality_Flag", "mgrs_tile",
  "BKR10_ID", "BKR_NAME", "land_cover",
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
    target_modis_sif = as.numeric(target_modis_sif),
    mgrs_tile_original = stringr::str_remove(
      normalize_mgrs_tile(mgrs_tile),
      "^T"
    ),
    mgrs_tile_t = normalize_mgrs_tile(mgrs_tile),
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
    is.finite(target_modis_sif),
    stringr::str_to_lower(final_check_modis_sif) == "accept",
    stringr::str_detect(mgrs_tile_t, "^T[0-9]{2}[A-Z]{3}$"),
    !is.na(land_cover),
    !is.na(land_cover_class),
    land_cover_class != "",
    if_all(all_of(corner_cols), ~ is.finite(.x))
  )

message("Reading assigned MGRS tile geometries: ", mgrs_geometry_path)
tile_result <- build_tile_table(
  mgrs_geometry_path,
  df_prepared$mgrs_tile_t
)
tile_info <- tile_result$table
tile_boundaries <- tile_result$boundaries
analysis_crs <- tile_result$crs

write_csv(
  tile_info,
  file.path(output_dir, "density_4000m_mgrs_geometry_reference.csv")
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
  st_transform(analysis_crs)

sif_centroids <- suppressWarnings(st_centroid(sif_sf))
centroid_xy <- st_coordinates(sif_centroids)

sif_sf <- sif_sf %>%
  mutate(
    centroid_x = centroid_xy[, "X"],
    centroid_y = centroid_xy[, "Y"],
    input_tile_available = mgrs_tile_t %in% tile_info$mgrs_tile_t,
    base_keep_for_grouping = input_tile_available
  )

# -----------------------------------------------------------------------------
# Remove edge footprints whose polygons are mostly outside their assigned tile.
# The tile ID already assigned to each SIF row is trusted; the matching polygon
# from mgrs_de.rds supplies only its boundary. A footprint is retained when at
# least 80% of its polygon area lies inside that boundary. This is calculated
# before density clustering.
# -----------------------------------------------------------------------------
tile_boundaries <- tile_boundaries %>%
  select(mgrs_tile_t, geometry)

edge_candidates <- sif_sf %>%
  filter(base_keep_for_grouping)

if (nrow(edge_candidates) == 0) {
  stop("No SIF polygons remain before the assigned-tile edge filter.")
}

tile_overlap_fractions <- edge_candidates %>%
  group_by(mgrs_tile_t) %>%
  group_split(.keep = TRUE) %>%
  map_dfr(function(tile_sif) {
    tile_id <- first(tile_sif$mgrs_tile_t)
    tile_boundary <- tile_boundaries %>%
      filter(mgrs_tile_t == tile_id) %>%
      st_geometry()

    full_areas <- tibble(
      sif_row_id = tile_sif$sif_row_id,
      footprint_area_m2 = as.numeric(st_area(tile_sif))
    )

    intersections <- suppressWarnings(
      st_intersection(
        tile_sif %>% select(sif_row_id),
        tile_boundary
      )
    )

    inside_areas <- if (nrow(intersections) == 0) {
      tibble(
        sif_row_id = integer(),
        footprint_inside_tile_m2 = numeric()
      )
    } else {
      intersections %>%
        mutate(
          footprint_inside_tile_m2 = as.numeric(st_area(geometry))
        ) %>%
        st_drop_geometry() %>%
        group_by(sif_row_id) %>%
        summarise(
          footprint_inside_tile_m2 = sum(footprint_inside_tile_m2),
          .groups = "drop"
        )
    }

    full_areas %>%
      left_join(inside_areas, by = "sif_row_id") %>%
      mutate(
        footprint_inside_tile_m2 = coalesce(
          footprint_inside_tile_m2,
          0
        ),
        tile_inside_fraction = case_when(
          footprint_area_m2 > 0 ~ pmin(
            1,
            footprint_inside_tile_m2 / footprint_area_m2
          ),
          TRUE ~ NA_real_
        ),
        tile_outside_fraction = 1 - tile_inside_fraction
      )
  })

sif_sf <- sif_sf %>%
  left_join(tile_overlap_fractions, by = "sif_row_id") %>%
  mutate(
    passes_tile_edge_filter =
      base_keep_for_grouping &
      is.finite(tile_inside_fraction) &
      tile_inside_fraction >= 0.80,
    keep_for_grouping = passes_tile_edge_filter
  )

sif_df <- sif_sf %>%
  filter(keep_for_grouping)

message(
  "SIF polygons before assigned-tile edge filter: ",
  format(nrow(edge_candidates), big.mark = ",")
)
message(
  "SIF polygons removed (<80% inside assigned tile): ",
  format(nrow(edge_candidates) - nrow(sif_df), big.mark = ",")
)
message(
  "SIF polygons retained for density clustering (nrow(sif_df)): ",
  format(nrow(sif_df), big.mark = ",")
)

row_audit <- sif_sf %>%
  st_drop_geometry() %>%
  count(
    mgrs_tile_t,
    input_tile_available,
    passes_tile_edge_filter,
    keep_for_grouping,
    name = "n_sif_rows"
  ) %>%
  arrange(mgrs_tile_t)

write_csv(
  row_audit,
  file.path(output_dir, "density_4000m_input_row_audit.csv")
)

study_tbl <- sif_df %>%
  st_drop_geometry() %>%
  inner_join(
    tile_info %>%
      select(
        mgrs_tile_t,
        window_crs,
        tile_xmin, tile_xmax, tile_ymin, tile_ymax,
        grid_xmin, grid_xmax, grid_ymin, grid_ymax
      ),
    by = "mgrs_tile_t"
  )

if (nrow(study_tbl) == 0) {
  stop("No rows remain after the assigned-tile edge filter.")
}

# A stratum represents one OCO-2 track subset within one assigned MGRS tile,
# date and measurement mode. Sentinel products are acquired later for the
# accepted windows and therefore play no role in this grouping stage.
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
    n_land_covers = n_distinct(land_cover),
    full_window_size =
      abs(first(cell_xmax) - first(cell_xmin) - window_size_m) < 1e-8 &&
      abs(first(cell_ymax) - first(cell_ymin) - window_size_m) < 1e-8,
    window_inside_tile_bounds =
      first(cell_xmin) >= first(tile_xmin) &&
      first(cell_xmax) <= first(tile_xmax) &&
      first(cell_ymin) >= first(tile_ymin) &&
      first(cell_ymax) <= first(tile_ymax),
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
    any(!group_integrity$full_window_size) ||
    any(!group_integrity$window_inside_tile_bounds) ||
    any(!group_integrity$centroids_inside_window)) {
  stop("At least one accepted density group failed its integrity checks.")
}

#---------------------
group_integrity %>%
  summarise(
    bad_dates = sum(n_dates != 1),
    bad_modes = sum(n_modes != 1),
    bad_tiles = sum(n_tiles != 1),
    bad_land_cover = sum(n_land_covers != 1),
    bad_window_size = sum(!full_window_size),
    bad_tile_bounds = sum(!window_inside_tile_bounds),
    bad_centroids = sum(!centroids_inside_window)
  )

bad_ids <- group_integrity %>%
  filter(!window_inside_tile_bounds) %>%
  pull(aggregation_id)

density_assignments_full %>%
  filter(aggregation_id %in% bad_ids) %>%
  group_by(aggregation_id, mgrs_tile_t) %>%
  summarise(
    left_excess_m   = max(first(tile_xmin) - first(cell_xmin), 0),
    right_excess_m  = max(first(cell_xmax) - first(tile_xmax), 0),
    bottom_excess_m = max(first(tile_ymin) - first(cell_ymin), 0),
    top_excess_m    = max(first(cell_ymax) - first(tile_ymax), 0),
    .groups = "drop"
  )


# Verify against the actual MGRS polygons as well as their numeric bounds.
# Windows failing this check are not cropped; the script stops instead.
window_geometry_checks <- density_assignments_full %>%
  distinct(
    aggregation_id, mgrs_tile_t,
    cell_xmin, cell_ymin, cell_xmax, cell_ymax
  ) %>%
  mutate(
    geometry = pmap(
      list(cell_xmin, cell_ymin, cell_xmax, cell_ymax),
      make_square_polygon
    )
  ) %>%
  st_as_sf(crs = analysis_crs) %>%
  group_by(mgrs_tile_t) %>%
  group_split(.keep = TRUE) %>%
  map_dfr(function(tile_windows) {
    tile_id <- first(tile_windows$mgrs_tile_t)
    tile_boundary <- tile_boundaries %>%
      filter(mgrs_tile_t == tile_id)

    tibble(
      aggregation_id = tile_windows$aggregation_id,
      window_inside_mgrs_geometry = lengths(
        st_covered_by(tile_windows, tile_boundary)
      ) == 1L
    )
  })

if (any(!window_geometry_checks$window_inside_mgrs_geometry)) {
  stop(
    "At least one complete 4 km window extends beyond its assigned MGRS ",
    "geometry. No windows were cropped."
  )
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
    window_crs,
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
    mean_tile_inside_fraction = safe_mean(tile_inside_fraction),
    min_tile_inside_fraction = safe_min(tile_inside_fraction),
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
    window_crs,
    target_modis_sif,
    final_check_modis_sif,
    Quality_Flag,
    BKR10_ID,
    BKR_NAME,
    land_cover,
    land_cover_class,
    tile_inside_fraction,
    tile_outside_fraction,
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
  st_as_sf(crs = analysis_crs) %>%
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

#-------------------------------------------------------------------------------

# summary(density_manifest$aggregated_target_modis_sif)
# 
# noise_df <- tibble(
#   n = 1:20,
#   noise_remaining = 1 / sqrt(n),
#   noise_reduction = 1 - noise_remaining
# ) %>%
#   pivot_longer(
#     cols = c(noise_remaining, noise_reduction),
#     names_to = "metric",
#     values_to = "value"
#   ) %>%
#   mutate(
#     metric = recode(
#       metric,
#       noise_remaining = "Noise remaining: 1/sqrt(n)",
#       noise_reduction = "Noise reduction: 1 - 1/sqrt(n)"
#     )
#   )
# 
# ggplot(noise_df, aes(x = n, y = value, color = metric)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.8) +
#   scale_x_continuous(breaks = 1:20) +
#   scale_y_continuous(
#     labels = scales::label_percent(),
#     limits = c(0, 1)
#   ) +
#   labs(
#     x = "Number of averaged soundings (n)",
#     y = "Percentage",
#     color = NULL
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )



