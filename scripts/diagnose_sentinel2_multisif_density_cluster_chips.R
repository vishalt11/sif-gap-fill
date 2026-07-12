library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(dbscan)
library(leaflet)
library(htmlwidgets)

# Sentinel-2 multi-footprint chip manifest/assignment builder.
#
# Goal:
#   Build same-date, same-MGRS-tile chip groups for later Sentinel-2 band-chip
#   extraction. This mirrors the MODIS multi-SIF DBSCAN diagnostic, but uses the
#   Sentinel-2 raster CRS/grid and 6 km x 6 km chips.
#
# Important:
#   The output script does not create image chips. It only writes:
#     1. one chip manifest row per 6 km chip
#     2. one assignment row per SIF footprint assigned to that chip
#
# Later Python/R chip extraction can use chip_xmin/ymin/xmax/ymax and the source
# Sentinel tile path to crop/rasterize bands and footprint masks.

sf::sf_use_s2(FALSE)

input_rds <- "data/9tiles_2_7_inrange_M01_cnn.rds"
mgrs_tif_dir <- "data/temp_data/mgrs_tifs"
output_dir <- "data/sentinel2_multisif_chip_diagnostics"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sentinel-2 is 10 m here. A 6 km chip corresponds to about 600 x 600 pixels.
chip_size_m <- 6000

# DBSCAN on footprint centroids, separately for each MGRS tile and date.
# min_sif_per_chip = 4 means keep clusters/chunks with at least 4 supervision
# footprints. max_sif_per_chip avoids making very long overpass-track chips.
dbscan_eps_m <- 1800
min_sif_per_chip <- 4
max_sif_per_chip <- 10



corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

target_cols <- c("Daily_SIF_757nm", "Daily_SIF_771nm", "target_modis_sif")
final_check_cols <- c("final_check_757", "final_check_771", "final_check_modis_sif")

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

normalize_mgrs_tile <- function(x) {
  tile <- as.character(x)
  tile <- str_trim(tile)
  tile <- str_remove(tile, "^T")
  if_else(is.na(tile) | tile == "", NA_character_, paste0("T", tile))
}

extract_mgrs_tile_from_path <- function(path) {
  str_extract(basename(path), "T\\d{2}[A-Z]{3}")
}

make_balanced_chunks <- function(n, min_n, max_n) {
  if (n < min_n) {
    return(list())
  }

  if (n <= max_n) {
    return(list(seq_len(n)))
  }

  n_chunks <- ceiling(n / max_n)

  while (floor(n / n_chunks) < min_n && n_chunks > 1) {
    n_chunks <- n_chunks - 1
  }

  chunk_sizes <- rep(floor(n / n_chunks), n_chunks)
  remainder <- n %% n_chunks

  if (remainder > 0) {
    chunk_sizes[seq_len(remainder)] <- chunk_sizes[seq_len(remainder)] + 1
  }

  chunk_ends <- cumsum(chunk_sizes)
  chunk_starts <- c(1, head(chunk_ends, -1) + 1)

  map2(chunk_starts, chunk_ends, seq)
}

add_pca_track_score <- function(cluster_tbl) {
  coords <- as.matrix(cluster_tbl[, c("centroid_x", "centroid_y")])

  if (nrow(coords) < 2) {
    return(cluster_tbl %>% mutate(track_score = centroid_x))
  }

  track_score <- tryCatch(
    {
      stats::prcomp(coords, center = TRUE, scale. = FALSE)$x[, 1]
    },
    error = function(e) {
      coords[, 1]
    }
  )

  cluster_tbl %>%
    mutate(track_score = track_score)
}

snap_chip_bounds_to_raster <- function(center_x, center_y, raster_info) {
  half_size <- chip_size_m / 2

  raw_xmin <- center_x - half_size
  raw_ymax <- center_y + half_size

  # Align the chip to the Sentinel raster grid so later extraction can use an
  # exact pixel window. We anchor x to xmin and y to ymax because north-up
  # rasters commonly index rows downward from ymax.
  col_off <- floor((raw_xmin - raster_info$xmin) / raster_info$xres)
  row_off <- floor((raster_info$ymax - raw_ymax) / raster_info$yres)

  chip_cols <- round(chip_size_m / raster_info$xres)
  chip_rows <- round(chip_size_m / raster_info$yres)

  xmin <- raster_info$xmin + col_off * raster_info$xres
  xmax <- xmin + chip_cols * raster_info$xres
  ymax <- raster_info$ymax - row_off * raster_info$yres
  ymin <- ymax - chip_rows * raster_info$yres

  tibble(
    chip_xmin = xmin,
    chip_ymin = ymin,
    chip_xmax = xmax,
    chip_ymax = ymax,
    chip_cols = chip_cols,
    chip_rows = chip_rows,
    col_off = col_off,
    row_off = row_off
  )
}

chip_inside_tile_status <- function(bounds_tbl, raster_info) {
  inside <- bounds_tbl$chip_xmin >= raster_info$xmin &&
    bounds_tbl$chip_xmax <= raster_info$xmax &&
    bounds_tbl$chip_ymin >= raster_info$ymin &&
    bounds_tbl$chip_ymax <= raster_info$ymax

  if (inside) {
    return("inside")
  }

  "edge_flag"
}

build_tile_table <- function(tif_dir) {
  tif_paths <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

  if (length(tif_paths) == 0) {
    stop("No .tif files found in: ", tif_dir)
  }

  tile_tbl <- tibble(
    mgrs_tile_t = extract_mgrs_tile_from_path(tif_paths),
    sentinel_tif_path = tif_paths
  ) %>%
    filter(!is.na(mgrs_tile_t)) %>%
    distinct(mgrs_tile_t, .keep_all = TRUE) %>%
    arrange(mgrs_tile_t)

  if (nrow(tile_tbl) == 0) {
    stop("Could not extract T-prefixed MGRS tile codes from tif names.")
  }

  tile_tbl %>%
    mutate(
      raster = map(sentinel_tif_path, terra::rast),
      raster_crs = map_chr(raster, terra::crs),
      raster_ext = map(raster, terra::ext),
      raster_res = map(raster, terra::res),
      xmin = map_dbl(raster_ext, ~ .x$xmin),
      xmax = map_dbl(raster_ext, ~ .x$xmax),
      ymin = map_dbl(raster_ext, ~ .x$ymin),
      ymax = map_dbl(raster_ext, ~ .x$ymax),
      xres = map_dbl(raster_res, 1),
      yres = map_dbl(raster_res, 2),
      ncol = map_int(raster, terra::ncol),
      nrow = map_int(raster, terra::nrow),
      geometry = pmap(
        list(xmin, ymin, xmax, ymax),
        make_square_polygon
      )
    ) %>%
    select(-raster, -raster_ext, -raster_res) %>%
    st_as_sf(crs = terra::crs(terra::rast(tile_tbl$sentinel_tif_path[[1]])))
}

split_dbscan_cluster <- function(cluster_tbl) {
  n_cluster <- nrow(cluster_tbl)

  if (n_cluster < min_sif_per_chip) {
    return(tibble())
  }

  cluster_tbl <- cluster_tbl %>%
    add_pca_track_score() %>%
    arrange(track_score)

  chunks <- make_balanced_chunks(
    n = nrow(cluster_tbl),
    min_n = min_sif_per_chip,
    max_n = max_sif_per_chip
  )

  if (length(chunks) == 0) {
    return(tibble())
  }

  date_string <- format(cluster_tbl$Delta_Date[[1]], "%Y%m%d")
  mgrs_tile_t <- cluster_tbl$mgrs_tile_t[[1]]
  dbscan_cluster <- cluster_tbl$dbscan_cluster[[1]]

  imap_dfr(chunks, function(chunk_rows, chunk_i) {
    chunk_tbl <- cluster_tbl[chunk_rows, ]
    chip_id <- paste0(
      "s2_6km_",
      mgrs_tile_t,
      "_",
      date_string,
      "_db",
      dbscan_cluster,
      "_chunk",
      chunk_i
    )

    chunk_tbl %>%
      mutate(
        chip_size_m = chip_size_m,
        chip_id = chip_id,
        dbscan_cluster_n = n_cluster,
        chunk_number = chunk_i,
        chunk_n_sif = nrow(chunk_tbl)
      )
  })
}

cluster_one_tile_date <- function(group_tbl) {
  if (nrow(group_tbl) < min_sif_per_chip) {
    return(tibble())
  }

  coords <- as.matrix(group_tbl[, c("centroid_x", "centroid_y")])

  dbscan_result <- dbscan::dbscan(
    coords,
    eps = dbscan_eps_m,
    minPts = min_sif_per_chip
  )

  clustered_tbl <- group_tbl %>%
    mutate(dbscan_cluster = dbscan_result$cluster) %>%
    filter(dbscan_cluster > 0)

  if (nrow(clustered_tbl) < min_sif_per_chip) {
    return(tibble())
  }

  clustered_tbl %>%
    group_by(dbscan_cluster) %>%
    group_split() %>%
    map_dfr(split_dbscan_cluster)
}

message("Reading Sentinel tile tifs from: ", mgrs_tif_dir)
tile_sf <- build_tile_table(mgrs_tif_dir)

sentinel_crs <- st_crs(tile_sf)
sentinel_crs_wkt <- sentinel_crs$wkt

tile_info <- tile_sf %>%
  st_drop_geometry() %>%
  select(
    mgrs_tile_t,
    sentinel_tif_path,
    raster_crs,
    xmin,
    xmax,
    ymin,
    ymax,
    xres,
    yres,
    ncol,
    nrow
  )

message("Sentinel tiles found: ", paste(tile_info$mgrs_tile_t, collapse = ", "))

message("Reading SIF file: ", input_rds)
df_raw <- readRDS(input_rds)

df <- if (inherits(df_raw, "sf")) {
  st_drop_geometry(df_raw)
} else {
  as_tibble(df_raw)
}

required_cols <- c("Delta_Date", "mgrs_tile", corner_cols)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!"target_modis_sif" %in% names(df) &&
    all(c("Daily_SIF_757nm", "Daily_SIF_771nm") %in% names(df))) {
  df <- df %>%
    mutate(target_modis_sif = (Daily_SIF_757nm + 1.5 * Daily_SIF_771nm) / 2)
}

message("Rows read: ", nrow(df))

sif_sf <- df %>%
  mutate(
    sif_row_id = row_number(),
    Delta_Date = as.Date(Delta_Date),
    sif_year = lubridate::year(Delta_Date),
    sif_doy = lubridate::yday(Delta_Date),
    mgrs_tile = str_remove(as.character(mgrs_tile), "^T"),
    mgrs_tile_t = normalize_mgrs_tile(mgrs_tile),
    across(all_of(corner_cols), as.numeric)
  ) %>%
  filter(
    !is.na(Delta_Date),
    !is.na(mgrs_tile_t),
    if_all(all_of(corner_cols), ~ !is.na(.x))
  ) %>%
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
  st_make_valid()

message("Valid SIF polygons: ", nrow(sif_sf))

missing_tiles <- setdiff(unique(sif_sf$mgrs_tile_t), tile_info$mgrs_tile_t)
if (length(missing_tiles) > 0) {
  warning(
    "No Sentinel tif found for these SIF MGRS tiles; rows will be dropped: ",
    paste(missing_tiles, collapse = ", ")
  )
}

sif_projected <- sif_sf %>%
  filter(mgrs_tile_t %in% tile_info$mgrs_tile_t) %>%
  st_transform(sentinel_crs)

sif_centroids <- suppressWarnings(st_centroid(sif_projected))
centroid_xy <- st_coordinates(sif_centroids)

sif_for_clustering <- sif_projected %>%
  st_drop_geometry() %>%
  mutate(
    centroid_x = centroid_xy[, 1],
    centroid_y = centroid_xy[, 2]
  ) %>%
  left_join(tile_info, by = "mgrs_tile_t")

message("Rows available for clustering: ", nrow(sif_for_clustering))
message("Running same-tile, same-date DBSCAN clustering...")

density_assignments <- sif_for_clustering %>%
  group_by(mgrs_tile_t, Delta_Date) %>%
  group_split() %>%
  map_dfr(cluster_one_tile_date)

if (nrow(density_assignments) == 0) {
  stop("No Sentinel-2 density-cluster chips were found. Try increasing dbscan_eps_m or lowering min_sif_per_chip.")
}

dropped_counts <- density_assignments %>%
  count(chunk_n_sif, name = "n_sif_rows") %>%
  arrange(chunk_n_sif)

write_csv(
  dropped_counts,
  file.path(output_dir, "sentinel2_density_cluster_chunk_counts.csv")
)

message("Assigned SIF rows before chunk-count filter: ", nrow(density_assignments))

density_assignments <- density_assignments %>%
  filter(chunk_n_sif >= min_sif_per_chip, chunk_n_sif <= max_sif_per_chip)

if (nrow(density_assignments) == 0) {
  stop("No Sentinel-2 chips remain after applying min/max SIF count filter.")
}

message("Assigned SIF rows after chunk-count filter: ", nrow(density_assignments))
message("Unique Sentinel-2 chips: ", n_distinct(density_assignments$chip_id))

write_csv(
  density_assignments,
  file.path(output_dir, "sentinel2_density_cluster_sif_assignments.csv")
)

has_target_modis_sif <- "target_modis_sif" %in% names(density_assignments)
has_sif_757 <- "Daily_SIF_757nm" %in% names(density_assignments)
has_sif_771 <- "Daily_SIF_771nm" %in% names(density_assignments)
has_state <- "state" %in% names(density_assignments)
has_hzs <- "hzs" %in% names(density_assignments)
has_source_file <- "source_file" %in% names(density_assignments)

chip_summary_base <- density_assignments %>%
  group_by(chip_id, mgrs_tile_t, mgrs_tile, Delta_Date) %>%
  summarise(
    sif_year = first(sif_year),
    sif_doy = first(sif_doy),
    chip_size_m = first(chip_size_m),
    chip_center_x = mean(centroid_x),
    chip_center_y = mean(centroid_y),
    n_sif = n(),
    dbscan_cluster_n = first(dbscan_cluster_n),
    chunk_number = first(chunk_number),
    chunk_n_sif = first(chunk_n_sif),
    n_state = if (has_state) n_distinct(state, na.rm = TRUE) else NA_integer_,
    states = if (has_state) paste(sort(unique(na.omit(state))), collapse = ";") else NA_character_,
    n_hzs = if (has_hzs) n_distinct(hzs, na.rm = TRUE) else NA_integer_,
    hzs_values = if (has_hzs) paste(sort(unique(na.omit(hzs))), collapse = ";") else NA_character_,
    n_source_file = if (has_source_file) n_distinct(source_file, na.rm = TRUE) else NA_integer_,
    mean_target_modis_sif = if (has_target_modis_sif) mean(target_modis_sif, na.rm = TRUE) else NA_real_,
    min_target_modis_sif = if (has_target_modis_sif) min(target_modis_sif, na.rm = TRUE) else NA_real_,
    max_target_modis_sif = if (has_target_modis_sif) max(target_modis_sif, na.rm = TRUE) else NA_real_,
    mean_Daily_SIF_757nm = if (has_sif_757) mean(Daily_SIF_757nm, na.rm = TRUE) else NA_real_,
    mean_Daily_SIF_771nm = if (has_sif_771) mean(Daily_SIF_771nm, na.rm = TRUE) else NA_real_,
    sif_row_ids = paste(sif_row_id, collapse = ","),
    .groups = "drop"
  ) %>%
  left_join(tile_info, by = "mgrs_tile_t")

chip_bounds <- chip_summary_base %>%
  pmap_dfr(function(...) {
    row <- list(...)
    raster_info <- list(
      xmin = row$xmin,
      xmax = row$xmax,
      ymin = row$ymin,
      ymax = row$ymax,
      xres = row$xres,
      yres = row$yres
    )

    bounds <- snap_chip_bounds_to_raster(
      center_x = row$chip_center_x,
      center_y = row$chip_center_y,
      raster_info = raster_info
    )

    bounds %>%
      mutate(
        chip_status = chip_inside_tile_status(bounds, raster_info),
        chip_area_m2 = (chip_xmax - chip_xmin) * (chip_ymax - chip_ymin)
      )
  })

chip_manifest <- bind_cols(chip_summary_base, chip_bounds) %>%
  mutate(
    chip_geometry = pmap(
      list(chip_xmin, chip_ymin, chip_xmax, chip_ymax),
      make_square_polygon
    )
  )

chip_manifest_sf <- st_as_sf(
  chip_manifest,
  sf_column_name = "chip_geometry",
  crs = sentinel_crs
)

chip_manifest_csv <- chip_manifest %>%
  select(
    chip_id,
    mgrs_tile,
    mgrs_tile_t,
    Delta_Date,
    sif_year,
    sif_doy,
    sentinel_tif_path,
    raster_crs,
    chip_size_m,
    chip_rows,
    chip_cols,
    row_off,
    col_off,
    chip_xmin,
    chip_ymin,
    chip_xmax,
    chip_ymax,
    chip_center_x,
    chip_center_y,
    chip_status,
    chip_area_m2,
    n_sif,
    dbscan_cluster_n,
    chunk_number,
    chunk_n_sif,
    n_state,
    states,
    n_hzs,
    hzs_values,
    n_source_file,
    mean_target_modis_sif,
    min_target_modis_sif,
    max_target_modis_sif,
    mean_Daily_SIF_757nm,
    mean_Daily_SIF_771nm,
    sif_row_ids
  )

write_csv(
  chip_manifest_csv,
  file.path(output_dir, "sentinel2_multi_sif_6km_chip_manifest.csv")
)

saveRDS(
  chip_manifest_sf,
  file.path(output_dir, "sentinel2_multi_sif_6km_chip_manifest.rds")
)

assignment_export_cols <- c(
  "chip_id", "sif_row_id", "Delta_Date", "sif_year", "sif_doy",
  "mgrs_tile", "mgrs_tile_t", "sentinel_tif_path",
  "chip_size_m", "chunk_number", "chunk_n_sif",
  "chip_status", "chip_rows", "chip_cols", "row_off", "col_off",
  "chip_xmin", "chip_ymin", "chip_xmax", "chip_ymax",
  "centroid_x", "centroid_y", "track_score",
  target_cols,
  final_check_cols,
  "Latitude", "Longitude","product_path",
  corner_cols,
  "state", "hzs", "Quality_Flag", "Metadata.MeasurementMode", "source_file"
)

assignment_chip_cols <- chip_manifest_csv %>%
  select(
    chip_id,
    chip_status,
    chip_rows,
    chip_cols,
    row_off,
    col_off,
    chip_xmin,
    chip_ymin,
    chip_xmax,
    chip_ymax
  )

density_assignments_for_export <- density_assignments %>%
  left_join(assignment_chip_cols, by = "chip_id")

assignment_export_cols <- intersect(assignment_export_cols, names(density_assignments_for_export))

python_chip_assignments <- density_assignments_for_export %>%
  select(all_of(assignment_export_cols)) %>%
  arrange(mgrs_tile_t, Delta_Date, chip_id, track_score, sif_row_id)

write_csv(
  python_chip_assignments,
  file.path(output_dir, "sentinel2_multi_sif_6km_chip_assignments.csv")
)

chip_count_distribution <- chip_manifest_csv %>%
  count(n_sif, chip_status, name = "n_chips") %>%
  group_by(chip_status) %>%
  mutate(
    pct_chips = n_chips / sum(n_chips),
    cum_chips_ge_n = rev(cumsum(rev(n_chips))),
    pct_chips_ge_n = cum_chips_ge_n / sum(n_chips)
  ) %>%
  ungroup() %>%
  arrange(chip_status, n_sif)

write_csv(
  chip_count_distribution,
  file.path(output_dir, "sentinel2_density_cluster_chip_count_distribution.csv")
)

tile_status_summary <- chip_manifest_csv %>%
  count(mgrs_tile_t, chip_status, name = "n_chips") %>%
  arrange(mgrs_tile_t, chip_status)

write_csv(
  tile_status_summary,
  file.path(output_dir, "sentinel2_chip_tile_status_summary.csv")
)

message("Wrote Sentinel-2 chip diagnostics to: ", output_dir)
message("  - sentinel2_multi_sif_6km_chip_manifest.csv")
message("  - sentinel2_multi_sif_6km_chip_assignments.csv")
message("  - sentinel2_multi_sif_6km_chip_manifest.rds")
message("  - sentinel2_density_cluster_chip_count_distribution.csv")
message("  - sentinel2_chip_tile_status_summary.csv")

message("Chip status summary:")
print(tile_status_summary)

message("Chip count distribution:")
print(chip_count_distribution, n = 100)

python_chip_assignments %>%
  group_by(chip_id) %>%
  summarise(
    n_dates = n_distinct(Delta_Date),
    n_tiles = n_distinct(mgrs_tile_t),
    n_sif_rows = n(),
    n_unique_sif_rows = n_distinct(sif_row_id),
    .groups = "drop"
  ) %>%
  summarise(
    max_dates = max(n_dates),
    max_tiles = max(n_tiles),
    any_duplicate_sif_inside_chip = any(n_sif_rows != n_unique_sif_rows)
  )

# ---------------------------------------------------------------------------
# Leaflet preview for one Sentinel tile/status

leaflet_tile <- "T32UPC"
leaflet_chip_status <- "edge_flag"

message(
  "Building Leaflet preview for ",
  leaflet_tile,
  " / ",
  leaflet_chip_status,
  "..."
)

leaflet_chip_ids <- chip_manifest_csv %>%
  filter(
    mgrs_tile_t == leaflet_tile,
    chip_status == leaflet_chip_status
  ) %>%
  pull(chip_id)




if (length(leaflet_chip_ids) == 0) {
  warning(
    "No chips found for leaflet_tile = ",
    leaflet_tile,
    " and leaflet_chip_status = ",
    leaflet_chip_status
  )
} else {
  leaflet_chips <- chip_manifest_sf %>%
    filter(chip_id %in% leaflet_chip_ids) %>%
    st_transform(4326)

  leaflet_assignments <- density_assignments_for_export %>%
    filter(chip_id %in% leaflet_chip_ids) %>%
    select(
      chip_id,
      sif_row_id,
      Delta_Date,
      mgrs_tile_t,
      chip_status
    )

  leaflet_sif <- sif_projected %>%
    filter(sif_row_id %in% leaflet_assignments$sif_row_id) %>%
    left_join(
      leaflet_assignments,
      by = c("sif_row_id", "Delta_Date", "mgrs_tile_t")
    ) %>%
    mutate(
      sif_label = paste0(
        "<b>SIF row: ", sif_row_id, "</b>",
        "<br>chip_id: ", chip_id,
        "<br>date: ", Delta_Date,
        "<br>tile: ", mgrs_tile_t,
        "<br>status: ", chip_status,
        if ("target_modis_sif" %in% names(.)) {
          paste0("<br>target_modis_sif: ", round(target_modis_sif, 4))
        } else {
          ""
        },
        if ("state" %in% names(.)) {
          paste0("<br>state: ", state)
        } else {
          ""
        },
        if ("hzs" %in% names(.)) {
          paste0("<br>hzs: ", hzs)
        } else {
          ""
        }
      )
    ) %>%
    st_transform(4326)

  leaflet_tile_boundary <- tile_sf %>%
    filter(mgrs_tile_t == leaflet_tile) %>%
    st_transform(4326)

  chip_label_data <- leaflet_chips %>%
    mutate(
      chip_label = paste0(
        "<b>Chip: ", chip_id, "</b>",
        "<br>date: ", Delta_Date,
        "<br>n_sif: ", n_sif,
        "<br>status: ", chip_status,
        "<br>row_off: ", row_off,
        "<br>col_off: ", col_off
      )
    )

  if ("target_modis_sif" %in% names(leaflet_sif)) {
    sif_pal <- leaflet::colorNumeric(
      palette = "viridis",
      domain = leaflet_sif$target_modis_sif,
      na.color = "#cccccc"
    )
    sif_fill <- ~sif_pal(target_modis_sif)
  } else {
    sif_pal <- NULL
    sif_fill <- "#fdae61"
  }

  leaflet_map <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron, group = "Basemap") %>%
    addPolygons(
      data = leaflet_tile_boundary,
      fill = FALSE,
      color = "#111111",
      weight = 3,
      opacity = 1,
      label = ~mgrs_tile_t,
      group = "Sentinel tile boundary"
    ) %>%
    addPolygons(
      data = chip_label_data,
      fillColor = "#f4a3a3",
      fillOpacity = 0.18,
      color = "#2b7bba",
      weight = 1.2,
      opacity = 0.9,
      label = ~lapply(chip_label, htmltools::HTML),
      group = paste0("Chips ", leaflet_chip_status)
    ) %>%
    addPolygons(
      data = leaflet_sif,
      fillColor = sif_fill,
      fillOpacity = 0.65,
      color = "#444444",
      weight = 0.7,
      opacity = 0.85,
      label = ~lapply(sif_label, htmltools::HTML),
      group = "SIF footprints"
    ) %>%
    addLayersControl(
      overlayGroups = c(
        "Sentinel tile boundary",
        paste0("Chips ", leaflet_chip_status),
        "SIF footprints"
      ),
      options = layersControlOptions(collapsed = FALSE)
    )

  if (!is.null(sif_pal)) {
    leaflet_map <- leaflet_map %>%
      addLegend(
        position = "bottomright",
        pal = sif_pal,
        values = leaflet_sif$target_modis_sif,
        title = "target_modis_sif",
        opacity = 0.85
      )
  }

  leaflet_out <- file.path(
    output_dir,
    paste0(
      "sentinel2_leaflet_",
      leaflet_tile,
      "_",
      leaflet_chip_status,
      "_chips.html"
    )
  )

  htmlwidgets::saveWidget(leaflet_map, leaflet_out, selfcontained = TRUE)
  message("Wrote Leaflet preview: ", leaflet_out)
}



