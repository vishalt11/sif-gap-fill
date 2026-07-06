library(tidyverse)
library(sf)
library(lubridate)
library(leaflet)
library(terra)

# Fixed-grid diagnostic for multi-footprint CNN chips.
#
# Option 1:
#   Build a regular projected grid over the SIF footprint extent.
#   Assign each SIF footprint centroid to exactly one chip cell.
#   Count how many same-date SIF footprints fall in each chip.
#
# This is only a grouping diagnostic. It does not read MODIS rasters or make NPZ
# files. The actual Python chip-prep script can later use these chip/date groups
# to rasterize multiple footprint masks inside each predictor chip.

sf::sf_use_s2(FALSE)

input_csv <- "data/extracted_modis_data/modis_2_7_bin_uncertainity_corrected_M0QF0_gt0_6area_cnn.csv"
output_dir <- "data/cnn_multisif_chip_diagnostics"
evi_dir <- "data/glass_evi_modis_250m"
ndvi_dir <- "data/glass_ndvi_modis_250m"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# MODIS 250 m is nominal. Your GLASS MODIS sinusoidal rasters are closer to
# 231.656 m. For a first grouping diagnostic, keep the user-facing 250 m chip
# size. Change this to 231.6564 if you want exact GLASS pixel dimensions.
pixel_size_m <- 250
chip_pixel_sizes <- c(16, 24, 32)
glass_tiles <- c("h18v03", "h18v04")
glass_start_doy <- 33
glass_end_doy <- 209
glass_step_days <- 8
modis_index_valid_range <- c(-0.2, 1)

# Use the same MODIS sinusoidal projection as the GLASS FAPAR/EVI/NDVI tiles.
grid_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

optional_summary_cols <- c(
  "state", "hzs", "source_file",
  "target_modis_sif", "Daily_SIF_757nm", "Daily_SIF_771nm",
  "sif_area_km2_evi"
)

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

interval_glass_doy <- function(dates) {
  glass_doys <- seq(glass_start_doy, glass_end_doy, by = glass_step_days)
  sif_doy <- yday(dates)
  sif_doy <- pmin(pmax(sif_doy, glass_start_doy), glass_end_doy)
  glass_idx <- floor((sif_doy - glass_start_doy) / glass_step_days) + 1
  glass_doys[glass_idx]
}

one_matching_glass_file <- function(base_dir, product_prefix, year, doy, tile) {
  year_dir <- file.path(base_dir, tile, as.character(year))

  if (!dir.exists(year_dir)) {
    return(NA_character_)
  }

  file_pattern <- sprintf(
    "^%s.*A%04d%03d\\.%s\\..*\\.(hdf|tif|tiff)$",
    product_prefix,
    year,
    doy,
    tile
  )

  matches <- list.files(
    year_dir,
    pattern = file_pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(matches) == 0) {
    return(NA_character_)
  }

  matches[[1]]
}

read_glass_mosaic <- function(base_dir, product_prefix, year, doy) {
  files <- map_chr(
    glass_tiles,
    ~ one_matching_glass_file(base_dir, product_prefix, year, doy, .x)
  ) %>%
    discard(is.na)

  if (length(files) == 0) {
    warning(
      "No ",
      product_prefix,
      " files found for ",
      year,
      " DOY ",
      sprintf("%03d", doy)
    )
    return(NULL)
  }

  raster_list <- map(files, terra::rast)

  if (length(raster_list) == 1) {
    return(raster_list[[1]])
  }

  do.call(terra::mosaic, c(raster_list, list(fun = "mean")))
}

clean_modis_index <- function(raster_layer) {
  terra::ifel(
    raster_layer < modis_index_valid_range[[1]] |
      raster_layer > modis_index_valid_range[[2]],
    NA,
    raster_layer
  )
}

crop_project_modis_for_leaflet <- function(raster_layer, chip_polygons_projected) {
  chip_vect <- terra::vect(chip_polygons_projected)

  raster_layer %>%
    terra::crop(chip_vect, snap = "out") %>%
    terra::mask(chip_vect) %>%
    clean_modis_index() %>%
    terra::project("EPSG:4326", method = "near")
}

add_modis_leaflet_layers <- function(
  leaflet_map,
  product_name,
  base_dir,
  product_prefix,
  raster_dates,
  chip_polygons_projected,
  raster_palette,
  raster_opacity,
  raster_max_bytes
) {
  added_groups <- character()

  for (i in seq_len(nrow(raster_dates))) {
    target_year <- raster_dates$sif_year[[i]]
    target_doy <- raster_dates$composite_doy[[i]]

    chips_for_raster <- chip_polygons_projected %>%
      filter(
        sif_year == target_year,
        composite_doy == target_doy
      )

    product_raster <- read_glass_mosaic(
      base_dir = base_dir,
      product_prefix = product_prefix,
      year = target_year,
      doy = target_doy
    )

    if (is.null(product_raster)) {
      next
    }

    product_raster_leaflet <- crop_project_modis_for_leaflet(
      product_raster,
      chips_for_raster
    )

    names(product_raster_leaflet) <- "value"

    group_name <- sprintf(
      "%s %d DOY %03d",
      product_name,
      target_year,
      target_doy
    )

    if (all(is.na(terra::values(product_raster_leaflet, mat = FALSE)))) {
      warning("Skipping empty raster layer: ", group_name)
      next
    }

    message("Adding raster layer: ", group_name)

    leaflet_map <- leaflet_map %>%
      addRasterImage(
        product_raster_leaflet,
        colors = raster_palette,
        opacity = raster_opacity,
        group = group_name,
        project = FALSE,
        maxBytes = raster_max_bytes
      )

    added_groups <- c(added_groups, group_name)
  }

  list(
    map = leaflet_map,
    groups = added_groups
  )
}

required_cols <- c("Delta_Date", corner_cols)

message("Reading: ", input_csv)
df <- read_csv(
  input_csv,
  col_types = cols(
    Delta_Date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

message("Rows read: ", nrow(df))

sif_sf <- df %>%
  mutate(
    sif_row_id = row_number(),
    Delta_Date = as.Date(Delta_Date),
    across(all_of(corner_cols), as.numeric)
  ) %>%
  filter(
    !is.na(Delta_Date),
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

sif_projected <- sif_sf %>%
  st_transform(grid_crs)

sif_centroids <- suppressWarnings(st_centroid(sif_projected))
centroid_xy <- st_coordinates(sif_centroids)

sif_centroids_tbl <- sif_centroids %>%
  st_drop_geometry() %>%
  mutate(
    centroid_x = centroid_xy[, "X"],
    centroid_y = centroid_xy[, "Y"]
  )

project_bbox <- st_bbox(sif_projected)

assign_to_fixed_grid <- function(chip_pixels) {
  chip_size_m <- chip_pixels * pixel_size_m

  # Snap the grid origin to an exact chip-size multiple so chip IDs are stable
  # and reproducible, independent of tiny bbox changes.
  grid_xmin <- floor(project_bbox[["xmin"]] / chip_size_m) * chip_size_m
  grid_ymin <- floor(project_bbox[["ymin"]] / chip_size_m) * chip_size_m

  sif_centroids_tbl %>%
    mutate(
      chip_pixels = chip_pixels,
      chip_size_m = chip_size_m,
      chip_col = floor((centroid_x - grid_xmin) / chip_size_m),
      chip_row = floor((centroid_y - grid_ymin) / chip_size_m),
      chip_xmin = grid_xmin + chip_col * chip_size_m,
      chip_ymin = grid_ymin + chip_row * chip_size_m,
      chip_xmax = chip_xmin + chip_size_m,
      chip_ymax = chip_ymin + chip_size_m,
      chip_id = paste0(
        "px", chip_pixels,
        "_r", chip_row,
        "_c", chip_col
      ),
      chip_date_id = paste0(chip_id, "_", format(Delta_Date, "%Y%m%d"))
    )
}

message("Assigning SIF centroids to fixed grids...")
sif_chip_assignments <- map_dfr(chip_pixel_sizes, assign_to_fixed_grid)

write_csv(
  sif_chip_assignments,
  file.path(output_dir, "fixed_grid_sif_assignments.csv")
)

summarise_chip_dates <- function(assignments) {
  has_state <- "state" %in% names(assignments)
  has_hzs <- "hzs" %in% names(assignments)
  has_source_file <- "source_file" %in% names(assignments)
  has_target_modis_sif <- "target_modis_sif" %in% names(assignments)
  has_sif_area <- "sif_area_km2_evi" %in% names(assignments)

  group_cols <- c(
    "chip_pixels", "chip_size_m", "Delta_Date", "chip_id", "chip_date_id",
    "chip_row", "chip_col", "chip_xmin", "chip_ymin", "chip_xmax", "chip_ymax"
  )

  assignments %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_sif = n(),
      n_state = if (has_state) n_distinct(state, na.rm = TRUE) else NA_integer_,
      states = if (has_state) paste(sort(unique(na.omit(state))), collapse = ";") else NA_character_,
      n_hzs = if (has_hzs) n_distinct(hzs, na.rm = TRUE) else NA_integer_,
      hzs_values = if (has_hzs) paste(sort(unique(na.omit(hzs))), collapse = ";") else NA_character_,
      n_source_file = if (has_source_file) n_distinct(source_file, na.rm = TRUE) else NA_integer_,
      mean_target_modis_sif = if (has_target_modis_sif) mean(target_modis_sif, na.rm = TRUE) else NA_real_,
      min_target_modis_sif = if (has_target_modis_sif) min(target_modis_sif, na.rm = TRUE) else NA_real_,
      max_target_modis_sif = if (has_target_modis_sif) max(target_modis_sif, na.rm = TRUE) else NA_real_,
      mean_sif_area_km2 = if (has_sif_area) mean(sif_area_km2_evi, na.rm = TRUE) else NA_real_,
      min_sif_area_km2 = if (has_sif_area) min(sif_area_km2_evi, na.rm = TRUE) else NA_real_,
      max_sif_area_km2 = if (has_sif_area) max(sif_area_km2_evi, na.rm = TRUE) else NA_real_,
      sif_row_ids = paste(sif_row_id, collapse = ","),
      .groups = "drop"
    )
}

chip_date_summary <- summarise_chip_dates(sif_chip_assignments) %>%
  arrange(chip_pixels, Delta_Date, desc(n_sif), chip_id)

write_csv(
  chip_date_summary,
  file.path(output_dir, "fixed_grid_chip_date_summary.csv")
)

chip_count_distribution <- chip_date_summary %>%
  count(chip_pixels, n_sif, name = "n_chip_dates") %>%
  group_by(chip_pixels) %>%
  mutate(
    pct_chip_dates = n_chip_dates / sum(n_chip_dates),
    cum_chip_dates_ge_n = rev(cumsum(rev(n_chip_dates))),
    pct_chip_dates_ge_n = cum_chip_dates_ge_n / sum(n_chip_dates)
  ) %>%
  ungroup() %>%
  arrange(chip_pixels, n_sif)

write_csv(
  chip_count_distribution,
  file.path(output_dir, "fixed_grid_chip_count_distribution.csv")
)

threshold_summary <- chip_date_summary %>%
  group_by(chip_pixels, chip_size_m) %>%
  summarise(
    n_chip_dates = n(),
    n_sif_rows_assigned = sum(n_sif),
    max_sif_per_chip_date = max(n_sif),
    mean_sif_per_chip_date = mean(n_sif),
    median_sif_per_chip_date = median(n_sif),
    chip_dates_ge_2 = sum(n_sif >= 2),
    chip_dates_ge_3 = sum(n_sif >= 3),
    chip_dates_ge_5 = sum(n_sif >= 5),
    chip_dates_ge_10 = sum(n_sif >= 10),
    sif_rows_in_chip_dates_ge_2 = sum(n_sif[n_sif >= 2]),
    sif_rows_in_chip_dates_ge_3 = sum(n_sif[n_sif >= 3]),
    sif_rows_in_chip_dates_ge_5 = sum(n_sif[n_sif >= 5]),
    sif_rows_in_chip_dates_ge_10 = sum(n_sif[n_sif >= 10]),
    pct_chip_dates_ge_2 = chip_dates_ge_2 / n_chip_dates,
    pct_chip_dates_ge_3 = chip_dates_ge_3 / n_chip_dates,
    pct_chip_dates_ge_5 = chip_dates_ge_5 / n_chip_dates,
    pct_chip_dates_ge_10 = chip_dates_ge_10 / n_chip_dates,
    pct_sif_rows_in_chip_dates_ge_2 = sif_rows_in_chip_dates_ge_2 / n_sif_rows_assigned,
    pct_sif_rows_in_chip_dates_ge_3 = sif_rows_in_chip_dates_ge_3 / n_sif_rows_assigned,
    pct_sif_rows_in_chip_dates_ge_5 = sif_rows_in_chip_dates_ge_5 / n_sif_rows_assigned,
    pct_sif_rows_in_chip_dates_ge_10 = sif_rows_in_chip_dates_ge_10 / n_sif_rows_assigned,
    .groups = "drop"
  ) %>%
  arrange(chip_pixels)

write_csv(
  threshold_summary,
  file.path(output_dir, "fixed_grid_threshold_summary.csv")
)

chip_date_sf <- chip_date_summary %>%
  mutate(
    geometry = pmap(
      list(chip_xmin, chip_ymin, chip_xmax, chip_ymax),
      make_square_polygon
    )
  ) %>%
  st_as_sf(crs = grid_crs)

saveRDS(
  chip_date_sf,
  file.path(output_dir, "fixed_grid_chip_date_polygons.rds")
)

message("\nThreshold summary:")
print(threshold_summary)

message("\nCount distribution:")
print(chip_count_distribution[chip_count_distribution$chip_pixels == 16,], n = 100)

message("\nWrote outputs to: ", output_dir)
message("Main files:")
message("  - fixed_grid_threshold_summary.csv")
message("  - fixed_grid_chip_count_distribution.csv")
message("  - fixed_grid_chip_date_summary.csv")
message("  - fixed_grid_sif_assignments.csv")
message("  - fixed_grid_chip_date_polygons.rds")

#-------------------------------------------------------------------------------
# Leaflet preview: sample chip-date cells and their same-date SIF footprints

leaflet_chip_pixels <- 16
leaflet_min_sif <- 5
leaflet_sample_n <- 15
leaflet_seed <- 42
leaflet_add_modis_rasters <- TRUE
leaflet_modis_raster_opacity <- 0.65
leaflet_raster_max_bytes <- 80 * 1024 * 1024

leaflet_candidates <- chip_date_sf %>%
  filter(
    chip_pixels == leaflet_chip_pixels,
    n_sif >= leaflet_min_sif
  )

if (nrow(leaflet_candidates) == 0) {
  warning(
    "No chip-date cells found for leaflet preview with chip_pixels = ",
    leaflet_chip_pixels,
    " and n_sif >= ",
    leaflet_min_sif
  )
} else {
  set.seed(leaflet_seed)

  leaflet_chips <- leaflet_candidates %>%
    slice_sample(n = min(leaflet_sample_n, nrow(leaflet_candidates))) %>%
    mutate(
      sif_year = year(Delta_Date),
      sif_doy = yday(Delta_Date),
      composite_doy = interval_glass_doy(Delta_Date),
      composite_date = make_date(sif_year, 1, 1) + days(composite_doy - 1),
      composite_id = paste0(
        sif_year,
        "_",
        sprintf("%03d", composite_doy)
      )
    )

  leaflet_assignments <- sif_chip_assignments %>%
    filter(
      chip_pixels == leaflet_chip_pixels,
      chip_date_id %in% leaflet_chips$chip_date_id
    ) %>%
    select(sif_row_id, chip_id, chip_date_id, chip_pixels)

  leaflet_sif <- sif_projected %>%
    inner_join(leaflet_assignments, by = "sif_row_id") %>%
    st_transform(4326)

  leaflet_chips_map <- leaflet_chips %>%
    st_transform(4326) %>%
    mutate(
      chip_popup = paste0(
        "<b>", chip_date_id, "</b>",
        "<br>Date: ", Delta_Date,
        "<br>Chip: ", chip_pixels, "x", chip_pixels,
        "<br>SIF DOY: ", sif_doy,
        "<br>EVI/NDVI composite DOY: ", sprintf("%03d", composite_doy),
        "<br>EVI/NDVI composite date: ", composite_date,
        "<br>SIF footprints: ", n_sif,
        "<br>States: ", states,
        "<br>HZS: ", hzs_values
      )
    )

  leaflet_raster_dates <- leaflet_chips %>%
    st_drop_geometry() %>%
    distinct(sif_year, composite_doy, composite_date) %>%
    arrange(sif_year, composite_doy)

  leaflet_sif_map <- leaflet_sif %>%
    mutate(
      sif_popup = paste0(
        "<b>SIF row: ", sif_row_id, "</b>",
        "<br>Date: ", Delta_Date,
        "<br>Chip-date: ", chip_date_id,
        if ("target_modis_sif" %in% names(.)) {
          paste0("<br>target_modis_sif: ", round(target_modis_sif, 4))
        } else {
          ""
        },
        if ("sif_area_km2_evi" %in% names(.)) {
          paste0("<br>area km2: ", round(sif_area_km2_evi, 4))
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
    )

  chip_pal <- colorNumeric(
    palette = "YlOrRd",
    domain = leaflet_chips_map$n_sif,
    na.color = "#cccccc"
  )

  if ("target_modis_sif" %in% names(leaflet_sif_map)) {
    sif_pal <- colorNumeric(
      palette = "RdYlGn",
      domain = leaflet_sif_map$target_modis_sif,
      na.color = "#999999"
    )
    sif_fill <- ~sif_pal(target_modis_sif)
  } else {
    sif_fill <- "#ff7f00"
  }

  index_pal <- colorNumeric(
    palette = "YlGn",
    domain = modis_index_valid_range,
    na.color = "transparent"
  )

  leaflet_map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Esri imagery") %>%
    addProviderTiles(providers$CartoDB.Positron, group = "CartoDB positron")

  raster_overlay_groups <- character()

  if (leaflet_add_modis_rasters) {
    evi_layers <- add_modis_leaflet_layers(
      leaflet_map = leaflet_map,
      product_name = "EVI",
      base_dir = evi_dir,
      product_prefix = "GLASS14D01",
      raster_dates = leaflet_raster_dates,
      chip_polygons_projected = leaflet_chips,
      raster_palette = index_pal,
      raster_opacity = leaflet_modis_raster_opacity,
      raster_max_bytes = leaflet_raster_max_bytes
    )

    leaflet_map <- evi_layers$map
    raster_overlay_groups <- c(raster_overlay_groups, evi_layers$groups)

    ndvi_layers <- add_modis_leaflet_layers(
      leaflet_map = leaflet_map,
      product_name = "NDVI",
      base_dir = ndvi_dir,
      product_prefix = "GLASS13D01",
      raster_dates = leaflet_raster_dates,
      chip_polygons_projected = leaflet_chips,
      raster_palette = index_pal,
      raster_opacity = leaflet_modis_raster_opacity,
      raster_max_bytes = leaflet_raster_max_bytes
    )

    leaflet_map <- ndvi_layers$map
    raster_overlay_groups <- c(raster_overlay_groups, ndvi_layers$groups)
  }

  leaflet_map <- leaflet_map %>%
    addPolygons(
      data = leaflet_chips_map,
      group = "sampled chip boxes",
      color = "#2b6cb0",
      weight = 2,
      fillColor = ~chip_pal(n_sif),
      fillOpacity = 0.18,
      popup = ~chip_popup,
      label = ~paste0(chip_date_id, " | n=", n_sif)
    ) %>%
    addPolygons(
      data = leaflet_sif_map,
      group = "SIF footprints",
      color = "#111111",
      weight = 1,
      fillColor = sif_fill,
      fillOpacity = 0.55,
      popup = ~sif_popup,
      label = ~paste0("SIF row ", sif_row_id, " | ", chip_date_id)
    ) %>%
    addLegend(
      position = "bottomright",
      pal = chip_pal,
      values = leaflet_chips_map$n_sif,
      title = "SIFs per chip-date",
      opacity = 0.8
    ) %>%
    addLegend(
      position = "bottomleft",
      pal = index_pal,
      values = modis_index_valid_range,
      title = "EVI / NDVI",
      opacity = leaflet_modis_raster_opacity
    )

  if (length(raster_overlay_groups) > 0) {
    leaflet_map <- leaflet_map %>%
      hideGroup(raster_overlay_groups)
  }

  leaflet_map <- leaflet_map %>%
    addLayersControl(
      baseGroups = c("Esri imagery", "CartoDB positron"),
      overlayGroups = c(
        "sampled chip boxes",
        "SIF footprints",
        raster_overlay_groups
      ),
      options = layersControlOptions(collapsed = FALSE)
    )

  leaflet_html <- file.path(
    output_dir,
    sprintf(
      "fixed_grid_leaflet_sample_%dpx_min%d_%dchips.html",
      leaflet_chip_pixels,
      leaflet_min_sif,
      min(leaflet_sample_n, nrow(leaflet_candidates))
    )
  )

  htmlwidgets::saveWidget(
    leaflet_map,
    file = leaflet_html,
    selfcontained = TRUE
  )

  message("\nLeaflet preview wrote: ", leaflet_html)
  message(
    "Sampled ",
    nrow(leaflet_chips),
    " chip-date cells and ",
    nrow(leaflet_sif_map),
    " SIF footprints."
  )
}

