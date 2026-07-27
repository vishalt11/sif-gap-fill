suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(future)
  library(furrr)
})

# Extract polygon-mean predictors for tabular SIF regression.
#
# Output predictors:
#   mean_ndmi, mean_ndvi, mean_evi, mean_nirv, mean_ndre
#   mean_fapar, mean_par, mean_apar
#   active_crop_fraction
#
# Sentinel-2 spectral indices follow the same definitions, 20 m grid and
# land-mask handling used for the CNN predictor chips. The older non-crop
# low-NDVI contamination mask is deliberately not applied here.

# -----------------------------------------------------------------------------
# Configuration

input_csv <- "data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.csv"

output_dir <- "data/sentinel2_polygon_means"
output_rds <- file.path(
  output_dir,
  "9tiles_2_7_M01_QF01_inoutrange_PARrm_polygon_means.rds"
)
output_csv <- file.path(
  output_dir,
  "9tiles_2_7_M01_QF01_inoutrange_PARrm_polygon_means.csv"
)

fapar_dir <- "data/glass_geotiff/fapar"
par_dir <- "data/viirs_vnp18a2_daily_mean_par_germany_native"
crop_type_dir <- "data/crop_type_tif"

parallel_workers <- 2L
buffer_cells <- 2L
sentinel_min_inside_fraction <- 0.99

reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4

evi_valid_min <- -1
evi_valid_max <- 1
evi_denominator_epsilon <- 1e-6

fapar_tiles <- c("h18v03", "h18v04")
fapar_available_doys <- seq(33L, 209L, by = 8L)

par_accepted_qa_codes <- c(1, 2)
par_fill_value <- -1
par_valid_min <- 0
par_valid_max <- 700

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(memfrac = 0.55, progress = 0)

# -----------------------------------------------------------------------------
# General helpers

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, call. = FALSE)
  }
  path
}

require_columns <- function(data, columns, description) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      description,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

normalize_mgrs_tile <- function(tile) {
  tile <- as.character(tile)
  if_else(str_starts(tile, "T"), tile, paste0("T", tile))
}

product_id_from_path <- function(product_path) {
  basename(str_remove(as.character(product_path), "[/\\\\]+$"))
}

wasp_file <- function(product_path, type, band_or_resolution) {
  product_id <- product_id_from_path(product_path)
  filename <- paste0(
    product_id,
    "_",
    type,
    "_",
    band_or_resolution,
    ".tif"
  )
  subdirectory <- if (type == "FLG") "MASKS" else ""
  check_file_exists(file.path(product_path, subdirectory, filename))
}

make_sif_polygon <- function(
  lon1, lat1,
  lon2, lat2,
  lon3, lat3,
  lon4, lat4
) {
  st_polygon(list(rbind(
    c(lon1, lat1),
    c(lon2, lat2),
    c(lon3, lat3),
    c(lon4, lat4),
    c(lon1, lat1)
  )))
}

as_spatvector <- function(polygons) {
  if (inherits(polygons, "SpatVector")) {
    return(polygons)
  }
  terra::vect(st_as_sf(polygons))
}

buffered_polygon_extent <- function(
  raster,
  polygons_same_crs,
  buffer_cells = 2L
) {
  polygons_vect <- as_spatvector(polygons_same_crs)
  polygon_extent <- terra::ext(polygons_vect)
  raster_extent <- terra::ext(raster)

  x_buffer <- abs(terra::xres(raster)) * buffer_cells
  y_buffer <- abs(terra::yres(raster)) * buffer_cells

  crop_xmin <- max(
    terra::xmin(polygon_extent) - x_buffer,
    terra::xmin(raster_extent)
  )
  crop_xmax <- min(
    terra::xmax(polygon_extent) + x_buffer,
    terra::xmax(raster_extent)
  )
  crop_ymin <- max(
    terra::ymin(polygon_extent) - y_buffer,
    terra::ymin(raster_extent)
  )
  crop_ymax <- min(
    terra::ymax(polygon_extent) + y_buffer,
    terra::ymax(raster_extent)
  )

  if (
    !all(is.finite(c(crop_xmin, crop_xmax, crop_ymin, crop_ymax))) ||
      crop_xmin >= crop_xmax ||
      crop_ymin >= crop_ymax
  ) {
    return(NULL)
  }

  terra::ext(crop_xmin, crop_xmax, crop_ymin, crop_ymax)
}

crop_raster_to_polygons <- function(
  raster,
  polygons_same_crs,
  buffer_cells = 2L
) {
  crop_extent <- buffered_polygon_extent(
    raster,
    polygons_same_crs,
    buffer_cells
  )
  if (is.null(crop_extent)) {
    return(NULL)
  }
  terra::crop(raster, crop_extent, snap = "out")
}

extract_layer_means <- function(raster, polygons_same_crs) {
  polygons_same_crs <- st_as_sf(polygons_same_crs)
  extracted <- terra::extract(
    raster,
    terra::vect(polygons_same_crs),
    fun = mean,
    na.rm = TRUE,
    touches = FALSE
  )

  if (nrow(extracted) != nrow(polygons_same_crs)) {
    stop(
      "terra::extract returned ",
      nrow(extracted),
      " rows for ",
      nrow(polygons_same_crs),
      " polygons.",
      call. = FALSE
    )
  }

  values <- as_tibble(extracted) %>%
    select(-ID) %>%
    mutate(
      across(
        everything(),
        ~ {
          value <- as.numeric(.x)
          value[is.nan(value) | is.infinite(value)] <- NA_real_
          value
        }
      )
    )

  bind_cols(
    tibble(sif_row_id = polygons_same_crs$sif_row_id),
    values
  )
}

assert_unique_ids <- function(data, description) {
  if (anyDuplicated(data$sif_row_id)) {
    stop(description, " contains duplicate sif_row_id values.", call. = FALSE)
  }
  invisible(data)
}

normalized_difference <- function(a, b) {
  denominator <- a + b
  terra::ifel(
    abs(denominator) > 1e-6,
    (a - b) / denominator,
    NA
  )
}

valid_fraction_layer <- function(raster) {
  terra::ifel(is.na(raster), 0, 1)
}

# -----------------------------------------------------------------------------
# Crop codes and active-growth calendar

crop_codes <- c(
  winter_wheat = 11L,
  winter_barley = 12L,
  winter_rye = 13L,
  other_winter_cereals = 14L,
  spring_wheat = 21L,
  spring_barley = 22L,
  spring_oat = 23L,
  maize = 30L,
  legumes = 40L,
  potato = 50L,
  sugar_beet = 60L,
  rapeseed = 71L,
  clover_alfalfa = 81L,
  arable_grass = 82L,
  permanent_grassland = 83L,
  vineyard = 90L,
  fruit_trees_and_other_woody_vegetation = 100L,
  hops = 110L,
  other_agricultural_use = 111L
)

active_growth_months <- list(
  `11` = c(2:7, 10:11),
  `12` = c(2:6, 10:11),
  `13` = c(2:7, 10:11),
  `14` = c(2:7, 10:11),
  `21` = 3:8,
  `22` = 3:8,
  `23` = 3:8,
  `30` = 5:10,
  `40` = 4:9,
  `50` = 4:9,
  `60` = 4:10,
  `71` = c(2:7, 9:11),
  `81` = 3:10,
  `82` = 3:11,
  `83` = 3:11,
  `90` = 4:10,
  `100` = 3:10,
  `110` = 4:9,
  `111` = 3:10
)

active_codes_for_month <- function(month) {
  as.integer(names(active_growth_months)[
    map_lgl(active_growth_months, ~ month %in% .x)
  ])
}

active_crop_raster <- function(crop_raster, month) {
  active_codes <- active_codes_for_month(month)
  if (length(active_codes) == 0) {
    return(terra::ifel(is.na(crop_raster), 0, 0))
  }

  active_condition <- crop_raster == active_codes[[1]]
  if (length(active_codes) > 1) {
    for (code in active_codes[-1]) {
      active_condition <- active_condition | crop_raster == code
    }
  }

  # The crop products use NA/background for non-crop in parts of Germany.
  # Those cells are intentionally part of the denominator and therefore zero.
  terra::ifel(
    is.na(crop_raster),
    0,
    terra::ifel(active_condition, 1, 0)
  )
}

# -----------------------------------------------------------------------------
# Input and SIF polygons

corner_columns <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

required_input_columns <- c(
  "Delta_Date",
  "mgrs_tile",
  "product_path",
  corner_columns
)

message("Reading SIF rows from ", input_csv)
sif_data <- readr::read_csv(
  check_file_exists(input_csv),
  show_col_types = FALSE
)
require_columns(sif_data, required_input_columns, "SIF input")

sif_data <- sif_data %>%
  mutate(
    across(all_of(corner_columns), as.numeric),
    source_csv_row = row_number(),
    sif_row_id = row_number(),
    Delta_Date = as.Date(Delta_Date),
    sif_year = lubridate::year(Delta_Date),
    sif_month = lubridate::month(Delta_Date),
    sif_doy = lubridate::yday(Delta_Date),
    mgrs_tile_t = normalize_mgrs_tile(mgrs_tile),
    product_path = as.character(product_path),
    fapar_composite_doy = 1L + 8L * ((sif_doy - 1L) %/% 8L)
  )

if (any(!complete.cases(sif_data[, corner_columns]))) {
  stop("Some SIF rows have missing corner coordinates.", call. = FALSE)
}
if (any(is.na(sif_data$Delta_Date))) {
  stop("Some SIF rows have an invalid Delta_Date.", call. = FALSE)
}
if (any(is.na(sif_data$product_path) | sif_data$product_path == "")) {
  stop("Some SIF rows have a missing product_path.", call. = FALSE)
}

bad_fapar_doys <- sif_data %>%
  filter(!fapar_composite_doy %in% fapar_available_doys) %>%
  distinct(Delta_Date, sif_doy, fapar_composite_doy)

if (nrow(bad_fapar_doys) > 0) {
  stop(
    "Some SIF dates map outside the downloaded FAPAR composite range:\n",
    paste(capture.output(print(bad_fapar_doys, n = Inf)), collapse = "\n"),
    call. = FALSE
  )
}

message("Building ", nrow(sif_data), " SIF footprint polygons.")
sif_geometry <- pmap(
  list(
    sif_data$Lon_corner1, sif_data$Lat_corner1,
    sif_data$Lon_corner2, sif_data$Lat_corner2,
    sif_data$Lon_corner3, sif_data$Lat_corner3,
    sif_data$Lon_corner4, sif_data$Lat_corner4
  ),
  make_sif_polygon
)

sif_sf <- st_sf(
  sif_data,
  geometry = do.call(st_sfc, c(sif_geometry, list(crs = 4326)))
) %>%
  st_make_valid()

sif_lookup <- sif_sf %>%
  st_drop_geometry() %>%
  select(
    sif_row_id,
    Delta_Date,
    sif_year,
    sif_month,
    sif_doy,
    mgrs_tile_t,
    product_path,
    fapar_composite_doy
  )

# -----------------------------------------------------------------------------
# Sentinel raster extent coverage

first_product <- sif_lookup$product_path[[1]]
sentinel_template <- rast(
  wasp_file(first_product, "FRC", "B5")
)
sentinel_crs <- terra::crs(sentinel_template)
rm(sentinel_template)

message("Transforming SIF polygons to the Sentinel CRS once.")
sif_sentinel_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(sentinel_crs)

sentinel_extent_geometry <- function(reference_raster) {
  raster_bbox <- st_bbox(
    c(
      xmin = terra::xmin(reference_raster),
      ymin = terra::ymin(reference_raster),
      xmax = terra::xmax(reference_raster),
      ymax = terra::ymax(reference_raster)
    ),
    crs = sf::st_crs(terra::crs(reference_raster))
  )
  st_as_sfc(raster_bbox)
}

coverage_for_product <- function(target_product) {
  group_rows <- which(sif_lookup$product_path == target_product)
  polygons <- sif_sentinel_crs[group_rows, ]

  reference <- rast(wasp_file(target_product, "FRC", "B5"))
  if (!isTRUE(
    sf::st_crs(terra::crs(reference)) == sf::st_crs(sentinel_crs)
  )) {
    stop(
      "Sentinel CRS differs from the template for ",
      target_product,
      call. = FALSE
    )
  }

  raster_extent <- sentinel_extent_geometry(reference)
  fully_covered <- lengths(
    st_covered_by(st_geometry(polygons), raster_extent)
  ) > 0

  inside_fraction <- rep(1, nrow(polygons))
  edge_rows <- which(!fully_covered)

  # Exact intersections are only needed for the small edge subset.
  if (length(edge_rows) > 0) {
    inside_fraction[edge_rows] <- map_dbl(edge_rows, function(index) {
      footprint <- st_geometry(polygons)[index]
      footprint_area <- as.numeric(st_area(footprint))
      if (!is.finite(footprint_area) || footprint_area <= 0) {
        return(NA_real_)
      }

      intersection <- suppressWarnings(
        st_intersection(footprint, raster_extent)
      )
      if (
        length(intersection) == 0 ||
          all(st_is_empty(intersection))
      ) {
        return(0)
      }

      intersection_area <- sum(as.numeric(st_area(intersection)))
      pmin(1, pmax(0, intersection_area / footprint_area))
    })
  }

  tibble(
    sif_row_id = polygons$sif_row_id,
    sentinel_inside_fraction = inside_fraction,
    sentinel_edge_incomplete = (
      is.na(inside_fraction) |
        inside_fraction < sentinel_min_inside_fraction
    ),
    sentinel_predictor_eligible = (
      !is.na(inside_fraction) &
        inside_fraction >= sentinel_min_inside_fraction
    )
  )
}

product_keys <- sif_lookup %>%
  distinct(product_path) %>%
  arrange(product_path) %>%
  pull(product_path)

message(
  "Checking actual Sentinel raster coverage for ",
  length(product_keys),
  " products."
)

future::plan(future::multisession, workers = parallel_workers)

sentinel_coverage <- furrr::future_map_dfr(
  product_keys,
  coverage_for_product,
  .progress = TRUE,
  .options = furrr::furrr_options(seed = FALSE, scheduling = 1)
)
assert_unique_ids(sentinel_coverage, "Sentinel coverage results")

sif_lookup <- sif_lookup %>%
  left_join(sentinel_coverage, by = "sif_row_id")

# -----------------------------------------------------------------------------
# Sentinel-2 spectral-index means

sentinel_band_resolution <- c(
  B2 = "R1",
  B4 = "R1",
  B8 = "R1",
  B5 = "R2",
  B8A = "R2",
  B11 = "R2"
)

load_masked_band_to_20m <- function(
  product_path,
  band_id,
  flag_r1,
  flag_r2,
  template_20m
) {
  band <- rast(wasp_file(product_path, "FRC", band_id))
  flag <- if (sentinel_band_resolution[[band_id]] == "R1") {
    flag_r1
  } else {
    flag_r2
  }

  band_crop <- terra::crop(
    band,
    terra::ext(template_20m),
    snap = "out"
  )
  flag_crop <- terra::crop(flag, band_crop, snap = "out")

  if (!terra::compareGeom(
    flag_crop,
    band_crop,
    stopOnError = FALSE
  )) {
    flag_crop <- terra::resample(
      flag_crop,
      band_crop,
      method = "near"
    )
  }

  masked_band <- terra::ifel(
    flag_crop == land_flag_value &
      band_crop != reflectance_nodata,
    band_crop / reflectance_quantification_value,
    NA
  )

  if (!terra::compareGeom(
    masked_band,
    template_20m,
    stopOnError = FALSE
  )) {
    masked_band <- terra::resample(
      masked_band,
      template_20m,
      method = "average"
    )
  }

  masked_band
}

extract_sentinel_product <- function(target_product) {
  group_rows <- which(
    sif_lookup$product_path == target_product &
      sif_lookup$sentinel_predictor_eligible
  )
  if (length(group_rows) == 0) {
    return(NULL)
  }

  polygons <- sif_sentinel_crs[group_rows, ]
  polygons_vect <- terra::vect(polygons)

  flag_r1 <- rast(wasp_file(target_product, "FLG", "R1"))
  flag_r2 <- rast(wasp_file(target_product, "FLG", "R2"))

  target_extent <- buffered_polygon_extent(
    flag_r2,
    polygons_vect,
    buffer_cells
  )
  if (is.null(target_extent)) {
    stop(
      "No Sentinel intersection for eligible polygons: ",
      target_product,
      call. = FALSE
    )
  }

  template_20m <- terra::crop(
    flag_r2,
    target_extent,
    snap = "out"
  )

  bands <- map(
    names(sentinel_band_resolution),
    ~ load_masked_band_to_20m(
      target_product,
      .x,
      flag_r1,
      flag_r2,
      template_20m
    )
  )
  names(bands) <- names(sentinel_band_resolution)

  ndmi <- normalized_difference(bands$B8, bands$B11)
  ndvi <- normalized_difference(bands$B8, bands$B4)

  evi_denominator <- (
    bands$B8 +
      6 * bands$B4 -
      7.5 * bands$B2 +
      1
  )
  evi <- terra::ifel(
    abs(evi_denominator) > evi_denominator_epsilon,
    2.5 * (bands$B8 - bands$B4) / evi_denominator,
    NA
  )
  evi <- terra::ifel(
    evi >= evi_valid_min & evi <= evi_valid_max,
    evi,
    NA
  )

  nirv <- bands$B8 * ndvi
  ndre <- normalized_difference(bands$B8A, bands$B5)

  all_indices_valid <- terra::ifel(
    !is.na(ndmi) &
      !is.na(ndvi) &
      !is.na(evi) &
      !is.na(nirv) &
      !is.na(ndre),
    1,
    0
  )

  index_stack <- c(
    ndmi,
    ndvi,
    evi,
    nirv,
    ndre,
    all_indices_valid
  )
  names(index_stack) <- c(
    "mean_ndmi",
    "mean_ndvi",
    "mean_evi",
    "mean_nirv",
    "mean_ndre",
    "sentinel_indices_valid_fraction"
  )

  result <- extract_layer_means(index_stack, polygons)

  rm(
    flag_r1,
    flag_r2,
    template_20m,
    bands,
    ndmi,
    ndvi,
    evi,
    nirv,
    ndre,
    all_indices_valid,
    index_stack
  )
  gc(verbose = FALSE)

  result
}

eligible_product_keys <- sif_lookup %>%
  filter(sentinel_predictor_eligible) %>%
  distinct(product_path) %>%
  arrange(product_path) %>%
  pull(product_path)

message(
  "Extracting five Sentinel indices from ",
  length(eligible_product_keys),
  " products."
)

sentinel_means <- furrr::future_map_dfr(
  eligible_product_keys,
  extract_sentinel_product,
  .progress = TRUE,
  .options = furrr::furrr_options(seed = FALSE, scheduling = 1)
)
assert_unique_ids(sentinel_means, "Sentinel spectral-index results")

# -----------------------------------------------------------------------------
# FAPAR, PAR and APAR means

one_matching_file <- function(
  directory,
  pattern,
  description
) {
  matches <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE
  )
  if (length(matches) != 1) {
    stop(
      "Expected one file for ",
      description,
      "; found ",
      length(matches),
      if (length(matches) > 0) {
        paste0(":\n", paste(matches, collapse = "\n"))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  matches[[1]]
}

fapar_file <- function(year, doy, tile) {
  directory <- file.path(fapar_dir, tile, year)
  pattern <- sprintf(
    "^GLASS09D01\\.V[0-9]+\\.A%d%03d\\.%s\\..*\\.tif$",
    year,
    doy,
    tile
  )
  one_matching_file(
    directory,
    pattern,
    sprintf("FAPAR %d DOY %03d tile %s", year, doy, tile)
  )
}

par_files <- function(date) {
  date_text <- format(as.Date(date), "%Y-%m-%d")
  year <- format(as.Date(date), "%Y")
  directory <- file.path(par_dir, year)

  list(
    par = check_file_exists(file.path(
      directory,
      paste0(
        "VNP18A2.002_",
        date_text,
        "_Daily_Mean_PAR_VIIRS_Sinusoidal_native.tif"
      )
    )),
    qa = check_file_exists(file.path(
      directory,
      paste0(
        "VNP18A2.002_",
        date_text,
        "_PAR_Quality_VIIRS_Sinusoidal_native.tif"
      )
    ))
  )
}

fapar_template <- rast(
  fapar_file(
    min(sif_lookup$sif_year),
    fapar_available_doys[[1]],
    fapar_tiles[[1]]
  )
)
fapar_crs <- terra::crs(fapar_template)
rm(fapar_template)

message("Transforming SIF polygons to the FAPAR CRS once.")
sif_fapar_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(fapar_crs)

read_fapar_mosaic <- function(
  year,
  doy,
  polygons_fapar_crs
) {
  rasters <- map(fapar_tiles, function(tile) {
    raster <- rast(fapar_file(year, doy, tile))
    cropped <- crop_raster_to_polygons(
      raster,
      polygons_fapar_crs,
      buffer_cells
    )
    if (is.null(cropped)) {
      return(NULL)
    }
    terra::ifel(
      cropped >= 0 & cropped <= 1,
      cropped,
      NA
    )
  }) %>%
    compact()

  if (length(rasters) == 0) {
    return(NULL)
  }
  if (length(rasters) == 1) {
    return(rasters[[1]])
  }

  do.call(
    terra::mosaic,
    c(rasters, fun = "mean")
  )
}

extract_radiation_group <- function(
  target_date,
  target_tile
) {
  target_date <- as.Date(target_date, origin = "1970-01-01")
  group_rows <- which(
    sif_lookup$Delta_Date == target_date &
      sif_lookup$mgrs_tile_t == target_tile
  )
  if (length(group_rows) == 0) {
    return(NULL)
  }

  polygons_fapar <- sif_fapar_crs[group_rows, ]
  first_row <- group_rows[[1]]
  target_year <- sif_lookup$sif_year[[first_row]]
  target_doy <- sif_lookup$fapar_composite_doy[[first_row]]

  fapar <- read_fapar_mosaic(
    target_year,
    target_doy,
    polygons_fapar
  )
  if (is.null(fapar)) {
    stop(
      "No FAPAR raster intersects ",
      target_tile,
      " on ",
      target_date,
      call. = FALSE
    )
  }

  paths <- par_files(target_date)
  par <- rast(paths$par)
  qa <- rast(paths$qa)

  if (
    !terra::same.crs(par, qa) ||
      !terra::compareGeom(par, qa, stopOnError = FALSE)
  ) {
    stop(
      "PAR and QA grids differ for ",
      target_date,
      call. = FALSE
    )
  }

  polygons_par <- terra::project(
    terra::vect(polygons_fapar),
    terra::crs(par)
  )
  par_crop <- crop_raster_to_polygons(
    par,
    polygons_par,
    buffer_cells
  )
  if (is.null(par_crop)) {
    stop(
      "No PAR raster intersects ",
      target_tile,
      " on ",
      target_date,
      call. = FALSE
    )
  }

  qa_crop <- terra::crop(qa, par_crop, snap = "out")
  if (!terra::compareGeom(
    par_crop,
    qa_crop,
    stopOnError = FALSE
  )) {
    qa_crop <- terra::resample(
      qa_crop,
      par_crop,
      method = "near"
    )
  }

  par_valid <- terra::ifel(
    (qa_crop == par_accepted_qa_codes[[1]] |
      qa_crop == par_accepted_qa_codes[[2]]) &
      par_crop != par_fill_value &
      par_crop >= par_valid_min &
      par_crop <= par_valid_max,
    par_crop,
    NA
  )

  # Use the cropped FAPAR grid as a common coarse grid. This avoids creating
  # artificial 20 m detail while still calculating mean(FAPAR * PAR), rather
  # than multiplying two polygon means.
  par_on_fapar <- terra::project(
    par_valid,
    fapar,
    method = "bilinear"
  )
  apar <- fapar * par_on_fapar

  radiation_stack <- c(
    fapar,
    par_on_fapar,
    apar,
    valid_fraction_layer(fapar),
    valid_fraction_layer(par_on_fapar),
    valid_fraction_layer(apar)
  )
  names(radiation_stack) <- c(
    "mean_fapar",
    "mean_par",
    "mean_apar",
    "fapar_valid_fraction",
    "par_valid_fraction",
    "apar_valid_fraction"
  )

  result <- extract_layer_means(
    radiation_stack,
    polygons_fapar
  )

  rm(
    fapar,
    par,
    qa,
    par_crop,
    qa_crop,
    par_valid,
    par_on_fapar,
    apar,
    radiation_stack
  )
  gc(verbose = FALSE)

  result
}

radiation_groups <- sif_lookup %>%
  distinct(Delta_Date, mgrs_tile_t) %>%
  arrange(Delta_Date, mgrs_tile_t)

message(
  "Extracting FAPAR, daily PAR and APAR for ",
  nrow(radiation_groups),
  " date/tile groups."
)

radiation_means <- furrr::future_pmap_dfr(
  list(
    radiation_groups$Delta_Date,
    radiation_groups$mgrs_tile_t
  ),
  extract_radiation_group,
  .progress = TRUE,
  .options = furrr::furrr_options(seed = FALSE, scheduling = 1)
)
assert_unique_ids(radiation_means, "Radiation results")

# -----------------------------------------------------------------------------
# Active-crop fraction

crop_file <- function(year) {
  check_file_exists(
    file.path(crop_type_dir, paste0("croptypes_", year, ".tif"))
  )
}

crop_template <- rast(crop_file(min(sif_lookup$sif_year)))
crop_crs <- terra::crs(crop_template)
rm(crop_template)

message("Transforming SIF polygons to the crop-raster CRS once.")
sif_crop_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(crop_crs)

extract_active_crop_group <- function(
  target_year,
  target_month,
  target_tile
) {
  group_rows <- which(
    sif_lookup$sif_year == target_year &
      sif_lookup$sif_month == target_month &
      sif_lookup$mgrs_tile_t == target_tile
  )
  if (length(group_rows) == 0) {
    return(NULL)
  }

  crop_raster <- rast(crop_file(target_year))
  polygons <- sif_crop_crs[group_rows, ]

  if (!isTRUE(
    sf::st_crs(terra::crs(crop_raster)) == sf::st_crs(crop_crs)
  )) {
    polygons <- st_transform(polygons, terra::crs(crop_raster))
  }

  crop_local <- crop_raster_to_polygons(
    crop_raster,
    polygons,
    buffer_cells
  )
  if (is.null(crop_local)) {
    stop(
      "No crop raster intersects ",
      target_tile,
      " for ",
      target_year,
      " month ",
      target_month,
      call. = FALSE
    )
  }

  active_crop <- active_crop_raster(
    crop_local,
    target_month
  )
  names(active_crop) <- "active_crop_fraction"

  result <- extract_layer_means(active_crop, polygons)

  rm(crop_raster, crop_local, active_crop)
  gc(verbose = FALSE)

  result
}

crop_groups <- sif_lookup %>%
  distinct(sif_year, sif_month, mgrs_tile_t) %>%
  arrange(sif_year, sif_month, mgrs_tile_t)

message(
  "Extracting active-crop fractions for ",
  nrow(crop_groups),
  " year/month/tile groups."
)

active_crop_means <- furrr::future_pmap_dfr(
  list(
    crop_groups$sif_year,
    crop_groups$sif_month,
    crop_groups$mgrs_tile_t
  ),
  extract_active_crop_group,
  .progress = TRUE,
  .options = furrr::furrr_options(seed = FALSE, scheduling = 1)
)
assert_unique_ids(active_crop_means, "Active-crop results")

future::plan(future::sequential)

# -----------------------------------------------------------------------------
# Join results and save

predictor_columns <- c(
  "mean_ndmi",
  "mean_ndvi",
  "mean_evi",
  "mean_nirv",
  "mean_ndre",
  "mean_fapar",
  "mean_par",
  "mean_apar",
  "active_crop_fraction"
)

model_sf <- sif_sf %>%
  left_join(sentinel_coverage, by = "sif_row_id") %>%
  left_join(sentinel_means, by = "sif_row_id") %>%
  left_join(radiation_means, by = "sif_row_id") %>%
  left_join(active_crop_means, by = "sif_row_id") %>%
  mutate(
    tabular_predictor_complete = (
      sentinel_predictor_eligible &
        if_all(
          all_of(predictor_columns),
          ~ !is.na(.x) & is.finite(.x)
        )
    )
  )

if (nrow(model_sf) != nrow(sif_data)) {
  stop(
    "Output row count changed from ",
    nrow(sif_data),
    " to ",
    nrow(model_sf),
    ".",
    call. = FALSE
  )
}
if (anyDuplicated(model_sf$sif_row_id)) {
  stop("Output contains duplicate sif_row_id values.", call. = FALSE)
}

saveRDS(model_sf, output_rds)
model_sf %>%
  st_drop_geometry() %>%
  readr::write_csv(output_csv)

message("Input SIF rows: ", nrow(sif_data))
message(
  "Sentinel-complete footprints: ",
  sum(model_sf$sentinel_predictor_eligible, na.rm = TRUE)
)
message(
  "Sentinel edge-incomplete footprints: ",
  sum(model_sf$sentinel_edge_incomplete, na.rm = TRUE)
)
message(
  "Rows with all nine tabular predictors: ",
  sum(model_sf$tabular_predictor_complete, na.rm = TRUE)
)

message("Predictor summaries:")
print(
  model_sf %>%
    st_drop_geometry() %>%
    select(all_of(predictor_columns)) %>%
    summary()
)

message("Saved RDS: ", output_rds)
message("Saved CSV: ", output_csv)
