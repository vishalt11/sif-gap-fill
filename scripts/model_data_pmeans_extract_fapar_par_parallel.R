suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(future)
  library(furrr)
})

# Add polygon-mean GLASS FAPAR, daily VIIRS PAR, APAR and active-crop
# percentage to the existing Sentinel-2 polygon-means table.

# -----------------------------------------------------------------------------
# Configuration

input_csv <- paste0(
  "data/pmeans_model_data/",
  "model_data_pmeans_6tiles_with_modis_sif.csv"
)

output_csv <- paste0(
  "data/pmeans_model_data/",
  "model_data_pmeans_6tiles_with_modis_sif_fapar_par.csv"
)
output_rds <- paste0(
  "data/pmeans_model_data/",
  "model_data_pmeans_6tiles_with_modis_sif_fapar_par.rds"
)

fapar_dir <- "data/glass_geotiff/fapar"
par_dir <- "data/viirs_vnp18a2_daily_mean_par_germany_native"

fapar_tiles <- c("h18v03", "h18v04")
fapar_available_doys <- seq(33L, 209L, by = 8L)

par_accepted_qa_codes <- c(1, 2)
par_fill_value <- -1
par_valid_min <- 0
par_valid_max <- 700

parallel_workers <- 2L
buffer_cells <- 2L

terraOptions(memfrac = 0.50, progress = 0)
options(future.globals.maxSize = 2 * 1024^3)

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

one_matching_file <- function(directory, pattern, description) {
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

buffered_polygon_extent <- function(
  raster,
  polygons_same_crs,
  buffer_cells = 2L
) {
  polygon_extent <- terra::ext(terra::vect(polygons_same_crs))
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

  crop_values <- c(crop_xmin, crop_xmax, crop_ymin, crop_ymax)
  if (
    !all(is.finite(crop_values)) ||
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

extract_polygon_mean <- function(raster, polygons_same_crs, value_name) {
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

  values <- as.numeric(extracted[[2]])
  values[is.nan(values) | is.infinite(values)] <- NA_real_

  tibble(
    sif_row_id = polygons_same_crs$sif_row_id,
    !!value_name := values
  )
}

assert_unique_ids <- function(data, description) {
  if (anyDuplicated(data$sif_row_id)) {
    stop(description, " contains duplicate sif_row_id values.", call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Crop activity calendar

active_growth_months <- list(
  crop_count_winter_wheat = c(2:7, 10:11),
  crop_count_winter_barley = c(2:6, 10:11),
  crop_count_winter_rye = c(2:7, 10:11),
  crop_count_other_winter_cereals = c(2:7, 10:11),
  crop_count_spring_wheat = 3:8,
  crop_count_spring_barley = 3:8,
  crop_count_spring_oat = 3:8,
  crop_count_maize = 5:10,
  crop_count_legumes = 4:9,
  crop_count_potato = 4:9,
  crop_count_sugar_beet = 4:10,
  crop_count_rapeseed = c(2:7, 9:11),
  crop_count_clover_alfalfa = 3:10,
  crop_count_arable_grass = 3:11,
  crop_count_permanent_grassland = 3:11,
  crop_count_vineyard = 4:10,
  crop_count_fruit_trees_and_other_woody_vegetation = 3:10,
  crop_count_hops = 4:9,
  crop_count_other_agricultural_use = 3:10
)

add_active_crop_percentage <- function(data) {
  require_columns(
    data,
    c("month", "crop_pixel_count", names(active_growth_months)),
    "Active-crop input"
  )

  active_growth_pixel_count <- rep(0, nrow(data))

  for (crop_column in names(active_growth_months)) {
    crop_count <- as.numeric(data[[crop_column]])
    active_month <- data$month %in% active_growth_months[[crop_column]]

    active_growth_pixel_count <- active_growth_pixel_count +
      if_else(active_month, crop_count, 0)
  }

  result <- data %>%
    mutate(
      active_growth_pixel_count = active_growth_pixel_count,
      active_crop_pct = if_else(
        is.finite(crop_pixel_count) & crop_pixel_count > 0,
        active_growth_pixel_count / crop_pixel_count,
        NA_real_
      )
    )

  invalid_fraction <- which(
    !is.na(result$active_crop_pct) &
      (
        result$active_crop_pct < -1e-8 |
          result$active_crop_pct > 1 + 1e-8
      )
  )

  if (length(invalid_fraction) > 0) {
    stop(
      length(invalid_fraction),
      " active_crop_pct values fall outside [0, 1]. ",
      "Check crop_pixel_count and the crop-count columns.",
      call. = FALSE
    )
  }

  result
}

# -----------------------------------------------------------------------------
# Raster file helpers

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
  date <- as.Date(date, origin = "1970-01-01")
  date_text <- format(date, "%Y-%m-%d")
  year <- format(date, "%Y")
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

read_fapar_mosaic <- function(year, doy, polygons_fapar_crs) {
  rasters <- map(fapar_tiles, function(tile) {
    raster <- rast(fapar_file(year, doy, tile))
    raster_crop <- crop_raster_to_polygons(
      raster,
      polygons_fapar_crs,
      buffer_cells
    )

    if (is.null(raster_crop)) {
      return(NULL)
    }

    terra::ifel(
      raster_crop >= 0 & raster_crop <= 1,
      raster_crop,
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

  do.call(terra::mosaic, c(rasters, fun = "mean"))
}

# -----------------------------------------------------------------------------
# Main workflow

main <- function() {
  corner_columns <- c(
    "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
    "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
  )

  required_columns <- c(
    "Delta_Date",
    "crop_pixel_count",
    corner_columns,
    names(active_growth_months)
  )

  message("Reading ", input_csv)
  sif_data <- read_csv(
    check_file_exists(input_csv),
    show_col_types = FALSE
  )
  require_columns(sif_data, required_columns, "Input CSV")

  sif_data <- sif_data %>%
    select(-any_of(c(
      "sif_row_id",
      "source_csv_row",
      "month",
      "fapar_composite_doy",
      "mean_fapar",
      "mean_par",
      "apar",
      "active_growth_pixel_count",
      "active_crop_pct"
    ))) %>%
    mutate(
      across(all_of(corner_columns), as.numeric),
      Delta_Date = as.Date(Delta_Date),
      source_csv_row = row_number(),
      sif_row_id = row_number(),
      month = lubridate::month(Delta_Date),
      extraction_year = lubridate::year(Delta_Date),
      extraction_doy = lubridate::yday(Delta_Date),
      fapar_composite_doy = (
        1L + 8L * ((extraction_doy - 1L) %/% 8L)
      )
    )

  if (any(is.na(sif_data$Delta_Date))) {
    stop("Some rows have an invalid Delta_Date.", call. = FALSE)
  }
  if (any(!complete.cases(sif_data[, corner_columns]))) {
    stop("Some rows have missing SIF corner coordinates.", call. = FALSE)
  }

  invalid_fapar_dates <- sif_data %>%
    filter(!fapar_composite_doy %in% fapar_available_doys) %>%
    distinct(Delta_Date, extraction_doy, fapar_composite_doy)

  if (nrow(invalid_fapar_dates) > 0) {
    stop(
      "Some dates map outside the downloaded FAPAR range:\n",
      paste(
        capture.output(print(invalid_fapar_dates, n = Inf)),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  message("Building ", nrow(sif_data), " SIF polygons.")
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
    sif_data %>% select(sif_row_id),
    geometry = do.call(
      st_sfc,
      c(sif_geometry, list(crs = 4326))
    )
  ) %>%
    st_make_valid()

  sif_lookup <- sif_data %>%
    select(
      sif_row_id,
      Delta_Date,
      extraction_year,
      fapar_composite_doy
    )

  first_year <- min(sif_lookup$extraction_year)
  fapar_template <- rast(fapar_file(
    first_year,
    fapar_available_doys[[1]],
    fapar_tiles[[1]]
  ))
  fapar_crs <- terra::crs(fapar_template)
  rm(fapar_template)

  first_par_paths <- par_files(min(sif_lookup$Delta_Date))
  par_template <- rast(first_par_paths$par)
  par_crs <- terra::crs(par_template)
  rm(par_template)

  message("Transforming all polygons to the FAPAR and PAR CRSs once.")
  sif_fapar_crs <- sif_sf %>%
    st_transform(fapar_crs)
  sif_par_crs <- sif_sf %>%
    st_transform(par_crs)

  extract_fapar_group <- function(target_year, target_doy) {
    target_year <- as.integer(target_year)
    target_doy <- as.integer(target_doy)

    group_rows <- which(
      sif_lookup$extraction_year == target_year &
        sif_lookup$fapar_composite_doy == target_doy
    )
    polygons <- sif_fapar_crs[group_rows, ]

    fapar <- read_fapar_mosaic(
      target_year,
      target_doy,
      polygons
    )
    if (is.null(fapar)) {
      stop(
        "No FAPAR raster intersects year ",
        target_year,
        " DOY ",
        sprintf("%03d", target_doy),
        ".",
        call. = FALSE
      )
    }

    result <- extract_polygon_mean(
      fapar,
      polygons,
      "mean_fapar"
    )

    rm(fapar)
    gc(verbose = FALSE)
    result
  }

  extract_par_group <- function(target_date) {
    target_date <- as.Date(target_date, origin = "1970-01-01")
    group_rows <- which(sif_lookup$Delta_Date == target_date)
    polygons <- sif_par_crs[group_rows, ]

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
        ".",
        call. = FALSE
      )
    }

    par_crop <- crop_raster_to_polygons(
      par,
      polygons,
      buffer_cells
    )
    if (is.null(par_crop)) {
      stop(
        "No PAR raster intersects SIF polygons on ",
        target_date,
        ".",
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

    accepted_qa <- qa_crop == par_accepted_qa_codes[[1]]
    for (qa_code in par_accepted_qa_codes[-1]) {
      accepted_qa <- accepted_qa | qa_crop == qa_code
    }

    par_valid <- terra::ifel(
      accepted_qa &
        par_crop != par_fill_value &
        par_crop >= par_valid_min &
        par_crop <= par_valid_max,
      par_crop,
      NA
    )

    result <- extract_polygon_mean(
      par_valid,
      polygons,
      "mean_par"
    )

    rm(par, qa, par_crop, qa_crop, accepted_qa, par_valid)
    gc(verbose = FALSE)
    result
  }

  fapar_groups <- sif_lookup %>%
    distinct(extraction_year, fapar_composite_doy) %>%
    arrange(extraction_year, fapar_composite_doy)

  par_dates <- sif_lookup %>%
    distinct(Delta_Date) %>%
    arrange(Delta_Date) %>%
    pull(Delta_Date)

  message(
    "Extracting FAPAR for ",
    nrow(fapar_groups),
    " year/composite groups."
  )

  future::plan(
    future::multisession,
    workers = parallel_workers
  )
  on.exit(future::plan(future::sequential), add = TRUE)

  fapar_means <- furrr::future_map2_dfr(
    fapar_groups$extraction_year,
    fapar_groups$fapar_composite_doy,
    extract_fapar_group,
    .progress = TRUE,
    .options = furrr::furrr_options(
      seed = FALSE,
      scheduling = 1
    )
  )
  assert_unique_ids(fapar_means, "FAPAR results")

  message(
    "Extracting daily PAR for ",
    length(par_dates),
    " dates."
  )
  par_means <- furrr::future_map_dfr(
    par_dates,
    extract_par_group,
    .progress = TRUE,
    .options = furrr::furrr_options(
      seed = FALSE,
      scheduling = 1
    )
  )
  assert_unique_ids(par_means, "PAR results")

  future::plan(future::sequential)

  output_data <- sif_data %>%
    left_join(fapar_means, by = "sif_row_id") %>%
    left_join(par_means, by = "sif_row_id") %>%
    mutate(apar = mean_fapar * mean_par) %>%
    add_active_crop_percentage() %>%
    select(-extraction_year, -extraction_doy)

  if (nrow(output_data) != nrow(sif_data)) {
    stop(
      "Output row count changed from ",
      nrow(sif_data),
      " to ",
      nrow(output_data),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(output_data$sif_row_id)) {
    stop("Output contains duplicate sif_row_id values.", call. = FALSE)
  }

  message("Input rows: ", nrow(sif_data))
  message(
    "Rows with FAPAR: ",
    sum(!is.na(output_data$mean_fapar))
  )
  message(
    "Rows with PAR: ",
    sum(!is.na(output_data$mean_par))
  )
  message(
    "Rows with APAR: ",
    sum(!is.na(output_data$apar))
  )
  message(
    "Rows with active_crop_pct: ",
    sum(!is.na(output_data$active_crop_pct))
  )

  print(summary(output_data$mean_fapar))
  print(summary(output_data$mean_par))
  print(summary(output_data$apar))
  print(summary(output_data$active_crop_pct))

  saveRDS(output_data, output_rds)
  write_csv(output_data, output_csv)

  message("Saved RDS: ", output_rds)
  message("Saved CSV: ", output_csv)
}

main()
