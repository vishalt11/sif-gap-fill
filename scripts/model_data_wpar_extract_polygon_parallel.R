library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(furrr)
library(future)


#df <- readRDS('data/sif_sf_months2_7_cleaned.rds')

df <- read_csv('data/338k_base_crop_hzs.csv')

# spectral_indices_dir <- "data/spectral_indices_means_nonveg_masked"
fapar_dir <- "data/glass_fapar_modis_250m"
par_dir <- "data/glass_par_modis_005d"

glass_days <- seq(33, 209, by = 8)
glass_years <- 2019:2024
fapar_tiles <- c("h18v03", "h18v04")
parallel_workers <- 2

match_glass_composite_doy <- function(doy) {
  if (is.na(doy)) {
    return(NA_integer_)
  }

  matched_doys <- glass_days[doy >= glass_days & doy <= glass_days + 7]

  if (length(matched_doys) == 0) {
    return(NA_integer_)
  }

  matched_doys[[1]]
}

match_closest_downloaded_doy <- function(doy) {
  if (is.na(doy)) {
    return(NA_integer_)
  }

  glass_days[which.min(abs(glass_days - doy))]
}

# 
# 
# csv_files <- list.files(
#   spectral_indices_dir,
#   pattern = "\\.csv$",
#   full.names = TRUE
# )
# 
# if (length(csv_files) == 0) {
#   stop("No CSV files found in: ", spectral_indices_dir)
# }
# 
# df <- csv_files %>%
#   map_dfr(
#     ~ read_csv(
#       .x,
#       col_types = cols(
#         sif_id = col_character(),
#         mgrs_tile = col_character(),
#         Delta_Date = col_date(),
#         .default = col_guess()
#       ),
#       show_col_types = FALSE
#     )
#   )
# 
# 
# 
# 
corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

#target_states <- c("SACHSEN-ANHALT", "BAYERN", "NIEDERSACHSEN")

df <- df %>%
  st_drop_geometry() %>%
  mutate(
    Delta_Date = as.Date(Delta_Date),
    state = str_to_upper(state)
  ) %>%
  filter(
    #state %in% target_states,
    year(Delta_Date) %in% glass_years
  )

missing_corner_cols <- setdiff(corner_cols, names(df))
if (length(missing_corner_cols) > 0) {
  stop("Missing corner columns: ", paste(missing_corner_cols, collapse = ", "))
}

sif_sf <- df %>%
  mutate(across(all_of(corner_cols), as.numeric)) %>%
  filter(if_all(all_of(corner_cols), ~ !is.na(.x))) %>%
  mutate(
    geometry = pmap(
      list(
        Lon_corner1, Lat_corner1,
        Lon_corner2, Lat_corner2,
        Lon_corner3, Lat_corner3,
        Lon_corner4, Lat_corner4
      ),
      function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
        st_polygon(list(matrix(
          c(
            lon1, lat1,
            lon2, lat2,
            lon3, lat3,
            lon4, lat4,
            lon1, lat1
          ),
          ncol = 2,
          byrow = TRUE
        )))
      }
    )
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid() %>%
  mutate(
    sif_row_id = row_number(),
    sif_year = year(Delta_Date),
    sif_doy = yday(Delta_Date),
    fapar_doy = map_int(sif_doy, match_glass_composite_doy),
    par_doy = map_int(sif_doy, match_closest_downloaded_doy),
    fapar_start_date = as.Date(sprintf("%d-01-01", sif_year)) + fapar_doy - 1,
    par_date = as.Date(sprintf("%d-01-01", sif_year)) + par_doy - 1,
    fapar_day_offset = as.integer(Delta_Date - fapar_start_date),
    par_day_diff = abs(as.integer(Delta_Date - par_date))
  )
# 
# missing_corner_cols <- setdiff(corner_cols, names(df))
# if (length(missing_corner_cols) > 0) {
#   stop("Missing corner columns: ", paste(missing_corner_cols, collapse = ", "))
# }
# 
# sif_sf <- df %>%
#   mutate(across(all_of(corner_cols), as.numeric)) %>%
#   filter(if_all(all_of(corner_cols), ~ !is.na(.x))) %>%
#   mutate(
#     geometry = pmap(
#       list(
#         Lon_corner1, Lat_corner1,
#         Lon_corner2, Lat_corner2,
#         Lon_corner3, Lat_corner3,
#         Lon_corner4, Lat_corner4
#       ),
#       function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
#         st_polygon(list(matrix(
#           c(
#             lon1, lat1,
#             lon2, lat2,
#             lon3, lat3,
#             lon4, lat4,
#             lon1, lat1
#           ),
#           ncol = 2,
#           byrow = TRUE
#         )))
#       }
#     )
#   ) %>%
#   st_as_sf(crs = 4326) %>%
#   st_make_valid() %>%
#   mutate(
#     sif_row_id = row_number(),
#     sif_year = year(Delta_Date),
#     sif_doy = yday(Delta_Date),
#     closest_doy = map_int(sif_doy,~ if (is.na(.x)) NA_integer_ else glass_days[which.min(abs(glass_days - .x))]),
#     closest_glass_date = as.Date(sprintf("%d-01-01", sif_year)) + closest_doy - 1,
#     glass_day_diff = abs(as.integer(Delta_Date - closest_glass_date))
#   )


one_matching_file <- function(dir_path, pattern, description) {
  matches <- list.files(
    dir_path,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(matches) == 0) {
    warning("No file found for ", description)
    return(NA_character_)
  }

  if (length(matches) > 1) {
    stop("Multiple files found for ", description, ":\n", paste(matches, collapse = "\n"))
  }

  matches
}

fapar_file_path <- function(year, doy, tile) {
  tile_dir <- file.path(fapar_dir, tile, year)
  pattern <- sprintf("^GLASS09D01\\.V60\\.A%d%03d\\.%s\\..*\\.hdf$", year, doy, tile)
  description <- sprintf("FAPAR %d DOY %03d tile %s", year, doy, tile)
  one_matching_file(tile_dir, pattern, description)
}

par_file_path <- function(year, doy) {
  year_dir <- file.path(par_dir, year)
  pattern <- sprintf("^GLASS04B01\\.V42\\.A%d%03d\\..*\\.hdf$", year, doy)
  description <- sprintf("PAR %d DOY %03d", year, doy)
  one_matching_file(year_dir, pattern, description)
}

buffered_polygon_extent <- function(raster, polygons_same_crs, buffer_cells = 2) {
  polygons_vect <- vect(st_as_sf(polygons_same_crs))
  polygon_extent <- ext(polygons_vect)
  x_buffer <- abs(xres(raster)) * buffer_cells
  y_buffer <- abs(yres(raster)) * buffer_cells
  raster_extent <- ext(raster)

  crop_xmin <- max(xmin(polygon_extent) - x_buffer, xmin(raster_extent))
  crop_xmax <- min(xmax(polygon_extent) + x_buffer, xmax(raster_extent))
  crop_ymin <- max(ymin(polygon_extent) - y_buffer, ymin(raster_extent))
  crop_ymax <- min(ymax(polygon_extent) + y_buffer, ymax(raster_extent))

  if (crop_xmin >= crop_xmax || crop_ymin >= crop_ymax) {
    return(NULL)
  }

  ext(crop_xmin, crop_xmax, crop_ymin, crop_ymax)
}

crop_raster_to_polygon_extent <- function(raster, polygons_same_crs, buffer_cells = 2) {
  crop_extent <- buffered_polygon_extent(raster, polygons_same_crs, buffer_cells = buffer_cells)

  if (is.null(crop_extent)) {
    return(NULL)
  }

  tryCatch(
    crop(raster, crop_extent),
    error = function(e) {
      warning("Raster crop failed; using full raster. Error: ", conditionMessage(e))
      raster
    }
  )
}

read_fapar_mosaic <- function(year, doy, crop_polygons = NULL, buffer_cells = 2) {
  fapar_rasters <- map(fapar_tiles, function(tile) {
    path <- fapar_file_path(year, doy, tile)

    if (is.na(path)) {
      return(NULL)
    }

    r <- rast(path)
    if (!is.null(crop_polygons)) {
      r <- crop_raster_to_polygon_extent(r, crop_polygons, buffer_cells = buffer_cells)
    }

    if (is.null(r)) {
      return(NULL)
    }

    ifel(r < 0 | r > 1, NA, r)
  })

  fapar_rasters <- compact(fapar_rasters)

  if (length(fapar_rasters) == 0) {
    return(NULL)
  }

  if (length(fapar_rasters) == 1) {
    return(fapar_rasters[[1]])
  }

  do.call(terra::mosaic, c(fapar_rasters, fun = "mean"))
}

read_par_raster <- function(year, doy, crop_polygons = NULL, buffer_cells = 2) {
  path <- par_file_path(year, doy)

  if (is.na(path)) {
    return(NULL)
  }

  r <- rast(path)
  if (!is.null(crop_polygons)) {
    r <- crop_raster_to_polygon_extent(r, crop_polygons, buffer_cells = buffer_cells)
  }

  if (is.null(r)) {
    return(NULL)
  }

  r
}

extract_polygon_mean <- function(raster, polygons_same_crs) {
  if (is.null(raster)) {
    return(rep(NA_real_, nrow(polygons_same_crs)))
  }

  polygons_same_crs <- st_as_sf(polygons_same_crs)

  extracted <- terra::extract(
    raster,
    vect(polygons_same_crs),
    fun = mean,
    na.rm = TRUE
  )

  extracted[[2]]
}

fapar_template <- rast(fapar_file_path(2019, 33, "h18v03"))
par_template <- rast(par_file_path(2019, 33))

message("Preparing polygon lookup once.")
sif_lookup <- sif_sf %>%
  st_drop_geometry() %>%
  select(sif_row_id, sif_year, fapar_doy, par_doy)

message("Transforming polygons to FAPAR CRS once.")
sif_sf_fapar_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(crs(fapar_template))

message("Transforming polygons to PAR CRS once.")
sif_sf_par_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(crs(par_template))

extract_groups <- sif_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(sif_year), !is.na(fapar_doy), !is.na(par_doy)) %>%
  distinct(sif_year, fapar_doy, par_doy) %>%
  arrange(sif_year, fapar_doy, par_doy)

plan(multisession, workers = parallel_workers)

extraction_results <- future_pmap_dfr(
  list(extract_groups$sif_year, extract_groups$fapar_doy, extract_groups$par_doy),
  function(target_year, target_fapar_doy, target_par_doy) {
    message(
      "Extracting GLASS values for ",
      target_year,
      " FAPAR DOY ",
      sprintf("%03d", target_fapar_doy),
      " PAR DOY ",
      sprintf("%03d", target_par_doy)
    )

    group_rows <- which(
      sif_lookup$sif_year == target_year &
        sif_lookup$fapar_doy == target_fapar_doy &
        sif_lookup$par_doy == target_par_doy
    )

    sif_group_fapar <- sif_sf_fapar_crs[group_rows, ]
    sif_group_par <- sif_sf_par_crs[group_rows, ]

    fapar_mosaic <- read_fapar_mosaic(target_year, target_fapar_doy, sif_group_fapar, buffer_cells = 2)
    par_raster <- read_par_raster(target_year, target_par_doy, sif_group_par, buffer_cells = 2)

    tibble(
      sif_row_id = sif_lookup$sif_row_id[group_rows],
      mean_fapar = extract_polygon_mean(fapar_mosaic, sif_group_fapar),
      mean_par = extract_polygon_mean(par_raster, sif_group_par)
    )
  },
  .progress = TRUE,
  .options = furrr_options(seed = FALSE, scheduling = 1)
)

plan(sequential)

sif_sf <- sif_sf %>%
  left_join(extraction_results, by = "sif_row_id") %>%
  mutate(apar = mean_fapar * mean_par)

message("Combined CSV rows: ", nrow(df))
message("SIF polygons built: ", nrow(sif_sf))
message("Rows with FAPAR: ", sum(!is.na(sif_sf$mean_fapar)))
message("Rows with PAR: ", sum(!is.na(sif_sf$mean_par)))

summary(sif_sf$mean_fapar)
summary(sif_sf$mean_par)
summary(sif_sf$apar)

saveRDS(sif_sf, "data/extracted_modis_data/338k_crop_hzs_wpar.rds")

sif_sf %>%
  st_drop_geometry() %>%
  write_csv("data/extracted_modis_data/338k_crop_hzs_wpar.csv")

#-------------------------------------------------------------------------------
# testing fapar and par values

set.seed(42)
test_sif_row_ids <- sif_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(mean_fapar), !is.na(mean_par)) %>%
  slice_sample(n = 3) %>%
  pull(sif_row_id)

test_rows <- which(sif_sf$sif_row_id %in% test_sif_row_ids)

test_results <- map_dfr(test_rows, function(row_index) {
  target_year <- sif_sf$sif_year[row_index]
  target_fapar_doy <- sif_sf$fapar_doy[row_index]
  target_par_doy <- sif_sf$par_doy[row_index]

  fapar_mosaic <- read_fapar_mosaic(target_year, target_fapar_doy)
  par_raster <- read_par_raster(target_year, target_par_doy)
  parallel_mean_fapar <- sif_sf$mean_fapar[row_index]
  parallel_mean_par <- sif_sf$mean_par[row_index]
  parallel_apar <- sif_sf$apar[row_index]
  serial_mean_fapar <- extract_polygon_mean(fapar_mosaic, sif_sf_fapar_crs[row_index, ])
  serial_mean_par <- extract_polygon_mean(par_raster, sif_sf_par_crs[row_index, ])
  serial_apar <- serial_mean_fapar * serial_mean_par

  tibble(
    sif_row_id = sif_sf$sif_row_id[row_index],
    sif_year = target_year,
    fapar_doy = target_fapar_doy,
    par_doy = target_par_doy,
    parallel_mean_fapar = parallel_mean_fapar,
    serial_mean_fapar = serial_mean_fapar,
    fapar_diff = serial_mean_fapar - parallel_mean_fapar,
    parallel_mean_par = parallel_mean_par,
    serial_mean_par = serial_mean_par,
    par_diff = serial_mean_par - parallel_mean_par,
    parallel_apar = parallel_apar,
    serial_apar = serial_apar,
    apar_diff = serial_apar - parallel_apar
  )
})

print(test_results)

#-------------------------------------------------------------------------------











