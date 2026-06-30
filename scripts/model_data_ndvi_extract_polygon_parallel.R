library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(furrr)
library(future)


#df <- readRDS("data/sif_sf_months2_7_cleaned.rds")

df <- read_csv('data/338k_base_crop_hzs.csv')

ndvi_dir <- "data/glass_ndvi_modis_250m"
output_dir <- "data/extracted_modis_data"

glass_days <- seq(33, 209, by = 8)
glass_years <- 2019:2024
ndvi_tiles <- c("h18v03", "h18v04")
parallel_workers <- 2
buffer_cells <- 2

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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
    sif_area_km2 = as.numeric(st_area(st_transform(geometry, 3035))) / 1e6,
    sif_row_id = row_number(),
    sif_year = year(Delta_Date),
    sif_doy = yday(Delta_Date),
    ndvi_doy = map_int(sif_doy, match_glass_composite_doy),
    ndvi_start_date = as.Date(sprintf("%d-01-01", sif_year)) + ndvi_doy - 1,
    ndvi_day_offset = as.integer(Delta_Date - ndvi_start_date)
  )


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

ndvi_file_path <- function(year, doy, tile) {
  tile_dir <- file.path(ndvi_dir, tile, year)
  pattern <- sprintf("^GLASS13D01\\.V10\\.A%d%03d\\.%s\\..*\\.hdf$", year, doy, tile)
  description <- sprintf("NDVI %d DOY %03d tile %s", year, doy, tile)
  one_matching_file(tile_dir, pattern, description)
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

read_ndvi_mosaic <- function(year, doy, crop_polygons = NULL, buffer_cells = 2) {
  ndvi_rasters <- map(ndvi_tiles, function(tile) {
    path <- ndvi_file_path(year, doy, tile)

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

    ifel(r < -1 | r > 1, NA, r)
  })

  ndvi_rasters <- compact(ndvi_rasters)

  if (length(ndvi_rasters) == 0) {
    return(NULL)
  }

  if (length(ndvi_rasters) == 1) {
    return(ndvi_rasters[[1]])
  }

  do.call(terra::mosaic, c(ndvi_rasters, fun = "mean"))
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

ndvi_template <- rast(ndvi_file_path(2019, 33, "h18v03"))

message("Preparing polygon lookup once.")
sif_lookup <- sif_sf %>%
  st_drop_geometry() %>%
  select(sif_row_id, sif_year, ndvi_doy)

message("Transforming polygons to NDVI CRS once.")
sif_sf_ndvi_crs <- sif_sf %>%
  select(sif_row_id) %>%
  st_transform(crs(ndvi_template))

extract_groups <- sif_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(sif_year), !is.na(ndvi_doy)) %>%
  distinct(sif_year, ndvi_doy) %>%
  arrange(sif_year, ndvi_doy)

plan(multisession, workers = parallel_workers)

extraction_results <- future_pmap_dfr(
  list(extract_groups$sif_year, extract_groups$ndvi_doy),
  function(target_year, target_doy) {
    message("Extracting NDVI values for ", target_year, " DOY ", sprintf("%03d", target_doy))

    group_rows <- which(
      sif_lookup$sif_year == target_year &
        sif_lookup$ndvi_doy == target_doy
    )

    sif_group_ndvi <- sif_sf_ndvi_crs[group_rows, ]
    ndvi_mosaic <- read_ndvi_mosaic(target_year, target_doy, sif_group_ndvi, buffer_cells = buffer_cells)

    tibble(
      sif_row_id = sif_lookup$sif_row_id[group_rows],
      mean_ndvi = extract_polygon_mean(ndvi_mosaic, sif_group_ndvi)
    )
  },
  .progress = TRUE,
  .options = furrr_options(seed = FALSE, scheduling = 1)
)

plan(sequential)

sif_sf <- sif_sf %>%
  left_join(extraction_results, by = "sif_row_id")

message("Input rows after state/year filter: ", nrow(df))
message("SIF polygons built: ", nrow(sif_sf))
message("Rows with NDVI: ", sum(!is.na(sif_sf$mean_ndvi)))

summary(sif_sf$mean_ndvi)

saveRDS(sif_sf, file.path(output_dir, "338k_crop_hzs_ndvi.rds"))

sif_sf %>%
  st_drop_geometry() %>%
  write_csv(file.path(output_dir, "338k_crop_hzs_ndvi.csv"))

#-------------------------------------------------------------------------------
# testing NDVI values

set.seed(42)
test_sif_row_ids <- sif_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(mean_ndvi)) %>%
  slice_sample(n = 5) %>%
  pull(sif_row_id)

test_rows <- which(sif_sf$sif_row_id %in% test_sif_row_ids)

test_results <- map_dfr(test_rows, function(row_index) {
  target_year <- sif_sf$sif_year[row_index]
  target_doy <- sif_sf$ndvi_doy[row_index]

  ndvi_mosaic <- read_ndvi_mosaic(target_year, target_doy)
  parallel_mean_ndvi <- sif_sf$mean_ndvi[row_index]
  serial_mean_ndvi <- extract_polygon_mean(ndvi_mosaic, sif_sf_ndvi_crs[row_index, ])

  tibble(
    sif_row_id = sif_sf$sif_row_id[row_index],
    sif_year = target_year,
    ndvi_doy = target_doy,
    parallel_mean_ndvi = parallel_mean_ndvi,
    serial_mean_ndvi = serial_mean_ndvi,
    ndvi_diff = serial_mean_ndvi - parallel_mean_ndvi
  )
})

print(test_results)

#-------------------------------------------------------------------------------
