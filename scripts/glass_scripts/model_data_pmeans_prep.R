library(tidyverse)
library(sf)
library(terra)
library(lubridate)

spectral_indices_dir <- "data/spectral_indices_means_nonveg_masked"

csv_files <- list.files(
  spectral_indices_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("No CSV files found in: ", spectral_indices_dir)
}

df <- csv_files %>%
  map_dfr(
    ~ read_csv(
      .x,
      col_types = cols(
        sif_id = col_character(),
        mgrs_tile = col_character(),
        Delta_Date = col_date(),
        .default = col_guess()
      ),
      show_col_types = FALSE
    )
  )

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
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
  st_make_valid()

#----- FAPAR and PAR extraction-----------

glass_days <- seq(33, 209, by = 8)
fapar_tiles <- c("h18v03", "h18v04")
fapar_dir <- "data/glass_fapar_modis_250m"
par_dir <- "data/glass_par_modis_005d"

closest_glass_doy <- function(doy) {
  if (is.na(doy)) {
    return(NA_integer_)
  }

  glass_days[which.min(abs(glass_days - doy))]
}

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

read_fapar_mosaic <- function(year, doy) {
  fapar_rasters <- map(fapar_tiles, function(tile) {
    tile_dir <- file.path(fapar_dir, tile, year)
    pattern <- sprintf("^GLASS09D01\\.V60\\.A%d%03d\\.%s\\..*\\.hdf$", year, doy, tile)
    description <- sprintf("FAPAR %d DOY %03d tile %s", year, doy, tile)
    path <- one_matching_file(tile_dir, pattern, description)

    if (is.na(path)) {
      return(NULL)
    }

    r <- rast(path)
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

read_par_raster <- function(year, doy) {
  year_dir <- file.path(par_dir, year)
  pattern <- sprintf("^GLASS04B01\\.V42\\.A%d%03d\\..*\\.hdf$", year, doy)
  description <- sprintf("PAR %d DOY %03d", year, doy)
  path <- one_matching_file(year_dir, pattern, description)

  if (is.na(path)) {
    return(NULL)
  }

  rast(path)
}

extract_polygon_mean <- function(raster, polygons_sf) {
  polygons_raster_crs <- st_transform(polygons_sf, crs(raster))

  extracted <- terra::extract(
    raster,
    vect(polygons_raster_crs),
    fun = mean,
    na.rm = TRUE
  )

  extracted[[2]]
}

sif_sf <- sif_sf %>%
  mutate(
    sif_row_id = row_number(),
    sif_year = year(Delta_Date),
    sif_doy = yday(Delta_Date),
    closest_doy = map_int(sif_doy, closest_glass_doy),
    closest_glass_date = as.Date(sprintf("%d-01-01", sif_year)) + closest_doy - 1,
    glass_day_diff = abs(as.integer(Delta_Date - closest_glass_date)),
    mean_fapar = NA_real_,
    mean_par = NA_real_,
    apar = NA_real_
  )

extract_groups <- sif_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(sif_year), !is.na(closest_doy)) %>%
  distinct(sif_year, closest_doy) %>%
  arrange(sif_year, closest_doy)

for (i in seq_len(nrow(extract_groups))) {
  target_year <- extract_groups$sif_year[i]
  target_doy <- extract_groups$closest_doy[i]

  message("Extracting GLASS values for ", target_year, " DOY ", sprintf("%03d", target_doy))

  group_rows <- which(
    sif_sf$sif_year == target_year &
      sif_sf$closest_doy == target_doy
  )
  sif_group <- sif_sf[group_rows, ]

  fapar_mosaic <- read_fapar_mosaic(target_year, target_doy)
  if (!is.null(fapar_mosaic)) {
    sif_sf$mean_fapar[group_rows] <- extract_polygon_mean(fapar_mosaic, sif_group)
  }

  par_raster <- read_par_raster(target_year, target_doy)
  if (!is.null(par_raster)) {
    sif_sf$mean_par[group_rows] <- extract_polygon_mean(par_raster, sif_group)
  }
}

sif_sf <- sif_sf %>%
  mutate(apar = mean_fapar * mean_par)

message("Combined CSV rows: ", nrow(df))
message("SIF polygons built: ", nrow(sif_sf))
message("Rows with FAPAR: ", sum(!is.na(sif_sf$mean_fapar)))
message("Rows with PAR: ", sum(!is.na(sif_sf$mean_par)))

sif_sf <- sif_sf %>%
  mutate(apar = mean_fapar * mean_par)

summary(sif_sf$mean_fapar)
summary(sif_sf$mean_par)
summary(sif_sf$apar)

saveRDS(sif_sf,"data/model_data_wpar.rds")

sif_sf %>%
  st_drop_geometry() %>%
  write_csv("data/model_data_wpar.csv")
