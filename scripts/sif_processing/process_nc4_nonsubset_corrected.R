library(ncdf4)
library(lubridate)
library(tidyverse)


nc_data <- nc_open('data/oco2_helper/oco2_2018_2024/oco2_LtSIF_171231_B11012Ar_221129012634s.nc4')

attributes(nc_data$var)
lon_corners <- ncvar_get(nc_data,"Daily_SIF_757nm")


#-------------------------------------------------------------------------------
# Input/output paths
input_dir <- 'data/oco2_helper/oco2_2018_2024/'
output_file <- 'data/SIF_1_12.rds'
chunk_dir <- 'data/SIF_chunks'

if (!dir.exists(chunk_dir)) {
  dir.create(chunk_dir, recursive = TRUE)
}

files <- list.files(input_dir, pattern = "\\.nc4$", full.names = TRUE)
files <- sort(files)

cat("Found", length(files), "NC4 files.\n")


#-------------------------------------------------------------------------------
# Germany bounding box and target timeframe
lat_min <- 47.302
lat_max <- 54.983
lon_min <- 5.989
lon_max <- 15.017

month_min <- 1
month_max <- 12

vars_to_extract <- c("Daily_SIF_757nm", "Daily_SIF_771nm", "Daily_SIF_740nm", 
                     "SIF_Uncertainty_740", "Science/SIF_Uncertainty_757nm", "Science/SIF_Uncertainty_771nm",
                     "SZA", "VZA", "VAz", "SAz", 
                     "Delta_Time", "Latitude", "Longitude",
                     "Latitude_Corners", "Longitude_Corners", "Quality_Flag",
                     "Meteo/specific_humidity", "Meteo/surface_pressure",
                     "Meteo/temperature_skin", "Meteo/temperature_two_meter",
                     "Meteo/vapor_pressure_deficit", "Metadata/MeasurementMode", 
                     "Science/daily_correction_factor")


#-------------------------------------------------------------------------------
# Helpers
safe_name <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

subset_vector <- function(vals, row_mask, n_selected, short_name, variable_name) {
  vals <- as.vector(vals)

  if (length(vals) == length(row_mask)) {
    vals[row_mask]
  } else if (length(vals) == 1) {
    rep(vals, n_selected)
  } else {
    cat("Unexpected length for", variable_name, "in", short_name, "- filling NA.\n")
    rep(NA, n_selected)
  }
}

extract_one_file <- function(file_path) {
  short_name <- safe_name(file_path)
  cat("Reading:", file_path, "\n")

  nc <- tryCatch(ncdf4::nc_open(file_path), error = function(e) {
    cat("Could not open:", file_path, "\n")
    NULL
  })

  if (is.null(nc)) {
    return(NULL)
  }

  on.exit(ncdf4::nc_close(nc), add = TRUE)

  if (length(nc$var) == 0) {
    cat("Skipping (no variables):", file_path, "\n")
    return(NULL)
  }

  if (!all(c("Latitude", "Longitude", "Delta_Time") %in% names(nc$var))) {
    cat("Missing lat/lon/time in file:", short_name, "\n")
    return(NULL)
  }

  # Read only the filter variables first.
  lat <- ncvar_get(nc, "Latitude")
  lon <- ncvar_get(nc, "Longitude")
  dtime_raw <- ncvar_get(nc, "Delta_Time")
  dtime <- as.POSIXct(dtime_raw, origin = "1990-01-01", tz = "UTC")

  row_mask <- lat >= lat_min & lat <= lat_max &
    lon >= lon_min & lon <= lon_max &
    month(dtime) >= month_min & month(dtime) <= month_max
  row_mask[is.na(row_mask)] <- FALSE

  n_selected <- sum(row_mask)

  if (n_selected == 0) {
    cat("No Germany soundings in February-July in", short_name, "\n")
    return(NULL)
  }

  vars <- list()

  for (v in vars_to_extract) {
    if (v %in% names(nc$var)) {
      vals <- ncvar_get(nc, v)

      if (v == "Latitude_Corners") {
        if (length(dim(vals)) == 1) {
          vars$Lat_corner1 <- vals[1]
          vars$Lat_corner2 <- vals[2]
          vars$Lat_corner3 <- vals[3]
          vars$Lat_corner4 <- vals[4]
        } else if (length(dim(vals)) == 2) {
          vars$Lat_corner1 <- as.vector(vals[1, row_mask])
          vars$Lat_corner2 <- as.vector(vals[2, row_mask])
          vars$Lat_corner3 <- as.vector(vals[3, row_mask])
          vars$Lat_corner4 <- as.vector(vals[4, row_mask])
        } else {
          cat("Unexpected shape for Latitude_Corners in", short_name, "\n")
          vars$Lat_corner1 <- rep(NA, n_selected)
          vars$Lat_corner2 <- rep(NA, n_selected)
          vars$Lat_corner3 <- rep(NA, n_selected)
          vars$Lat_corner4 <- rep(NA, n_selected)
        }

      } else if (v == "Longitude_Corners") {
        if (length(dim(vals)) == 1) {
          vars$Lon_corner1 <- vals[1]
          vars$Lon_corner2 <- vals[2]
          vars$Lon_corner3 <- vals[3]
          vars$Lon_corner4 <- vals[4]
        } else if (length(dim(vals)) == 2) {
          vars$Lon_corner1 <- as.vector(vals[1, row_mask])
          vars$Lon_corner2 <- as.vector(vals[2, row_mask])
          vars$Lon_corner3 <- as.vector(vals[3, row_mask])
          vars$Lon_corner4 <- as.vector(vals[4, row_mask])
        } else {
          cat("Unexpected shape for Longitude_Corners in", short_name, "\n")
          vars$Lon_corner1 <- rep(NA, n_selected)
          vars$Lon_corner2 <- rep(NA, n_selected)
          vars$Lon_corner3 <- rep(NA, n_selected)
          vars$Lon_corner4 <- rep(NA, n_selected)
        }

      } else {
        vars[[v]] <- subset_vector(vals, row_mask, n_selected, short_name, v)
      }
    } else {
      cat("Variable", v, "missing in", short_name, "\n")
      vars[[v]] <- rep(NA, n_selected)
    }
  }

  df <- as.data.frame(vars)
  df$file_id <- short_name
  df$source_file <- basename(file_path)

  df
}


#-------------------------------------------------------------------------------
# Stream through files: open one file, subset it, write filtered chunk, close file.
chunk_files <- character()

for (i in seq_along(files)) {
  cat("Processing file", i, "of", length(files), "\n")

  df <- tryCatch(
    extract_one_file(files[i]),
    error = function(e) {
      cat("Error processing", files[i], ":", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(df) || nrow(df) == 0) {
    next
  }

  chunk_file <- file.path(chunk_dir, paste0(sprintf("%04d", i), "_", safe_name(files[i]), ".rds"))
  saveRDS(df, chunk_file)
  chunk_files <- c(chunk_files, chunk_file)

  cat("Saved", nrow(df), "filtered rows to", chunk_file, "\n")
  rm(df)
  gc(verbose = FALSE)
}


#-------------------------------------------------------------------------------
# Combine only the filtered chunks at the end and write one final RDS file.
if (length(chunk_files) == 0) {
  combined_df <- tibble()
  cat("No matching rows found. Saving empty data frame.\n")
} else {
  cat("Combining", length(chunk_files), "filtered chunks.\n")
  combined_df <- map_dfr(chunk_files, readRDS)
  combined_df$Delta_Time <- as.POSIXct(combined_df$Delta_Time, origin = "1990-01-01", tz = "UTC")
}

saveRDS(combined_df, output_file)
cat("Saved final combined data frame with", nrow(combined_df), "rows to", output_file, "\n")

combined_df <- readRDS('data/SIF_1_12.rds')

summary(combined_df$Daily_SIF_740nm)
summary(combined_df$Latitude)
summary(combined_df$Longitude)
summary(combined_df$Delta_Time)
sort(unique(lubridate::month(combined_df$Delta_Time)))

colSums(is.na(combined_df))
combined_df <- combined_df %>% select(-c(SIF_Uncertainty_740))
combined_df <- drop_na(combined_df)
combined_df <- combined_df %>% filter(Metadata.MeasurementMode %in% c(0,1))
combined_df <- combined_df %>% filter(Quality_Flag %in% c(0,1))
combined_df$Delta_Date <- as.Date(combined_df$Delta_Time)

library(sf)
library(giscoR)

germany_states <- giscoR::gisco_get_nuts(
  country = "DE",
  nuts_level = 1,
  resolution = "01",
  epsg = 4326
) %>%
  select(state = NUTS_NAME, geometry)


# Convert soundings to points using SIF footprint polygon centroids
sif_polygons <- lapply(seq_len(nrow(combined_df)), function(i) {
  coords <- matrix(
    c(
      combined_df$Lon_corner1[i], combined_df$Lat_corner1[i],
      combined_df$Lon_corner2[i], combined_df$Lat_corner2[i],
      combined_df$Lon_corner3[i], combined_df$Lat_corner3[i],
      combined_df$Lon_corner4[i], combined_df$Lat_corner4[i],
      combined_df$Lon_corner1[i], combined_df$Lat_corner1[i]
    ),
    ncol = 2,
    byrow = TRUE
  )

  st_polygon(list(coords))
})

combined_points <- st_sf(
  combined_df,
  geometry = st_centroid(st_transform(st_sfc(sif_polygons, crs = 4326), 3035)) %>%
    st_transform(4326)
)

# Keep only points inside Germany and add state name
combined_df_germany <- combined_points %>%
  st_join(germany_states, join = st_within, left = FALSE)

table(combined_df_germany$state)
table(month(combined_df_germany$Delta_Date))

saveRDS(combined_df_germany, file = 'data/sif_sf_1_12_cleaned.rds')
 

combined_df_germany_27 <- combined_df_germany[month(combined_df_germany$Delta_Date) %in% 2:7,]

saveRDS(combined_df_germany_27, file = 'data/sif_sf_2_7_cleaned.rds')
#-------------------------------------------------------------------------------



# #-----------------------------------------------------------------------------
# 
# #rast('../data/bavaria_ww_gt40.tif')
# #rast('../data/bavaria_ww_gt40.tif')
# #rast('../data/bavaria_ww_gt40.tif')
# 
# high_wheat <- rast('../data/sachsen-Anhalt_ww_gt40.tif')
# high_wheat_ll <- project(high_wheat, "EPSG:4326", method = "near")
# high_wheat_df <- as.data.frame(high_wheat_ll, xy = TRUE, na.rm = TRUE) %>% filter(crop == 1)
# 
# 
# # ggplot() +
# #   geom_sf(data = germany_states, fill = "grey95", color = "grey40", linewidth = 0.25) +
# #   geom_tile(
# #     data = high_wheat_df,
# #     aes(x = x, y = y),
# #     fill = "red",
# #     #alpha = 0.45
# #   ) +
# #   geom_sf(
# #     data = combined_df_germany %>% filter(state == 'BAYERN'),
# #     size = 0.05,
# #     #alpha = 0.2,
# #     color = "darkgreen"
# #   ) +
# #   coord_sf(xlim = c(9, 14), ylim = c(47.45, 50.85),expand = FALSE) +
# #   theme_minimal()
# 
# # Bavaria boundary in UTM 32N
# #NIEDERSACHSEN
# #BAYERN
# bavaria_utm <- germany_states %>%
#   filter(state == "SACHSEN-ANHALT") %>%
#   st_transform(32632)
# 
# bb <- st_bbox(bavaria_utm)
# 
# # 100 km UTM grid spacing
# grid_step <- 100000
# 
# e_seq <- seq(
#   floor(bb["xmin"] / grid_step) * grid_step,
#   ceiling(bb["xmax"] / grid_step) * grid_step,
#   by = grid_step
# )
# 
# n_seq <- seq(
#   floor(bb["ymin"] / grid_step) * grid_step,
#   ceiling(bb["ymax"] / grid_step) * grid_step,
#   by = grid_step
# )
# 
# utm_vertical <- st_sfc(
#   lapply(e_seq, function(e) {
#     st_linestring(matrix(c(e, bb["ymin"], e, bb["ymax"]), ncol = 2, byrow = TRUE))
#   }),
#   crs = 32632
# )
# 
# utm_horizontal <- st_sfc(
#   lapply(n_seq, function(n) {
#     st_linestring(matrix(c(bb["xmin"], n, bb["xmax"], n), ncol = 2, byrow = TRUE))
#   }),
#   crs = 32632
# )
# 
# utm_grid <- st_sf(
#   type = c(rep("easting", length(utm_vertical)), rep("northing", length(utm_horizontal))),
#   geometry = c(utm_vertical, utm_horizontal)
# )
# 
# ggplot() +
#   geom_sf(data = germany_states, fill = "grey95", color = "grey40", linewidth = 0.25) +
#   geom_tile(
#     data = high_wheat_df,
#     aes(x = x, y = y),
#     fill = "red",
#     #alpha = 0.45
#   ) +
#   geom_sf(
#     data = utm_grid,
#     color = "blue",
#     linewidth = 0.35,
#     #alpha = 0.8
#   ) +
#   geom_sf(
#     data = combined_df_germany %>% filter(state == "SACHSEN-ANHALT"),
#     size = 0.05,
#     color = "darkgreen"
#   ) +
#   coord_sf(
#     #xlim = c(9, 14), ylim = c(47.45, 50.85),
#     #xlim = c(6.5, 11.8), ylim = c(51.2, 54),
#     xlim = c(10, 13.6), ylim = c(51, 53.5),
#     expand = FALSE
#   ) +
#   theme_minimal()
