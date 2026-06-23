library(tidyverse)
library(sf)
library(terra)
library(arrow)

state_prefix <- "ba"
sif_file <- file.path("data", paste0(state_prefix, "_sif_mgrs_crop_composition.csv"))
wasp_tile_root <- "data/geodes_wasp_zips"
crop_type_dir <- "data/crop_type_tif"

target_tile <- "32UNA"

output_dir <- file.path(
  "data/area_to_point_nonveg_masked",
  paste0(state_prefix, "_sif_", target_tile, "_wasp_area_to_point_20m")
)
parquet_dir <- file.path(output_dir, "parquet")
polygon_output_parquet <- file.path(parquet_dir, "polygon_targets.parquet")
manifest_output_parquet <- file.path(parquet_dir, "pixel_table_manifest.parquet")

reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4
non_crop_code <- 0
winter_wheat_code <- 11
non_crop_low_ndvi_threshold <- 0.45
contaminated_20m_fraction_threshold <- 0.5

target_bands_10m <- c("B2", "B3", "B4", "B8")
target_bands_20m <- c("B5", "B6", "B7", "B8A", "B11", "B12")

wasp_root <- file.path(wasp_tile_root, target_tile)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parquet_dir, recursive = TRUE, showWarnings = FALSE)

wasp_file <- function(product_dir, product_id, type, band_or_res = NULL) {
  filename <- if (is.null(band_or_res)) {
    paste0(product_id, "_", type, ".tif")
  } else {
    paste0(product_id, "_", type, "_", band_or_res, ".tif")
  }

  subdir <- if (type %in% c("FLG", "WGT", "DTS")) "MASKS" else ""
  file.path(product_dir, subdir, filename)
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, call. = FALSE)
  }

  path
}

check_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    stop("Missing directory: ", path, call. = FALSE)
  }

  path
}

crop_type_file <- function(year) {
  file.path(crop_type_dir, paste0("croptypes_", year, ".tif"))
}

product_month <- function(product_id) {
  product_date <- str_match(product_id, "^SENTINEL2[A-Z]?_([0-9]{8})-")[, 2]

  if (is.na(product_date)) {
    stop("Could not parse product date from: ", product_id, call. = FALSE)
  }

  format(as.Date(product_date, format = "%Y%m%d"), "%Y-%m")
}

find_wasp_products <- function(wasp_root, target_tile) {
  check_dir_exists(wasp_root)

  year_dirs <- tibble(year_dir = list.dirs(wasp_root, recursive = FALSE, full.names = TRUE)) %>%
    mutate(year = basename(year_dir)) %>%
    filter(str_detect(year, "^[0-9]{4}$"))

  year_dirs %>%
    mutate(product_dir_outer = map(year_dir, ~ list.dirs(.x, recursive = FALSE, full.names = TRUE))) %>%
    select(-year_dir) %>%
    unnest(product_dir_outer) %>%
    transmute(
      product_year = as.integer(year),
      product_id = basename(product_dir_outer),
      product_dir_outer,
      product_dir_nested = file.path(product_dir_outer, product_id),
      product_dir = if_else(dir.exists(product_dir_nested), product_dir_nested, product_dir_outer),
      year_month = map_chr(product_id, product_month)
    ) %>%
    filter(
      str_detect(product_id, paste0("_T", target_tile, "_")),
      dir.exists(product_dir)
    ) %>%
    arrange(year_month)
}

make_sif_polygon <- function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
  st_polygon(list(rbind(
    c(lon1, lat1),
    c(lon2, lat2),
    c(lon3, lat3),
    c(lon4, lat4),
    c(lon1, lat1)
  )))
}

make_sif_polygons <- function(sif_df) {
  sif_polygon_geometry <- pmap(
    list(
      sif_df$Lon_corner1, sif_df$Lat_corner1,
      sif_df$Lon_corner2, sif_df$Lat_corner2,
      sif_df$Lon_corner3, sif_df$Lat_corner3,
      sif_df$Lon_corner4, sif_df$Lat_corner4
    ),
    make_sif_polygon
  )

  st_sf(
    sif_df,
    geometry = do.call(st_sfc, c(sif_polygon_geometry, list(crs = 4326)))
  ) %>%
    st_make_valid()
}

load_scaled_land_band <- function(product_dir, product_id, band_id, flg, sif_vect) {
  band <- rast(check_file_exists(wasp_file(product_dir, product_id, "FRC", band_id)))
  band <- crop(band, sif_vect)

  flg_crop <- crop(flg, band)

  band <- ifel(
    flg_crop == land_flag_value & band != reflectance_nodata,
    band / reflectance_quantification_value,
    NA
  )

  rm(flg_crop)
  gc(verbose = FALSE)

  band
}

align_categorical_raster <- function(categorical_raster, template) {
  if (!terra::same.crs(categorical_raster, template)) {
    return(project(categorical_raster, template, method = "near"))
  }

  categorical_raster <- crop(categorical_raster, template)
  resample(categorical_raster, template, method = "near")
}

downsample_to_20m <- function(raster_10m, template_20m) {
  resample(raster_10m, template_20m, method = "average")
}

normalized_difference <- function(a, b) {
  denominator <- a + b
  ifel(denominator != 0, (a - b) / denominator, NA)
}

safe_ratio <- function(numerator, denominator) {
  ifel(denominator != 0, numerator / denominator, NA)
}

build_contamination_masks <- function(product_year, product_dir, product_id, flg_r1, flg_r2, sif_vect) {
  message("Building non-crop low-NDVI contamination masks for ", product_year, ".")

  b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
  b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
  ndvi <- normalized_difference(b8, b4)

  crop_raster <- rast(check_file_exists(crop_type_file(product_year)))
  crop_10m <- align_categorical_raster(crop_raster, ndvi)

  non_crop <- is.na(crop_10m) | crop_10m == non_crop_code
  contamination_10m <- ifel(
    non_crop & !is.na(ndvi) & ndvi <= non_crop_low_ndvi_threshold,
    1,
    0
  )

  template_20m <- crop(flg_r2, sif_vect)
  contamination_fraction_20m <- resample(contamination_10m, template_20m, method = "average")
  contamination_20m <- ifel(
    contamination_fraction_20m > contaminated_20m_fraction_threshold,
    1,
    0
  )

  rm(
    b8, b4, ndvi, crop_raster, crop_10m, non_crop,
    template_20m, contamination_fraction_20m
  )
  gc(verbose = FALSE)

  list(
    mask_10m = contamination_10m,
    mask_20m = contamination_20m
  )
}

apply_contamination_mask <- function(raster_layer, contamination_mask) {
  mask_for_layer <- crop(contamination_mask, raster_layer)

  if (!compareGeom(mask_for_layer, raster_layer, stopOnError = FALSE)) {
    mask_for_layer <- resample(mask_for_layer, raster_layer, method = "near")
  }

  raster_layer <- ifel(mask_for_layer == 1, NA, raster_layer)

  rm(mask_for_layer)
  gc(verbose = FALSE)

  raster_layer
}

build_feature_stack_20m <- function(product_year, product_dir, product_id, flg_r1, flg_r2, contamination_masks, sif_vect) {
  message("Building 20 m feature stack for ", product_id, ".")

  template_20m <- crop(flg_r2, sif_vect)

  b2 <- downsample_to_20m(load_scaled_land_band(product_dir, product_id, "B2", flg_r1, sif_vect), template_20m)
  b3 <- downsample_to_20m(load_scaled_land_band(product_dir, product_id, "B3", flg_r1, sif_vect), template_20m)
  b4 <- downsample_to_20m(load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect), template_20m)
  b8 <- downsample_to_20m(load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect), template_20m)

  b5 <- load_scaled_land_band(product_dir, product_id, "B5", flg_r2, sif_vect)
  b6 <- load_scaled_land_band(product_dir, product_id, "B6", flg_r2, sif_vect)
  b7 <- load_scaled_land_band(product_dir, product_id, "B7", flg_r2, sif_vect)
  b8a <- load_scaled_land_band(product_dir, product_id, "B8A", flg_r2, sif_vect)
  b11 <- load_scaled_land_band(product_dir, product_id, "B11", flg_r2, sif_vect)
  b12 <- load_scaled_land_band(product_dir, product_id, "B12", flg_r2, sif_vect)

  band_stack <- c(b2, b3, b4, b8, b5, b6, b7, b8a, b11, b12)
  names(band_stack) <- paste0("pixel_", str_to_lower(c(target_bands_10m, target_bands_20m)))

  ndvi <- normalized_difference(b8, b4)
  ndre <- normalized_difference(b8, b5)
  ndre8a <- normalized_difference(b8a, b6)
  psri <- safe_ratio(b4 - b2, b6)
  osavi <- ifel((b8 + b4 + 0.16) != 0, 1.16 * (b8 - b4) / (b8 + b4 + 0.16), NA)
  ndwi <- normalized_difference(b3, b8)
  nirv <- ndvi * b8
  tcari <- ifel(b4 != 0, 3 * ((b5 - b4) - 0.2 * (b5 - b3) * (b5 / b4)), NA)
  ndmi <- normalized_difference(b8, b11)
  msi <- safe_ratio(b11, b8)
  ndmi_swir2 <- normalized_difference(b8, b12)
  msi_swir2 <- safe_ratio(b12, b8)
  nmdi <- ifel((b8 + (b11 - b12)) != 0, (b8 - (b11 - b12)) / (b8 + (b11 - b12)), NA)

  index_stack <- c(ndvi, ndre, ndre8a, psri, osavi, ndwi, nirv, tcari, ndmi, msi, ndmi_swir2, msi_swir2, nmdi)
  names(index_stack) <- paste0(
    "pixel_",
    c("ndvi", "ndre", "ndre8a", "psri", "osavi", "ndwi", "nirv", "tcari", "ndmi", "msi", "ndmi_swir2", "msi_swir2", "nmdi")
  )

  crop_raster <- rast(check_file_exists(crop_type_file(product_year)))
  crop_20m <- align_categorical_raster(crop_raster, template_20m)
  names(crop_20m) <- "pixel_crop_code"

  feature_stack <- c(band_stack, index_stack, crop_20m)
  feature_stack <- apply_contamination_mask(feature_stack, contamination_masks$mask_20m)

  rm(
    template_20m, b2, b3, b4, b8, b5, b6, b7, b8a, b11, b12,
    ndvi, ndre, ndre8a, psri, osavi, ndwi, nirv, tcari,
    ndmi, msi, ndmi_swir2, msi_swir2, nmdi,
    band_stack, index_stack, crop_raster, crop_20m
  )
  gc(verbose = FALSE)

  feature_stack
}

extract_pixel_membership <- function(feature_stack, sif_polygons_utm) {
  sif_vect <- vect(sif_polygons_utm)

  pixel_values <- terra::extract(
    feature_stack,
    sif_vect,
    cells = TRUE,
    xy = TRUE,
    touches = FALSE
  ) %>%
    as_tibble() %>%
    rename(
      extract_id = ID,
      raster_cell = cell,
      pixel_x = x,
      pixel_y = y
    ) %>%
    filter(if_any(matches("^pixel_(b|nd|nirv|osavi|msi|nmdi|psri|tcari)"), ~ !is.na(.x))) %>%
    group_by(extract_id) %>%
    mutate(
      pixel_index_in_polygon = row_number(),
      polygon_pixel_count_20m = n(),
      pixel_weight_equal = 1 / polygon_pixel_count_20m
    ) %>%
    ungroup() %>%
    mutate(
      pixel_is_winter_wheat = !is.na(pixel_crop_code) & pixel_crop_code == winter_wheat_code,
      pixel_is_crop = !is.na(pixel_crop_code) & pixel_crop_code != non_crop_code
    )

  rm(sif_vect)
  gc(verbose = FALSE)

  pixel_values
}

process_product <- function(product_id, product_dir, product_year, year_month) {
  sif_rows_month <- sif_rows %>%
    filter(year_month == .env$year_month) %>%
    mutate(extract_id = row_number())

  if (nrow(sif_rows_month) == 0) {
    message("Skipping ", year_month, ": no matching SIF rows.")
    return(NULL)
  }

  message("Processing ", year_month, " / ", product_id, " with ", nrow(sif_rows_month), " SIF rows.")

  sif_polygons <- make_sif_polygons(sif_rows_month)

  flg_r1 <- rast(check_file_exists(wasp_file(product_dir, product_id, "FLG", "R1")))
  flg_r2 <- rast(check_file_exists(wasp_file(product_dir, product_id, "FLG", "R2")))
  target_crs <- crs(flg_r2)

  sif_polygons_utm <- sif_polygons %>%
    st_transform(target_crs)

  sif_vect <- vect(sif_polygons_utm)

  contamination_masks <- build_contamination_masks(
    product_year = product_year,
    product_dir = product_dir,
    product_id = product_id,
    flg_r1 = flg_r1,
    flg_r2 = flg_r2,
    sif_vect = sif_vect
  )

  feature_stack <- build_feature_stack_20m(
    product_year = product_year,
    product_dir = product_dir,
    product_id = product_id,
    flg_r1 = flg_r1,
    flg_r2 = flg_r2,
    contamination_masks = contamination_masks,
    sif_vect = sif_vect
  )

  pixel_membership <- extract_pixel_membership(feature_stack, sif_polygons_utm) %>%
    left_join(
      sif_rows_month %>%
        select(
          extract_id,
          sif_extract_id,
          sif_id,
          Delta_Date,
          Latitude,
          Longitude,
          Daily_SIF_740nm,
          Quality_Flag,
          Metadata.MeasurementMode,
          state,
          mgrs_tile,
          crop_pixel_count,
          starts_with("crop_count_"),
          ww_pct
        ),
      by = "extract_id"
    ) %>%
    mutate(
      wasp_product_id = product_id,
      wasp_year_month = year_month,
      sif_month = as.integer(format(Delta_Date, "%m")),
      sif_year = as.integer(format(Delta_Date, "%Y")),
      non_crop_low_ndvi_threshold = non_crop_low_ndvi_threshold,
      contaminated_20m_fraction_threshold = contaminated_20m_fraction_threshold
    ) %>%
    select(-extract_id)

  output_parquet <- file.path(
    parquet_dir,
    paste0(state_prefix, "_sif_", target_tile, "_", year_month, "_area_to_point_pixels_20m.parquet")
  )

  write_parquet(pixel_membership, output_parquet)

  manifest_row <- tibble(
    wasp_year_month = year_month,
    wasp_product_id = product_id,
    product_year = product_year,
    n_sif_polygons = nrow(sif_rows_month),
    n_pixel_rows = nrow(pixel_membership),
    n_polygons_with_pixels = n_distinct(pixel_membership$sif_extract_id),
    pixel_parquet = output_parquet
  )

  rm(
    flg_r1, flg_r2, sif_polygons, sif_polygons_utm, sif_vect,
    contamination_masks, feature_stack, pixel_membership
  )
  gc(verbose = FALSE)

  manifest_row
}

wasp_products <- find_wasp_products(wasp_root, target_tile)

if (nrow(wasp_products) == 0) {
  stop("No local WASP products found under ", wasp_root, call. = FALSE)
}

sif_rows <- read_csv(check_file_exists(sif_file), show_col_types = FALSE) %>%
  mutate(
    Delta_Date = as.Date(Delta_Date),
    year_month = format(Delta_Date, "%Y-%m")
  ) %>%
  filter(mgrs_tile == target_tile) %>%
  mutate(sif_extract_id = row_number())

if (nrow(sif_rows) == 0) {
  stop(
    "No SIF rows found for tile ", target_tile,
    " in ", sif_file,
    call. = FALSE
  )
}

message("Selected ", nrow(sif_rows), " SIF rows for ", target_tile, ".")
message("Found ", nrow(wasp_products), " local WASP products under ", wasp_root, ".")

missing_sif_months <- setdiff(unique(sif_rows$year_month), wasp_products$year_month)

if (length(missing_sif_months) > 0) {
  warning(
    "No local WASP product found for SIF month(s): ",
    paste(sort(missing_sif_months), collapse = ", ")
  )
}

matched_sif_rows <- sif_rows %>%
  semi_join(wasp_products, by = "year_month") %>%
  mutate(
    month = as.integer(format(Delta_Date, "%m")),
    year = as.integer(format(Delta_Date, "%Y"))
  )

write_parquet(matched_sif_rows, polygon_output_parquet)

pixel_table_manifest <- pmap(
  wasp_products %>% select(product_id, product_dir, product_year, year_month),
  process_product
) %>%
  compact() %>%
  bind_rows()

if (nrow(pixel_table_manifest) == 0) {
  stop("No SIF rows matched the local WASP products under ", wasp_root, call. = FALSE)
}

write_parquet(pixel_table_manifest, manifest_output_parquet)

message("Saved polygon targets parquet to ", polygon_output_parquet)
message("Saved pixel table manifest parquet to ", manifest_output_parquet)
message("Pixel parquet tables were written under ", parquet_dir)

