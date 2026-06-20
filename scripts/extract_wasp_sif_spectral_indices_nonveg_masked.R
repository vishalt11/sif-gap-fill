library(tidyverse)
library(sf)
library(terra)

state_prefix <- "ba"
sif_file <- file.path("data", paste0(state_prefix, "_sif_mgrs_crop_composition.rds"))
wasp_tile_root <- "data/geodes_wasp_zips"
crop_type_dir <- "data/crop_type_tif"

target_tile <- "32UNA"
target_indices <- c("ndvi", "ndre", "ndre8a", "psri", "osavi", "ndwi", "nirv", "tcari")

output_rds <- file.path(
  "data/spectral_indices_means_nonveg_masked",
  paste0(state_prefix, "_sif_", target_tile, "_wasp_spectral_indices_nonveg_masked.rds")
)
output_csv <- file.path(
  "data/spectral_indices_means_nonveg_masked",
  paste0(state_prefix, "_sif_", target_tile, "_wasp_spectral_indices_nonveg_masked.csv")
)

reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4
non_crop_code <- 0
non_crop_low_ndvi_threshold <- 0.45
contaminated_20m_fraction_threshold <- 0.5

wasp_root <- file.path(wasp_tile_root, target_tile)

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

downsample_to_20m <- function(band_10m, template_20m) {
  resample(band_10m, template_20m, method = "average")
}

normalized_difference <- function(a, b) {
  denominator <- a + b
  ifel(denominator != 0, (a - b) / denominator, NA)
}

apply_contamination_mask <- function(index_raster, contamination_mask) {
  mask_for_index <- crop(contamination_mask, index_raster)

  if (!compareGeom(mask_for_index, index_raster, stopOnError = FALSE)) {
    mask_for_index <- resample(mask_for_index, index_raster, method = "near")
  }

  index_raster <- ifel(mask_for_index == 1, NA, index_raster)

  rm(mask_for_index)
  gc(verbose = FALSE)

  index_raster
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

extract_index_summary <- function(
  index_id,
  product_dir,
  product_id,
  flg_r1,
  flg_r2,
  contamination_masks,
  sif_polygons_utm
) {
  message("Extracting ", index_id, " means for ", nrow(sif_polygons_utm), " SIF polygons.")

  sif_vect <- vect(sif_polygons_utm)

  if (index_id == "ndvi") {
    b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
    b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
    index_raster <- normalized_difference(b8, b4)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_10m)
    rm(b8, b4)
  } else if (index_id == "osavi") {
    b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
    b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
    denominator <- b8 + b4 + 0.16
    index_raster <- ifel(denominator != 0, 1.16 * (b8 - b4) / denominator, NA)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_10m)
    rm(b8, b4, denominator)
  } else if (index_id == "ndwi") {
    b3 <- load_scaled_land_band(product_dir, product_id, "B3", flg_r1, sif_vect)
    b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
    index_raster <- normalized_difference(b3, b8)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_10m)
    rm(b3, b8)
  } else if (index_id == "nirv") {
    b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
    b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
    ndvi <- normalized_difference(b8, b4)
    index_raster <- ndvi * b8
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_10m)
    rm(b8, b4, ndvi)
  } else if (index_id == "ndre") {
    b5 <- load_scaled_land_band(product_dir, product_id, "B5", flg_r2, sif_vect)
    b8 <- load_scaled_land_band(product_dir, product_id, "B8", flg_r1, sif_vect)
    b8_20m <- downsample_to_20m(b8, b5)
    index_raster <- normalized_difference(b8_20m, b5)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_20m)
    rm(b5, b8, b8_20m)
  } else if (index_id == "ndre8a") {
    b8a <- load_scaled_land_band(product_dir, product_id, "B8A", flg_r2, sif_vect)
    b6 <- load_scaled_land_band(product_dir, product_id, "B6", flg_r2, sif_vect)
    index_raster <- normalized_difference(b8a, b6)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_20m)
    rm(b8a, b6)
  } else if (index_id == "psri") {
    b6 <- load_scaled_land_band(product_dir, product_id, "B6", flg_r2, sif_vect)
    b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
    b2 <- load_scaled_land_band(product_dir, product_id, "B2", flg_r1, sif_vect)
    b4_20m <- downsample_to_20m(b4, b6)
    b2_20m <- downsample_to_20m(b2, b6)
    index_raster <- ifel(b6 != 0, (b4_20m - b2_20m) / b6, NA)
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_20m)
    rm(b6, b4, b2, b4_20m, b2_20m)
  } else if (index_id == "tcari") {
    b5 <- load_scaled_land_band(product_dir, product_id, "B5", flg_r2, sif_vect)
    b4 <- load_scaled_land_band(product_dir, product_id, "B4", flg_r1, sif_vect)
    b3 <- load_scaled_land_band(product_dir, product_id, "B3", flg_r1, sif_vect)
    b4_20m <- downsample_to_20m(b4, b5)
    b3_20m <- downsample_to_20m(b3, b5)
    index_raster <- ifel(
      b4_20m != 0,
      3 * ((b5 - b4_20m) - 0.2 * (b5 - b3_20m) * (b5 / b4_20m)),
      NA
    )
    index_raster <- apply_contamination_mask(index_raster, contamination_masks$mask_20m)
    rm(b5, b4, b3, b4_20m, b3_20m)
  } else {
    stop("Unsupported spectral index: ", index_id, call. = FALSE)
  }

  gc(verbose = FALSE)

  mean_values <- terra::extract(
    index_raster,
    sif_vect,
    fun = mean,
    na.rm = TRUE,
    touches = FALSE
  )

  value_col <- setdiff(names(mean_values), "ID")[1]

  index_summary <- mean_values %>%
    transmute(
      extract_id = ID,
      !!paste0("mean_", index_id) := .data[[value_col]]
    )

  rm(index_raster, mean_values, sif_vect)
  gc(verbose = FALSE)

  index_summary
}

wasp_products <- find_wasp_products(wasp_root, target_tile)

if (nrow(wasp_products) == 0) {
  stop("No local WASP products found under ", wasp_root, call. = FALSE)
}

sif_rows <- readRDS(check_file_exists(sif_file)) %>%
  st_drop_geometry() %>%
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
  target_crs <- crs(flg_r1)

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

  index_summaries <- target_indices %>%
    map(
      extract_index_summary,
      product_dir = product_dir,
      product_id = product_id,
      flg_r1 = flg_r1,
      flg_r2 = flg_r2,
      contamination_masks = contamination_masks,
      sif_polygons_utm = sif_polygons_utm
    ) %>%
    reduce(left_join, by = "extract_id")

  product_result <- sif_polygons %>%
    left_join(index_summaries, by = "extract_id") %>%
    mutate(
      wasp_product_id = product_id,
      wasp_year_month = year_month,
      non_crop_low_ndvi_threshold = non_crop_low_ndvi_threshold,
      contaminated_20m_fraction_threshold = contaminated_20m_fraction_threshold
    )

  rm(
    flg_r1, flg_r2, sif_polygons, sif_polygons_utm, sif_vect,
    contamination_masks, index_summaries
  )
  gc(verbose = FALSE)

  product_result
}

sif_wasp_spectral_indices_nonveg_masked <- pmap(
  wasp_products %>% select(product_id, product_dir, product_year, year_month),
  process_product
) %>%
  compact() %>%
  bind_rows() %>%
  select(-any_of(c("year_month", "extract_id")))

if (nrow(sif_wasp_spectral_indices_nonveg_masked) == 0) {
  stop("No SIF rows matched the local WASP products under ", wasp_root, call. = FALSE)
}

saveRDS(sif_wasp_spectral_indices_nonveg_masked, output_rds)

sif_wasp_spectral_indices_nonveg_masked %>%
  st_drop_geometry() %>%
  write_csv(output_csv)

message("Saved RDS to ", output_rds)
message("Saved CSV to ", output_csv)
