library(tidyverse)
library(sf)
library(terra)

state_prefix <- "ba"
sif_file <- file.path("data", paste0(state_prefix, "_sif_mgrs_crop_composition.rds"))
wasp_tile_root <- "data/geodes_wasp_zips"

target_tile <- "32UNA"
target_bands_10m <- c("B2", "B3", "B4", "B8")
target_bands_20m <- c("B5", "B6", "B7", "B8A", "B11", "B12")
target_bands <- c(target_bands_10m, target_bands_20m)

output_rds <- file.path(
  "data",
  paste0(state_prefix, "_sif_", target_tile, "_wasp_band_means.rds")
)
output_csv <- file.path(
  "data",
  paste0(state_prefix, "_sif_", target_tile, "_wasp_band_means.csv")
)

reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4

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

band_resolution_group <- function(band_id) {
  case_when(
    band_id %in% target_bands_10m ~ "R1",
    band_id %in% target_bands_20m ~ "R2",
    TRUE ~ NA_character_
  )
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

extract_band_summary <- function(band_id, product_dir, product_id, flg_r1, flg_r2, sif_polygons_utm) {
  message("Extracting ", band_id, " means for ", nrow(sif_polygons_utm), " SIF polygons.")

  band_group <- band_resolution_group(band_id)

  if (is.na(band_group)) {
    stop("Unsupported band: ", band_id, call. = FALSE)
  }

  flg <- if (band_group == "R1") flg_r1 else flg_r2

  band <- rast(check_file_exists(wasp_file(product_dir, product_id, "FRC", band_id)))
  band <- crop(band, vect(sif_polygons_utm))

  flg_crop <- crop(flg, band)

  band <- ifel(
    flg_crop == land_flag_value & band != reflectance_nodata,
    band / reflectance_quantification_value,
    NA
  )

  mean_values <- terra::extract(
    band,
    vect(sif_polygons_utm),
    fun = mean,
    na.rm = TRUE,
    touches = FALSE
  )

  value_col <- setdiff(names(mean_values), "ID")[1]
  band_lower <- str_to_lower(band_id)

  band_summary <- mean_values %>%
    transmute(
      extract_id = ID,
      !!paste0("mean_", band_lower) := .data[[value_col]]
    )

  # To audit extraction quality later, re-enable a valid-pixel-count extract here.
  rm(band, flg_crop, mean_values)
  gc(verbose = FALSE)

  band_summary
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

process_product <- function(product_id, product_dir, year_month) {
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

  band_summaries <- target_bands %>%
    map(
      extract_band_summary,
      product_dir = product_dir,
      product_id = product_id,
      flg_r1 = flg_r1,
      flg_r2 = flg_r2,
      sif_polygons_utm = sif_polygons_utm
    ) %>%
    reduce(left_join, by = "extract_id")

  product_result <- sif_polygons %>%
    left_join(band_summaries, by = "extract_id") %>%
    mutate(
      wasp_product_id = product_id,
      wasp_year_month = year_month
    )

  rm(flg_r1, flg_r2, sif_polygons, sif_polygons_utm, band_summaries)
  gc(verbose = FALSE)

  product_result
}

sif_wasp_band_means <- pmap(
  wasp_products %>% select(product_id, product_dir, year_month),
  process_product
) %>%
  compact() %>%
  bind_rows() %>%
  select(-any_of(c("year_month", "extract_id")))

if (nrow(sif_wasp_band_means) == 0) {
  stop("No SIF rows matched the local WASP products under ", wasp_root, call. = FALSE)
}

saveRDS(sif_wasp_band_means, output_rds)

sif_wasp_band_means %>%
  st_drop_geometry() %>%
  write_csv(output_csv)

message("Saved RDS to ", output_rds)
message("Saved CSV to ", output_csv)
