library(tidyverse)
library(arrow)

area_to_point_dir <- file.path(
  "data/area_to_point_nonveg_masked",
  "ba_sif_32UQV_wasp_area_to_point_20m"
)

parquet_dir <- file.path(area_to_point_dir, "parquet")
dir.create(parquet_dir, recursive = TRUE, showWarnings = FALSE)

check_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    stop("Missing directory: ", path, call. = FALSE)
  }

  path
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, call. = FALSE)
  }

  path
}

safe_parquet_name <- function(path) {
  paste0(tools::file_path_sans_ext(basename(path)), ".parquet")
}

convert_rds_to_parquet <- function(input_rds, output_parquet) {
  message("Converting ", input_rds)

  object <- readRDS(input_rds)

  if (inherits(object, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to drop geometry from ", input_rds, call. = FALSE)
    }

    object <- sf::st_drop_geometry(object)
  }

  object <- as_tibble(object)
  write_parquet(object, output_parquet)

  tibble(
    input_rds = input_rds,
    output_parquet = output_parquet,
    n_rows = nrow(object),
    n_cols = ncol(object)
  )
}

check_dir_exists(area_to_point_dir)

pixel_rds_files <- list.files(
  area_to_point_dir,
  pattern = "_area_to_point_pixels_20m\\.rds$",
  full.names = TRUE
)

if (length(pixel_rds_files) == 0) {
  stop("No pixel RDS files found in ", area_to_point_dir, call. = FALSE)
}

pixel_conversion_manifest <- map_dfr(
  pixel_rds_files,
  ~ convert_rds_to_parquet(
    input_rds = .x,
    output_parquet = file.path(parquet_dir, safe_parquet_name(.x))
  )
)

polygon_targets_rds <- check_file_exists(file.path(area_to_point_dir, "polygon_targets.rds"))
polygon_targets_parquet <- file.path(parquet_dir, "polygon_targets.parquet")
polygon_conversion <- convert_rds_to_parquet(polygon_targets_rds, polygon_targets_parquet)

pixel_table_manifest_csv <- file.path(area_to_point_dir, "pixel_table_manifest.csv")

if (file.exists(pixel_table_manifest_csv)) {
  pixel_table_manifest <- read_csv(pixel_table_manifest_csv, show_col_types = FALSE) %>%
    mutate(
      pixel_parquet = file.path(
        parquet_dir,
        safe_parquet_name(pixel_rds)
      )
    )

  write_parquet(
    pixel_table_manifest,
    file.path(parquet_dir, "pixel_table_manifest.parquet")
  )

  write_csv(
    pixel_table_manifest,
    file.path(parquet_dir, "pixel_table_manifest.csv")
  )
}

conversion_manifest <- bind_rows(
  pixel_conversion_manifest,
  polygon_conversion
)

write_csv(
  conversion_manifest,
  file.path(parquet_dir, "rds_to_parquet_conversion_manifest.csv")
)

message("Saved parquet files to ", parquet_dir)
message("Saved conversion manifest to ", file.path(parquet_dir, "rds_to_parquet_conversion_manifest.csv"))
