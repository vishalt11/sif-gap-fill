# Convert downloaded GLASS HDF files to GeoTIFFs for rasterio-based Python work.
#
# Output structure mirrors the input structure:
#   data/glass_fapar_modis_250m/h18v03/2019/file.hdf
# becomes:
#   data/glass_geotiff/fapar/h18v03/2019/file.tif

library(terra)

#-------------------------------------------------------------------------------
# Config

OUTPUT_ROOT <- "data/glass_geotiff"

OVERWRITE <- FALSE

# Keep TRUE if you want the GeoTIFFs to already have obvious fill/invalid values
# masked out. Set FALSE if you want a direct HDF -> GeoTIFF copy.
APPLY_VALUE_MASKS <- TRUE

PRODUCTS <- list(
  fapar = list(
    input_dir = "data/glass_fapar_modis_250m",
    output_dir = file.path(OUTPUT_ROOT, "fapar"),
    valid_min = 0,
    valid_max = 1
  ),
  evi = list(
    input_dir = "data/glass_evi_modis_250m",
    output_dir = file.path(OUTPUT_ROOT, "evi"),
    valid_min = -1,
    valid_max = 1
  ),
  ndvi = list(
    input_dir = "data/glass_ndvi_modis_250m",
    output_dir = file.path(OUTPUT_ROOT, "ndvi"),
    valid_min = -1,
    valid_max = 1
  ),
  par = list(
    input_dir = "data/glass_par_modis_005d",
    output_dir = file.path(OUTPUT_ROOT, "par"),
    valid_min = 0,
    valid_max = Inf
  )
)

# Edit this if you only want to convert one product first, e.g. c("fapar").
PRODUCTS_TO_CONVERT <- names(PRODUCTS)

GDAL_OPTIONS <- c(
  "COMPRESS=LZW",
  "TILED=YES",
  "BIGTIFF=IF_SAFER"
)


#-------------------------------------------------------------------------------
# Helpers

relative_path <- function(path, root) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  substring(path_norm, nchar(root_norm) + 2)
}

output_path_for_hdf <- function(hdf_file, input_dir, output_dir) {
  rel <- relative_path(hdf_file, input_dir)
  rel_tif <- sub("\\.hdf$", ".tif", rel, ignore.case = TRUE)
  file.path(output_dir, rel_tif)
}

mask_invalid_values <- function(r, valid_min, valid_max) {
  if (!APPLY_VALUE_MASKS) {
    return(r)
  }

  if (is.finite(valid_max)) {
    terra::ifel(r < valid_min | r > valid_max, NA, r)
  } else {
    terra::ifel(r < valid_min, NA, r)
  }
}

convert_one_hdf <- function(hdf_file, product_name, product_cfg, index, total) {
  tif_file <- output_path_for_hdf(
    hdf_file = hdf_file,
    input_dir = product_cfg$input_dir,
    output_dir = product_cfg$output_dir
  )

  if (file.exists(tif_file) && !OVERWRITE) {
    message(sprintf("[%s %d/%d] exists, skipping: %s", product_name, index, total, tif_file))
    return(invisible(tif_file))
  }

  dir.create(dirname(tif_file), recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[%s %d/%d] converting: %s", product_name, index, total, hdf_file))

  tryCatch(
    {
      r <- terra::rast(hdf_file)
      r <- mask_invalid_values(
        r,
        valid_min = product_cfg$valid_min,
        valid_max = product_cfg$valid_max
      )

      terra::writeRaster(
        r,
        tif_file,
        overwrite = OVERWRITE,
        NAflag = -9999,
        gdal = GDAL_OPTIONS
      )

      invisible(tif_file)
    },
    error = function(e) {
      message(sprintf("  FAILED: %s", conditionMessage(e)))
      invisible(NA_character_)
    }
  )
}

convert_product <- function(product_name) {
  product_cfg <- PRODUCTS[[product_name]]

  if (is.null(product_cfg)) {
    stop("Unknown product: ", product_name)
  }

  if (!dir.exists(product_cfg$input_dir)) {
    message(sprintf("[%s] input directory does not exist, skipping: %s", product_name, product_cfg$input_dir))
    return(invisible(NULL))
  }

  hdf_files <- list.files(
    product_cfg$input_dir,
    pattern = "\\.hdf$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(hdf_files) == 0) {
    message(sprintf("[%s] no HDF files found in %s", product_name, product_cfg$input_dir))
    return(invisible(NULL))
  }

  message(sprintf("[%s] found %d HDF files", product_name, length(hdf_files)))

  for (i in seq_along(hdf_files)) {
    convert_one_hdf(
      hdf_file = hdf_files[[i]],
      product_name = product_name,
      product_cfg = product_cfg,
      index = i,
      total = length(hdf_files)
    )
  }

  invisible(NULL)
}


#-------------------------------------------------------------------------------
# Run

for (product_name in PRODUCTS_TO_CONVERT) {
  convert_product(product_name)
}

message("Done converting GLASS HDF files to GeoTIFF.")
