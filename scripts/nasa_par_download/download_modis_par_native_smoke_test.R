# Native-grid smoke test for MCD18A2.062 PAR at 12:00 UTC.
#
# This script randomly selects 10 distinct dates by shuffling valid rows from
# data/sif_dates.csv. It mosaics the downloaded MODIS tiles in their original
# sinusoidal CRS, then crops/masks them to the Germany bounding polygon.
# No raster reprojection or resampling is performed.
#
# It reuses helper functions from download_modis_par_sif_dates.R, but does not
# execute that script's full 201-date workflow.


# -----------------------------------------------------------------------------
# Load only the function/configuration section of the full downloader
# -----------------------------------------------------------------------------

helper_file <- "download_modis_par_sif_dates.R"
helper_lines <- readLines(helper_file, warn = FALSE)
main_heading <- grep("^# Main workflow$", helper_lines)

if (length(main_heading) != 1 || main_heading[[1]] < 2) {
  stop("Could not locate the Main workflow marker in: ", helper_file)
}

helper_code <- helper_lines[seq_len(main_heading[[1]] - 2)]
eval(parse(text = helper_code), envir = environment())


# -----------------------------------------------------------------------------
# Smoke-test configuration
# -----------------------------------------------------------------------------

smoke_date_count <- 10L
smoke_seed <- 20260712L

output_dir <- "data/modis_par_mcd18a2_germany_1200utc_native_smoke_test"
checklist_path <- file.path(
  output_dir,
  "mcd18a2_native_smoke_test_checklist.csv"
)
selection_path <- file.path(output_dir, "sampled_sif_rows_and_dates.csv")

keep_source_hdf <- FALSE
overwrite_outputs <- FALSE


# -----------------------------------------------------------------------------
# Select 10 distinct dates by randomly ordering valid SIF rows
# -----------------------------------------------------------------------------

sample_sif_rows_for_distinct_dates <- function(path, sample_count, seed) {
  input <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!(date_column %in% names(input))) {
    stop("Missing '", date_column, "' column in: ", path)
  }

  raw_dates <- trimws(as.character(input[[date_column]]))
  parsed_dates <- as.Date(substr(raw_dates, 1, 10), format = "%Y-%m-%d")
  valid_rows <- which(!is.na(parsed_dates))

  if (length(valid_rows) == 0) {
    stop("No valid dates found in: ", path)
  }

  set.seed(seed)
  shuffled_rows <- sample(valid_rows, length(valid_rows), replace = FALSE)
  first_row_for_each_date <- shuffled_rows[
    !duplicated(parsed_dates[shuffled_rows])
  ]

  if (length(first_row_for_each_date) < sample_count) {
    stop(
      "Requested ", sample_count, " distinct dates, but only ",
      length(first_row_for_each_date), " are available"
    )
  }

  selected_rows <- first_row_for_each_date[seq_len(sample_count)]

  data.frame(
    csv_data_row = selected_rows,
    csv_file_row = selected_rows + 1L,
    Delta_Date = parsed_dates[selected_rows],
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Native MODIS output paths and raster processing
# -----------------------------------------------------------------------------

output_paths <- function(target_date) {
  date_string <- format(target_date, "%Y-%m-%d")
  year_dir <- file.path(output_dir, format(target_date, "%Y"))
  prefix <- paste0(product_name, "_", date_string)

  list(
    par = file.path(
      year_dir,
      paste0(prefix, "_", par_layer, "_MODIS_Sinusoidal_native.tif")
    ),
    quality = file.path(
      year_dir,
      paste0(prefix, "_", quality_layer, "_MODIS_Sinusoidal_native.tif")
    )
  )
}

native_mosaic_layer <- function(hdf_paths, layer_name, germany_wgs84) {
  tiles <- lapply(
    hdf_paths,
    read_hdf_layer,
    layer_name = layer_name
  )

  if (layer_name == par_layer) {
    tiles <- lapply(tiles, valid_par_tile)
  } else if (layer_name == quality_layer) {
    tiles <- lapply(tiles, valid_quality_tile)
  }

  source_crs <- crs(tiles[[1]])
  if (is.na(source_crs) || !nzchar(source_crs)) {
    fatal_stop("The native MODIS tile has no CRS")
  }

  same_crs <- vapply(
    tiles,
    function(tile) same.crs(tile, tiles[[1]]),
    logical(1)
  )
  if (!all(same_crs)) {
    fatal_stop("Downloaded MODIS tiles do not share one native CRS")
  }

  # Native MODIS tiles already share resolution, origin, and sinusoidal CRS.
  # merge() places them together without reprojection or resampling.
  native_mosaic <- merge(
    sprc(tiles),
    first = TRUE,
    na.rm = TRUE,
    algo = 1
  )

  # Only the vector boundary is reprojected. Raster cells and values stay on
  # their original MODIS sinusoidal grid.
  germany_native <- project(germany_wgs84, crs(native_mosaic))
  native_subset <- crop(native_mosaic, germany_native, snap = "out")
  native_subset <- mask(native_subset, germany_native)
  names(native_subset) <- layer_name
  native_subset
}

process_native_date <- function(target_date, urls, token, germany_wgs84) {
  paths <- output_paths(target_date)
  date_string <- format(target_date, "%Y-%m-%d")

  if (keep_source_hdf) {
    download_dir <- file.path(
      output_dir,
      "raw_hdf",
      format(target_date, "%Y"),
      date_string
    )
    dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
    temporary_root <- NULL
  } else {
    temporary_root <- file.path(output_dir, "temporary_hdf")
    dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)
    download_dir <- tempfile(
      pattern = paste0(date_string, "_"),
      tmpdir = temporary_root
    )
    dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(
      safe_remove_temporary_directory(download_dir, temporary_root),
      add = TRUE
    )
  }

  hdf_paths <- vapply(
    urls,
    function(url) {
      destination <- file.path(download_dir, url_filename(url))
      download_hdf(url, destination, token)
    },
    character(1)
  )

  par <- native_mosaic_layer(
    hdf_paths,
    par_layer,
    germany_wgs84
  )
  quality <- native_mosaic_layer(
    hdf_paths,
    quality_layer,
    germany_wgs84
  )

  valid_pixels <- count_valid_pixels(par)
  if (is.na(valid_pixels) || valid_pixels == 0) {
    stop("No valid ", par_layer, " pixels found for ", date_string)
  }

  write_par_raster(par, paths$par)
  write_quality_raster(quality, paths$quality)

  list(
    par = paths$par,
    quality = paths$quality,
    valid_pixels = valid_pixels,
    resolution = res(par),
    crs = crs(par)
  )
}


# -----------------------------------------------------------------------------
# Run the smoke test
# -----------------------------------------------------------------------------

validate_bounds(germany_bounds)
token <- earthdata_token()
germany_wgs84 <- make_bounds_polygon(germany_bounds)

selection <- sample_sif_rows_for_distinct_dates(
  date_csv,
  smoke_date_count,
  smoke_seed
)
dates <- selection$Delta_Date

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(selection, selection_path, row.names = FALSE)
checklist <- load_checklist(checklist_path)

message("Native MODIS smoke-test dates: ", length(dates))
message("Random seed: ", smoke_seed)
message("Selected dates: ", paste(format(dates, "%Y-%m-%d"), collapse = ", "))
message("Output directory: ", output_dir)
message("No raster reprojection will be performed.")

for (index in seq_along(dates)) {
  target_date <- dates[index]
  date_string <- format(target_date, "%Y-%m-%d")
  paths <- output_paths(target_date)

  message("[", index, "/", length(dates), "] ", date_string)

  if (
    !overwrite_outputs &&
      file_is_complete(paths$par) &&
      file_is_complete(paths$quality)
  ) {
    message("  Existing native outputs are complete; skipped.")
    next
  }

  tryCatch(
    {
      urls <- search_granule_urls(target_date)
      granule_count <- length(urls)
      message("  Matching Germany tiles: ", granule_count)

      if (granule_count == 0) {
        checklist <- upsert_checklist(
          checklist,
          make_checklist_row(
            target_date,
            "no_granules",
            granules_found = 0,
            par_output = paths$par,
            quality_output = paths$quality,
            error = "No MCD18A2.062 granules intersected the date and bounds"
          )
        )
      } else {
        processed <- process_native_date(
          target_date,
          urls,
          token,
          germany_wgs84
        )

        checklist <- upsert_checklist(
          checklist,
          make_checklist_row(
            target_date,
            "completed",
            granules_found = granule_count,
            par_output = processed$par,
            quality_output = processed$quality,
            valid_par_pixels = processed$valid_pixels
          )
        )

        message("  Native resolution: ", paste(processed$resolution, collapse = " x "))
        message("  Wrote: ", processed$par)
        message("  Wrote: ", processed$quality)
      }

      save_checklist(checklist, checklist_path)
    },
    error = function(error) {
      checklist <<- upsert_checklist(
        checklist,
        make_checklist_row(
          target_date,
          "error",
          par_output = paths$par,
          quality_output = paths$quality,
          error = conditionMessage(error)
        )
      )
      save_checklist(checklist, checklist_path)
      message("  ERROR: ", conditionMessage(error))
    }
  )
}

completed <- sum(checklist$status == "completed", na.rm = TRUE)
message("Finished native smoke test: ", completed, "/", length(dates))
message("Sample manifest: ", selection_path)
message("Checklist: ", checklist_path)
