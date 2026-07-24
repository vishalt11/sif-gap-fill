# Download native-grid VIIRS VNP18A2.002 daily-mean PAR for all SIF dates.
#
# The script reads unique Delta_Date values from data/sif_dates.csv, downloads
# the VNP18A2 Version 2 tiles intersecting Germany, extracts Daily_Mean_PAR and
# PAR_Quality, mosaics the tiles in their original VIIRS sinusoidal grid, and
# crops/masks them with the Germany bounding polygon.
#
# No raster reprojection or resampling is performed. Only the Germany vector
# boundary is transformed into the native raster CRS for crop/mask operations.
#
# Before running, expose a NASA Earthdata user token to R:
#
#   Sys.setenv(EARTHDATA_TOKEN = "your-token")
#
# Required packages: terra, httr2


# -----------------------------------------------------------------------------
# Load only reusable helpers from the existing MODIS R downloader
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
# VNP18A2 Version 2 configuration
# -----------------------------------------------------------------------------

date_csv <- "data/main_sif_data/sif_dates_11tiles.csv"
date_column <- "Delta_Date"

collection_concept_id <- "C2631841566-LPCLOUD"
product_name <- "VNP18A2.002"
par_layer <- "Daily_Mean_PAR"
quality_layer <- "PAR_Quality"

# Padded German national extent, including offshore islands.
# Order: west, south, east, north in WGS84 longitude/latitude.
germany_bounds <- c(5.5, 47.0, 15.5, 55.2)

# The padded Germany bounds intersect these four native VIIRS sinusoidal tiles.
# A date is considered complete only after exactly one latest-production
# granule has been selected for each tile.
expected_germany_tiles <- c("h18v03", "h18v04", "h19v03", "h19v04")

output_dir <- "data/viirs_vnp18a2_daily_mean_par_germany_native"
checklist_path <- file.path(
  output_dir,
  "vnp18a2_daily_mean_par_download_checklist.csv"
)

keep_source_hdf <- FALSE
overwrite_outputs <- FALSE
request_timeout_seconds <- 300
request_retries <- 4


# -----------------------------------------------------------------------------
# Override MODIS-specific URL selection for VIIRS HDF-EOS5 (.h5) granules
# -----------------------------------------------------------------------------

cmr_entry_data_url <- function(entry) {
  links <- entry$links %||% list()
  if (length(links) == 0) {
    return(NA_character_)
  }

  candidates <- vapply(
    links,
    function(link) {
      href <- link$href %||% ""
      rel <- link$rel %||% ""
      inherited_value <- link$inherited %||% FALSE
      inherited <- isTRUE(inherited_value) ||
        identical(tolower(as.character(inherited_value)), "true")
      clean_href <- sub("\\?.*$", "", href)

      is_h5 <- grepl("\\.h5$", clean_href, ignore.case = TRUE)
      is_data <- grepl("data#", rel, fixed = TRUE)
      is_https <- startsWith(href, "https://")

      if (!inherited && is_h5 && is_data && is_https) href else NA_character_
    },
    character(1)
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  lpdaac_cloud <- grepl(
    "data.lpdaac.earthdatacloud.nasa.gov",
    candidates,
    fixed = TRUE
  )
  if (any(lpdaac_cloud)) {
    return(candidates[which(lpdaac_cloud)[[1]]])
  }

  candidates[[1]]
}

url_filename <- function(url) {
  clean_url <- sub("\\?.*$", "", url)
  filename <- basename(clean_url)

  if (!nzchar(filename) || !grepl("\\.h5$", filename, ignore.case = TRUE)) {
    stop("Could not determine VIIRS HDF5 filename from URL: ", url)
  }

  filename
}

vnp18a2_url_metadata <- function(url) {
  filename <- url_filename(url)
  pattern <- paste0(
    "^VNP18A2\\.A([0-9]{7})\\.",
    "(h[0-9]{2}v[0-9]{2})\\.002\\.",
    "([0-9]{13})\\.h5$"
  )
  match <- regmatches(filename, regexec(pattern, filename))[[1]]

  if (length(match) == 0) {
    stop("Unexpected VNP18A2.002 granule filename: ", filename)
  }

  data.frame(
    url = url,
    filename = filename,
    acquisition = paste0("A", match[[2]]),
    tile = match[[3]],
    production_time = match[[4]],
    stringsAsFactors = FALSE
  )
}

target_acquisition_token <- function(target_date) {
  paste0("A", format(target_date, "%Y"), format(target_date, "%j"))
}

validate_vnp18a2_filename_parser <- function() {
  test_filename <- "VNP18A2.A2022086.h18v03.002.2025209191452.h5"
  test_url <- paste0("https://example.invalid/", test_filename)
  parsed <- vnp18a2_url_metadata(test_url)

  checks <- c(
    identical(parsed$filename[[1]], test_filename),
    identical(parsed$acquisition[[1]], "A2022086"),
    identical(parsed$tile[[1]], "h18v03"),
    identical(parsed$production_time[[1]], "2025209191452"),
    identical(
      parsed$acquisition[[1]],
      target_acquisition_token(as.Date("2022-03-27"))
    )
  )

  if (!all(checks)) {
    stop(
      "Internal VNP18A2 filename-parser self-test failed. ",
      "No CMR search or raster download was attempted."
    )
  }

  invisible(TRUE)
}

# Override the generic MODIS CMR search. CMR can retain multiple reprocessed
# VNP18A2 files for the same acquisition date and tile. This implementation
# first requires the filename acquisition date to equal the requested date,
# then retains the file with the latest 13-digit production timestamp per tile.
search_granule_urls <- function(target_date) {
  date_string <- format(target_date, "%Y-%m-%d")
  temporal <- paste0(
    date_string, "T00:00:00Z,",
    date_string, "T23:59:59Z"
  )

  cmr_request <- request(cmr_granules_url) |>
    req_url_query(
      concept_id = collection_concept_id,
      temporal = temporal,
      bounding_box = paste(germany_bounds, collapse = ","),
      downloadable = "true",
      page_size = 200
    ) |>
    req_headers(Accept = "application/json") |>
    req_user_agent("VNP18A2-Germany-daily-mean-PAR-R-downloader") |>
    req_timeout(request_timeout_seconds) |>
    req_retry(max_tries = request_retries)

  response <- req_perform(cmr_request)
  body <- resp_body_json(response, simplifyVector = FALSE)
  entries <- body$feed$entry %||% list()

  raw_urls <- vapply(entries, cmr_entry_data_url, character(1))
  raw_urls <- unique(raw_urls[!is.na(raw_urls) & nzchar(raw_urls)])

  if (length(raw_urls) == 0) {
    message("  CMR candidates: 0")
    return(character(0))
  }

  metadata <- do.call(
    rbind,
    lapply(raw_urls, vnp18a2_url_metadata)
  )
  expected_acquisition <- target_acquisition_token(target_date)
  date_matched <- metadata[
    metadata$acquisition == expected_acquisition,
    ,
    drop = FALSE
  ]

  message(
    "  CMR candidates: ", nrow(metadata),
    "; requested-date matches: ", nrow(date_matched)
  )

  if (nrow(date_matched) == 0) {
    stop(
      "CMR returned granules, but none has acquisition token ",
      expected_acquisition,
      " for ", date_string
    )
  }

  unexpected_tiles <- setdiff(unique(date_matched$tile), expected_germany_tiles)
  if (length(unexpected_tiles) > 0) {
    stop(
      "Unexpected VNP18A2 tile(s) returned for Germany: ",
      paste(unexpected_tiles, collapse = ", ")
    )
  }

  # A 13-digit YYYYDDDHHMMSS timestamp is safely representable as an exact
  # integer in an R double. Sort newest first within each spatial tile.
  date_matched$production_number <- as.numeric(date_matched$production_time)
  date_matched <- date_matched[
    order(date_matched$tile, -date_matched$production_number),
    ,
    drop = FALSE
  ]
  selected <- date_matched[!duplicated(date_matched$tile), , drop = FALSE]

  missing_tiles <- setdiff(expected_germany_tiles, selected$tile)
  if (length(missing_tiles) > 0) {
    stop(
      "Missing expected VNP18A2 Germany tile(s) for ", date_string, ": ",
      paste(missing_tiles, collapse = ", ")
    )
  }

  selected <- selected[
    match(expected_germany_tiles, selected$tile),
    ,
    drop = FALSE
  ]

  duplicate_count <- nrow(date_matched) - nrow(selected)
  if (duplicate_count > 0) {
    message(
      "  Removed ", duplicate_count,
      " older reprocessed duplicate(s); retained latest production per tile:"
    )
  } else {
    message("  One requested-date granule was available per tile:")
  }
  for (index in seq_len(nrow(selected))) {
    message(
      "    ", selected$tile[[index]], " -> ",
      selected$production_time[[index]], "  ",
      selected$filename[[index]]
    )
  }

  selected$url
}


# -----------------------------------------------------------------------------
# Native VIIRS output paths and raster processing
# -----------------------------------------------------------------------------

output_paths <- function(target_date) {
  date_string <- format(target_date, "%Y-%m-%d")
  year_dir <- file.path(output_dir, format(target_date, "%Y"))
  prefix <- paste0(product_name, "_", date_string)

  list(
    par = file.path(
      year_dir,
      paste0(prefix, "_", par_layer, "_VIIRS_Sinusoidal_native.tif")
    ),
    quality = file.path(
      year_dir,
      paste0(prefix, "_", quality_layer, "_VIIRS_Sinusoidal_native.tif")
    )
  )
}

# Preserve every documented raw quality code from 0 through 4. The quality
# raster is categorical and is never interpolated.
valid_quality_tile <- function(raster) {
  ifel(raster >= 0 & raster <= 4, raster, NA)
}

native_mosaic_layer <- function(hdf_paths, layer_name, germany_wgs84) {
  tiles <- lapply(
    hdf_paths,
    read_hdf_layer,
    layer_name = layer_name
  )

  if (length(tiles) == 0) {
    stop("No VIIRS tiles were available for layer: ", layer_name)
  }

  if (layer_name == par_layer) {
    tiles <- lapply(tiles, valid_par_tile)
  } else if (layer_name == quality_layer) {
    tiles <- lapply(tiles, valid_quality_tile)
  }

  source_crs <- crs(tiles[[1]])
  if (is.na(source_crs) || !nzchar(source_crs)) {
    fatal_stop("The native VIIRS tile has no CRS")
  }

  same_native_crs <- vapply(
    tiles,
    function(tile) same.crs(tile, tiles[[1]]),
    logical(1)
  )
  if (!all(same_native_crs)) {
    fatal_stop("Downloaded VNP18A2 tiles do not share one native CRS")
  }

  # Adjacent VNP18A2 tiles already share one sinusoidal grid. merge() places
  # them together without changing their cell resolution, origin, or values.
  native_mosaic <- merge(
    sprc(tiles),
    first = TRUE,
    na.rm = TRUE,
    algo = 1
  )

  # Transform only the vector boundary; the raster remains untouched.
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
      "raw_h5",
      format(target_date, "%Y"),
      date_string
    )
    dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
    temporary_root <- NULL
  } else {
    temporary_root <- file.path(output_dir, "temporary_h5")
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

  daily_mean_par <- native_mosaic_layer(
    hdf_paths,
    par_layer,
    germany_wgs84
  )
  quality <- native_mosaic_layer(
    hdf_paths,
    quality_layer,
    germany_wgs84
  )

  valid_pixels <- count_valid_pixels(daily_mean_par)
  if (is.na(valid_pixels) || valid_pixels == 0) {
    stop("No valid ", par_layer, " pixels found for ", date_string)
  }

  write_par_raster(daily_mean_par, paths$par)
  write_quality_raster(quality, paths$quality)

  list(
    par = paths$par,
    quality = paths$quality,
    valid_pixels = valid_pixels,
    resolution = res(daily_mean_par),
    crs = crs(daily_mean_par)
  )
}


# -----------------------------------------------------------------------------
# Main workflow: every unique SIF date
# -----------------------------------------------------------------------------

validate_vnp18a2_filename_parser()
validate_bounds(germany_bounds)
dates <- load_unique_dates(date_csv)
token <- earthdata_token()
germany_wgs84 <- make_bounds_polygon(germany_bounds)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checklist <- load_checklist(checklist_path)

message("Unique SIF dates: ", length(dates))
message(
  "First/last date: ", format(dates[[1]], "%Y-%m-%d"), " / ",
  format(dates[[length(dates)]], "%Y-%m-%d")
)
message("Earthdata product/layer: ", product_name, " / ", par_layer)
message("Germany WGS84 bounds: ", paste(germany_bounds, collapse = ", "))
message("Output directory: ", output_dir)
message("The VIIRS raster grid will remain in its native sinusoidal CRS.")

for (index in seq_along(dates)) {
  target_date <- dates[index]
  date_string <- format(target_date, "%Y-%m-%d")
  paths <- output_paths(target_date)

  message("[", index, "/", length(dates), "] ", date_string)

  previous <- checklist[checklist$date == date_string, , drop = FALSE]
  previous_is_verified <-
    nrow(previous) == 1 &&
    identical(previous$status[[1]], "completed") &&
    identical(
      suppressWarnings(as.integer(previous$granules_found[[1]])),
      length(expected_germany_tiles)
    )
  outputs_exist <-
    file_is_complete(paths$par) &&
    file_is_complete(paths$quality)

  if (!overwrite_outputs && outputs_exist && previous_is_verified) {
    message(
      "  Existing outputs have a completed four-tile checklist record; skipped."
    )
    next
  }

  if (!overwrite_outputs && outputs_exist && !previous_is_verified) {
    previous_count <- if (nrow(previous) == 1) {
      previous$granules_found[[1]]
    } else {
      "missing"
    }
    message(
      "  Existing outputs are legacy/unverified (checklist granules: ",
      previous_count,
      "); regenerating with latest-production granules."
    )
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
            error = "No VNP18A2.002 granules intersected the date and bounds"
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

        message(
          "  Native resolution: ",
          paste(processed$resolution, collapse = " x ")
        )
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
message(
  "Finished. Completed dates in checklist: ",
  completed, "/", length(dates)
)
message("Checklist: ", checklist_path)
