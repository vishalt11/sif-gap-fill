# Download MCD18A2.062 PAR at 12:00 UTC for the unique dates in sif_dates.csv.
#
# Native R implementation using terra for HDF-EOS2 reading and raster work.
# It does not use Python, Rasterio, or reticulate.
#
# Before running, create a NASA Earthdata user token and make it available to R:
#
#   Sys.setenv(EARTHDATA_TOKEN = "your-token")
#
# Do not paste the token into this script. A persistent token can instead be
# placed in the user's .Renviron file as EARTHDATA_TOKEN=your-token.
#
# Required R packages: terra, httr2

library(terra)
library(httr2)


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

date_csv <- "data/main_sif_data/sif_dates.csv"
date_column <- "Delta_Date"
output_dir <- "data/modis_par_mcd18a2_germany_1200utc"

collection_concept_id <- "C2486282714-LPCLOUD"
product_name <- "MCD18A2.062"
par_layer <- "GMT_1200_PAR"
quality_layer <- "PAR_Quality"

# Padded German national extent, including offshore islands.
# Order: west, south, east, north in WGS84 longitude/latitude.
germany_bounds <- c(5.5, 47.0, 15.5, 55.2)

output_crs <- "EPSG:3035"
output_resolution_m <- 1000

keep_source_hdf <- FALSE
overwrite_outputs <- FALSE
request_timeout_seconds <- 300
request_retries <- 4

cmr_granules_url <- "https://cmr.earthdata.nasa.gov/search/granules.json"
checklist_path <- file.path(
  output_dir,
  "mcd18a2_1200utc_download_checklist.csv"
)

checklist_columns <- c(
  "date",
  "status",
  "granules_found",
  "par_output",
  "quality_output",
  "valid_par_pixels",
  "error",
  "updated_at_utc"
)


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0) fallback else x
}

utc_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

fatal_stop <- function(...) {
  condition <- simpleError(paste0(...))
  class(condition) <- c("fatal_processing_error", class(condition))
  stop(condition)
}

file_is_complete <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

validate_bounds <- function(bounds) {
  if (length(bounds) != 4 || any(!is.finite(bounds))) {
    stop("germany_bounds must contain four finite numbers")
  }

  west <- bounds[[1]]
  south <- bounds[[2]]
  east <- bounds[[3]]
  north <- bounds[[4]]

  if (west < -180 || east > 180 || west >= east) {
    stop("Invalid west/east Germany bounds")
  }
  if (south < -90 || north > 90 || south >= north) {
    stop("Invalid south/north Germany bounds")
  }
}

load_unique_dates <- function(path) {
  if (!file.exists(path)) {
    stop("Date CSV does not exist: ", path)
  }

  dates_df <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!(date_column %in% names(dates_df))) {
    stop(
      "Date CSV must contain column '", date_column,
      "'. Found: ", paste(names(dates_df), collapse = ", ")
    )
  }

  raw_dates <- trimws(as.character(dates_df[[date_column]]))
  raw_dates <- raw_dates[nzchar(raw_dates)]
  parsed_dates <- as.Date(substr(raw_dates, 1, 10), format = "%Y-%m-%d")

  invalid <- is.na(parsed_dates)
  if (any(invalid)) {
    stop(
      "Invalid ", date_column, " values: ",
      paste(unique(raw_dates[invalid]), collapse = ", ")
    )
  }

  dates <- sort(unique(parsed_dates))
  if (length(dates) == 0) {
    stop("No valid dates found in: ", path)
  }

  dates
}

make_bounds_polygon <- function(bounds) {
  west <- bounds[[1]]
  south <- bounds[[2]]
  east <- bounds[[3]]
  north <- bounds[[4]]

  coordinates <- matrix(
    c(
      west, south,
      east, south,
      east, north,
      west, north,
      west, south
    ),
    ncol = 2,
    byrow = TRUE
  )

  vect(coordinates, type = "polygons", crs = "EPSG:4326")
}

make_target_grid <- function(bounds_polygon_wgs84) {
  bounds_polygon_projected <- project(bounds_polygon_wgs84, output_crs)
  projected_extent <- ext(bounds_polygon_projected)

  left <- floor(xmin(projected_extent) / output_resolution_m) * output_resolution_m
  right <- ceiling(xmax(projected_extent) / output_resolution_m) * output_resolution_m
  bottom <- floor(ymin(projected_extent) / output_resolution_m) * output_resolution_m
  top <- ceiling(ymax(projected_extent) / output_resolution_m) * output_resolution_m

  width <- as.integer(round((right - left) / output_resolution_m))
  height <- as.integer(round((top - bottom) / output_resolution_m))

  template <- rast(
    ncols = width,
    nrows = height,
    xmin = left,
    xmax = right,
    ymin = bottom,
    ymax = top,
    crs = output_crs
  )

  list(
    template = template,
    mask = bounds_polygon_projected
  )
}


# -----------------------------------------------------------------------------
# Checklist helpers
# -----------------------------------------------------------------------------

empty_checklist <- function() {
  result <- as.data.frame(
    setNames(
      replicate(length(checklist_columns), character(), simplify = FALSE),
      checklist_columns
    ),
    stringsAsFactors = FALSE
  )
  result
}

load_checklist <- function(path) {
  if (!file.exists(path)) {
    return(empty_checklist())
  }

  checklist <- read.csv(
    path,
    stringsAsFactors = FALSE,
    colClasses = "character",
    check.names = FALSE
  )

  missing_columns <- setdiff(checklist_columns, names(checklist))
  for (column in missing_columns) {
    checklist[[column]] <- ""
  }

  checklist[, checklist_columns, drop = FALSE]
}

make_checklist_row <- function(
    target_date,
    status,
    granules_found = "",
    par_output = "",
    quality_output = "",
    valid_par_pixels = "",
    error = "") {
  data.frame(
    date = format(target_date, "%Y-%m-%d"),
    status = as.character(status),
    granules_found = as.character(granules_found),
    par_output = as.character(par_output),
    quality_output = as.character(quality_output),
    valid_par_pixels = as.character(valid_par_pixels),
    error = as.character(error),
    updated_at_utc = utc_now(),
    stringsAsFactors = FALSE
  )
}

upsert_checklist <- function(checklist, row) {
  checklist <- checklist[checklist$date != row$date[[1]], , drop = FALSE]
  checklist <- rbind(checklist, row)
  checklist[order(checklist$date), , drop = FALSE]
}

save_checklist <- function(checklist, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- paste0(path, ".part")
  write.csv(checklist, temporary_path, row.names = FALSE, na = "")

  if (file.exists(path)) {
    file.remove(path)
  }
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically replace checklist: ", path)
  }
}

output_paths <- function(target_date) {
  date_string <- format(target_date, "%Y-%m-%d")
  year_dir <- file.path(output_dir, format(target_date, "%Y"))
  prefix <- paste0(product_name, "_", date_string)

  list(
    par = file.path(
      year_dir,
      paste0(prefix, "_", par_layer, "_EPSG3035_1km.tif")
    ),
    quality = file.path(
      year_dir,
      paste0(prefix, "_", quality_layer, "_EPSG3035_1km.tif")
    )
  )
}


# -----------------------------------------------------------------------------
# NASA CMR search and Earthdata download
# -----------------------------------------------------------------------------

earthdata_token <- function() {
  token <- trimws(Sys.getenv("EARTHDATA_TOKEN", unset = ""))
  if (!nzchar(token)) {
    stop(
      "EARTHDATA_TOKEN is not set. Generate a NASA Earthdata user token, ",
      "then run Sys.setenv(EARTHDATA_TOKEN = 'your-token') before this script."
    )
  }
  token
}

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

      is_hdf <- grepl("\\.hdf$", clean_href, ignore.case = TRUE)
      is_data <- grepl("data#", rel, fixed = TRUE)
      is_https <- startsWith(href, "https://")

      if (!inherited && is_hdf && is_data && is_https) href else NA_character_
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

search_granule_urls <- function(target_date) {
  date_string <- format(target_date, "%Y-%m-%d")
  temporal <- paste0(
    date_string, "T00:00:00Z,",
    date_string, "T23:59:59Z"
  )

  request <- request(cmr_granules_url) |>
    req_url_query(
      concept_id = collection_concept_id,
      temporal = temporal,
      bounding_box = paste(germany_bounds, collapse = ","),
      downloadable = "true",
      page_size = 20
    ) |>
    req_headers(Accept = "application/json") |>
    req_user_agent("MCD18A2-Germany-PAR-R-downloader") |>
    req_timeout(request_timeout_seconds) |>
    req_retry(max_tries = request_retries)

  response <- req_perform(request)
  body <- resp_body_json(response, simplifyVector = FALSE)
  entries <- body$feed$entry %||% list()

  urls <- vapply(entries, cmr_entry_data_url, character(1))
  unique(urls[!is.na(urls) & nzchar(urls)])
}

url_filename <- function(url) {
  clean_url <- sub("\\?.*$", "", url)
  filename <- basename(clean_url)
  if (!nzchar(filename) || !grepl("\\.hdf$", filename, ignore.case = TRUE)) {
    stop("Could not determine HDF filename from URL: ", url)
  }
  filename
}

download_hdf <- function(url, destination, token) {
  if (file_is_complete(destination)) {
    return(destination)
  }

  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- paste0(destination, ".part")

  request <- request(url) |>
    req_headers(Authorization = paste("Bearer", token)) |>
    req_user_agent("MCD18A2-Germany-PAR-R-downloader") |>
    req_timeout(request_timeout_seconds) |>
    req_retry(max_tries = request_retries)

  req_perform(request, path = temporary_path)

  if (!file_is_complete(temporary_path)) {
    stop("Downloaded file is empty: ", url)
  }

  if (file.exists(destination)) {
    file.remove(destination)
  }
  if (!file.rename(temporary_path, destination)) {
    stop("Could not move completed download to: ", destination)
  }

  destination
}


# -----------------------------------------------------------------------------
# HDF subdataset extraction and raster processing with terra
# -----------------------------------------------------------------------------

normalize_layer_name <- function(x) {
  gsub("[[:space:]\"]", "", x)
}

read_hdf_layer <- function(hdf_path, layer_name) {
  direct <- try(
    rast(hdf_path, subds = layer_name),
    silent = TRUE
  )

  if (!inherits(direct, "try-error") && inherits(direct, "SpatRaster")) {
    if (nlyr(direct) >= 1) {
      return(direct[[1]])
    }
  }

  datasets <- try(sds(hdf_path), silent = TRUE)
  if (inherits(datasets, "try-error")) {
    fatal_stop(
      "terra could not open the HDF-EOS2 subdatasets in ", hdf_path,
      ". Underlying error: ", as.character(datasets)
    )
  }

  dataset_names <- names(datasets)
  normalized_names <- normalize_layer_name(dataset_names)
  normalized_target <- normalize_layer_name(layer_name)

  matches <- which(
    normalized_names == normalized_target |
      endsWith(normalized_names, paste0(":", normalized_target))
  )

  if (length(matches) != 1) {
    fatal_stop(
      "Expected exactly one '", layer_name, "' subdataset in ", hdf_path,
      "; found ", length(matches), ". Available subdatasets: ",
      paste(dataset_names, collapse = ", ")
    )
  }

  selected <- datasets[matches[[1]]]
  if (!inherits(selected, "SpatRaster")) {
    selected <- rast(selected)
  }
  selected[[1]]
}

valid_par_tile <- function(raster) {
  ifel(raster >= 0 & raster <= 700, raster, NA)
}

valid_quality_tile <- function(raster) {
  valid <- raster == 0 | raster == 1 | raster == 2 | raster == 4
  ifel(valid, raster, NA)
}

project_and_merge <- function(
    hdf_paths,
    layer_name,
    target_template,
    bounds_mask,
    method) {
  source_tiles <- lapply(
    hdf_paths,
    read_hdf_layer,
    layer_name = layer_name
  )

  if (layer_name == par_layer) {
    source_tiles <- lapply(source_tiles, valid_par_tile)
  } else if (layer_name == quality_layer) {
    source_tiles <- lapply(source_tiles, valid_quality_tile)
  }

  projected_tiles <- lapply(
    source_tiles,
    function(tile) {
      project(
        tile,
        target_template,
        method = method,
        mask = TRUE,
        threads = TRUE
      )
    }
  )

  merged <- merge(
    sprc(projected_tiles),
    first = TRUE,
    na.rm = TRUE,
    algo = 1
  )

  mask(merged, bounds_mask)
}

write_par_raster <- function(raster, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeRaster(
    raster,
    path,
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c(
      "COMPRESS=DEFLATE",
      "PREDICTOR=3",
      "TILED=YES",
      "BIGTIFF=IF_SAFER"
    )
  )
}

write_quality_raster <- function(raster, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeRaster(
    raster,
    path,
    overwrite = TRUE,
    datatype = "INT1U",
    NAflag = 255,
    gdal = c(
      "COMPRESS=DEFLATE",
      "PREDICTOR=2",
      "TILED=YES",
      "BIGTIFF=IF_SAFER"
    )
  )
}

count_valid_pixels <- function(raster) {
  count <- global(
    !is.na(raster),
    fun = "sum",
    na.rm = TRUE
  )
  as.integer(count[[1, 1]])
}

safe_remove_temporary_directory <- function(path, temporary_root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  normalized_root <- normalizePath(
    temporary_root,
    winslash = "/",
    mustWork = TRUE
  )

  root_prefix <- paste0(sub("/+$", "", normalized_root), "/")
  if (!startsWith(paste0(normalized_path, "/"), root_prefix)) {
    stop("Refusing to remove path outside temporary HDF root: ", path)
  }

  unlink(path, recursive = TRUE, force = TRUE)
}

process_date <- function(
    target_date,
    urls,
    token,
    target_grid) {
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

  par <- project_and_merge(
    hdf_paths = hdf_paths,
    layer_name = par_layer,
    target_template = target_grid$template,
    bounds_mask = target_grid$mask,
    method = "bilinear"
  )

  quality <- project_and_merge(
    hdf_paths = hdf_paths,
    layer_name = quality_layer,
    target_template = target_grid$template,
    bounds_mask = target_grid$mask,
    method = "near"
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
    valid_pixels = valid_pixels
  )
}


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

validate_bounds(germany_bounds)
dates <- load_unique_dates(date_csv)
token <- earthdata_token()

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checklist <- load_checklist(checklist_path)
bounds_polygon <- make_bounds_polygon(germany_bounds)
target_grid <- make_target_grid(bounds_polygon)

message("Unique SIF dates: ", length(dates))
message(
  "First/last date: ", format(dates[[1]], "%Y-%m-%d"), " / ",
  format(dates[[length(dates)]], "%Y-%m-%d")
)
message("Earthdata product/layer: ", product_name, " / ", par_layer)
message("UTC time: 12:00")
message("Germany WGS84 bounds: ", paste(germany_bounds, collapse = ", "))
message(
  "Output grid: ", ncol(target_grid$template), " x ",
  nrow(target_grid$template), " pixels, ", output_resolution_m, " m, ",
  output_crs
)
message("Output directory: ", output_dir)

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
    previous <- checklist[checklist$date == date_string, , drop = FALSE]
    previous_granules <- if (nrow(previous) > 0) {
      previous$granules_found[[1]]
    } else {
      ""
    }
    previous_pixels <- if (nrow(previous) > 0) {
      previous$valid_par_pixels[[1]]
    } else {
      ""
    }

    checklist <- upsert_checklist(
      checklist,
      make_checklist_row(
        target_date,
        "completed",
        granules_found = previous_granules,
        par_output = paths$par,
        quality_output = paths$quality,
        valid_par_pixels = previous_pixels
      )
    )
    save_checklist(checklist, checklist_path)
    message("  Existing outputs are complete; skipped.")
    next
  }

  result <- tryCatch(
    {
      urls <- search_granule_urls(target_date)
      granule_count <- length(urls)
      message("  Matching Germany tiles: ", granule_count)

      if (granule_count == 0) {
        checklist <<- upsert_checklist(
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
        save_checklist(checklist, checklist_path)
        NULL
      } else {
        processed <- process_date(
          target_date = target_date,
          urls = urls,
          token = token,
          target_grid = target_grid
        )

        checklist <<- upsert_checklist(
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
        save_checklist(checklist, checklist_path)
        message("  Wrote: ", processed$par)
        message("  Wrote: ", processed$quality)
        processed
      }
    },
    fatal_processing_error = function(error) {
      checklist <<- upsert_checklist(
        checklist,
        make_checklist_row(
          target_date,
          "fatal_processing_error",
          par_output = paths$par,
          quality_output = paths$quality,
          error = conditionMessage(error)
        )
      )
      save_checklist(checklist, checklist_path)
      stop(error)
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
      NULL
    }
  )

  invisible(result)
}

completed <- sum(checklist$status == "completed", na.rm = TRUE)
message(
  "Finished. Completed dates in checklist: ",
  completed, "/", length(dates)
)
message("Checklist: ", checklist_path)
