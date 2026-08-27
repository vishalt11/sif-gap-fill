# Diagnose why the VNP18A2 CMR search returns more granules than expected.
#
# This script is metadata-only. It does not download HDF5 files, create rasters,
# or modify the existing VNP18A2 outputs. Run it and share the complete console
# output when finished.

suppressPackageStartupMessages({
  library(httr2)
})


# -----------------------------------------------------------------------------
# Test configuration
# -----------------------------------------------------------------------------

target_date <- as.Date("2022-03-27")
collection_concept_id <- "C2631841566-LPCLOUD"
germany_bounds <- c(5.5, 47.0, 15.5, 55.2) # west, south, east, north
cmr_granules_url <- "https://cmr.earthdata.nasa.gov/search/granules.json"

request_timeout_seconds <- 60
request_retries <- 4


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0) right else left
}

scalar_text <- function(value, default = NA_character_) {
  if (is.null(value) || length(value) == 0) {
    return(default)
  }
  as.character(value[[1]])
}

acquisition_token <- function(date) {
  paste0("A", format(date, "%Y"), format(date, "%j"))
}

cmr_entry_data_url <- function(entry) {
  links <- entry$links %||% list()
  if (length(links) == 0) {
    return(NA_character_)
  }

  candidates <- vapply(
    links,
    function(link) {
      href <- scalar_text(link$href, "")
      rel <- scalar_text(link$rel, "")
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

  cloud_link <- grepl(
    "data.lpdaac.earthdatacloud.nasa.gov",
    candidates,
    fixed = TRUE
  )
  if (any(cloud_link)) {
    return(candidates[which(cloud_link)[[1]]])
  }

  candidates[[1]]
}

parse_vnp18a2_filename <- function(url) {
  filename <- if (is.na(url) || !nzchar(url)) {
    NA_character_
  } else {
    basename(sub("\\?.*$", "", url))
  }

  pattern <- paste0(
    "^(VNP18A2)\\.A([0-9]{4})([0-9]{3})\\.",
    "(h[0-9]{2}v[0-9]{2})\\.([0-9]{3})\\.",
    "([0-9]{13})\\.h5$"
  )
  match <- regmatches(filename, regexec(pattern, filename))[[1]]

  if (length(match) == 0) {
    return(list(
      filename = filename,
      acquisition = NA_character_,
      tile = NA_character_,
      collection = NA_character_,
      production_time = NA_character_
    ))
  }

  list(
    filename = filename,
    acquisition = paste0("A", match[[3]], match[[4]]),
    tile = match[[5]],
    collection = match[[6]],
    production_time = match[[7]]
  )
}

search_cmr <- function(temporal, query_label) {
  response <- request(cmr_granules_url) |>
    req_url_query(
      concept_id = collection_concept_id,
      temporal = temporal,
      bounding_box = paste(germany_bounds, collapse = ","),
      downloadable = "true",
      page_size = 200
    ) |>
    req_headers(Accept = "application/json") |>
    req_user_agent("VNP18A2-granule-date-diagnostic/1.0") |>
    req_timeout(request_timeout_seconds) |>
    req_retry(max_tries = request_retries) |>
    req_perform()

  body <- resp_body_json(response, simplifyVector = FALSE)
  entries <- body$feed$entry %||% list()
  cmr_hits <- resp_header(response, "CMR-Hits") %||% as.character(length(entries))

  if (length(entries) == 0) {
    return(list(
      label = query_label,
      temporal = temporal,
      cmr_hits = cmr_hits,
      entries_received = 0L,
      data = data.frame()
    ))
  }

  rows <- lapply(seq_along(entries), function(index) {
    entry <- entries[[index]]
    url <- cmr_entry_data_url(entry)
    parsed <- parse_vnp18a2_filename(url)

    data.frame(
      entry_number = index,
      cmr_id = scalar_text(entry$id),
      title = scalar_text(entry$title),
      producer_granule_id = scalar_text(entry$producer_granule_id),
      time_start = scalar_text(entry$time_start),
      time_end = scalar_text(entry$time_end),
      updated = scalar_text(entry$updated),
      filename = parsed$filename,
      acquisition = parsed$acquisition,
      tile = parsed$tile,
      collection = parsed$collection,
      production_time = parsed$production_time,
      url = url,
      stringsAsFactors = FALSE
    )
  })

  data <- do.call(rbind, rows)
  data <- data[!duplicated(data$url), , drop = FALSE]
  rownames(data) <- NULL

  list(
    label = query_label,
    temporal = temporal,
    cmr_hits = cmr_hits,
    entries_received = length(entries),
    data = data
  )
}

classify_acquisition <- function(acquisition, target, previous, next_day) {
  ifelse(
    is.na(acquisition),
    "UNPARSED",
    ifelse(
      acquisition == target,
      "TARGET_DATE",
      ifelse(
        acquisition == previous,
        "PREVIOUS_DATE",
        ifelse(acquisition == next_day, "NEXT_DATE", "OTHER_DATE")
      )
    )
  )
}

print_query_report <- function(result, expected, previous, next_day) {
  cat("\n", paste(rep("=", 88), collapse = ""), "\n", sep = "")
  cat("QUERY: ", result$label, "\n", sep = "")
  cat("Temporal parameter: ", result$temporal, "\n", sep = "")
  cat("CMR-Hits header: ", result$cmr_hits, "\n", sep = "")
  cat("Entries received: ", result$entries_received, "\n", sep = "")
  cat("Unique selected HDF5 URLs: ", nrow(result$data), "\n", sep = "")

  if (nrow(result$data) == 0) {
    cat("No downloadable VNP18A2 HDF5 granules were returned.\n")
    return(invisible(NULL))
  }

  data <- result$data
  data$date_class <- classify_acquisition(
    data$acquisition,
    expected,
    previous,
    next_day
  )

  cat("\nCounts by acquisition-date classification:\n")
  print(table(data$date_class, useNA = "ifany"))

  cat("\nCounts by acquisition token and spatial tile:\n")
  print(
    as.data.frame(table(data$acquisition, data$tile, useNA = "ifany")),
    row.names = FALSE
  )

  detail_columns <- c(
    "entry_number",
    "date_class",
    "acquisition",
    "tile",
    "collection",
    "production_time",
    "time_start",
    "time_end",
    "filename"
  )
  cat("\nGranule details in CMR return order:\n")
  print(data[, detail_columns, drop = FALSE], row.names = FALSE, right = FALSE)

  cat("\nFull selected data URLs:\n")
  for (index in seq_len(nrow(data))) {
    cat(sprintf("[%02d] %s\n", index, data$url[[index]]))
  }

  target_data <- data[data$acquisition == expected, , drop = FALSE]
  cat("\nTarget-date diagnostic:\n")
  cat("  Target-date granules: ", nrow(target_data), "\n", sep = "")
  cat(
    "  Distinct target-date tiles: ",
    length(unique(target_data$tile[!is.na(target_data$tile)])),
    "\n",
    sep = ""
  )

  if (nrow(target_data) > 0) {
    duplicate_key <- paste(target_data$acquisition, target_data$tile, sep = "|")
    duplicate_rows <- target_data[
      duplicated(duplicate_key) | duplicated(duplicate_key, fromLast = TRUE),
      ,
      drop = FALSE
    ]

    if (nrow(duplicate_rows) == 0) {
      cat("  No duplicate target-date tile IDs were found.\n")
    } else {
      cat("  DUPLICATE target-date tile IDs were found:\n")
      print(
        duplicate_rows[, c("tile", "production_time", "filename"), drop = FALSE],
        row.names = FALSE,
        right = FALSE
      )
    }
  }

  invisible(data)
}


# -----------------------------------------------------------------------------
# Run two metadata searches
# -----------------------------------------------------------------------------

date_string <- format(target_date, "%Y-%m-%d")
expected <- acquisition_token(target_date)
previous <- acquisition_token(target_date - 1)
next_day <- acquisition_token(target_date + 1)

whole_day_temporal <- paste0(
  date_string, "T00:00:00Z,",
  date_string, "T23:59:59Z"
)

# A narrow midday query should remove granules that only touch the requested
# date at 00:00:00 or otherwise overlap a UTC-day boundary.
midday_temporal <- paste0(
  date_string, "T12:00:00Z,",
  date_string, "T12:00:01Z"
)

cat(paste(rep("#", 88), collapse = ""), "\n", sep = "")
cat("VNP18A2.002 CMR GRANULE DIAGNOSTIC\n")
cat(paste(rep("#", 88), collapse = ""), "\n", sep = "")
cat("Target date: ", date_string, "\n", sep = "")
cat("Expected filename acquisition token: ", expected, "\n", sep = "")
cat("Previous-day token: ", previous, "\n", sep = "")
cat("Next-day token: ", next_day, "\n", sep = "")
cat("Collection concept ID: ", collection_concept_id, "\n", sep = "")
cat("Germany bounds: ", paste(germany_bounds, collapse = ", "), "\n", sep = "")
cat("This script performs CMR metadata requests only; it downloads no rasters.\n")

whole_day <- search_cmr(whole_day_temporal, "Downloader's current whole-UTC-day query")
midday <- search_cmr(midday_temporal, "Narrow midday query")

whole_day_data <- print_query_report(
  whole_day,
  expected,
  previous,
  next_day
)
midday_data <- print_query_report(
  midday,
  expected,
  previous,
  next_day
)

cat("\n", paste(rep("=", 88), collapse = ""), "\n", sep = "")
cat("AUTOMATIC INTERPRETATION\n")
cat(paste(rep("=", 88), collapse = ""), "\n", sep = "")

if (is.null(whole_day_data) || nrow(whole_day_data) == 0) {
  cat("No whole-day granules were available to interpret.\n")
} else {
  whole_target <- whole_day_data$acquisition == expected
  whole_target[is.na(whole_target)] <- FALSE
  non_target_count <- sum(!whole_target)

  target_rows <- whole_day_data[whole_target, , drop = FALSE]
  target_key <- paste(target_rows$acquisition, target_rows$tile, sep = "|")
  duplicate_target_tiles <- anyDuplicated(target_key) > 0

  if (non_target_count > 0) {
    cat(
      "The whole-day search returned ", non_target_count,
      " granule(s) whose filename acquisition day is not ", expected, ".\n",
      sep = ""
    )
    cat(
      "This supports the adjacent/other-date temporal-overlap explanation.\n"
    )
  } else {
    cat("Every whole-day result has the expected acquisition token ", expected, ".\n", sep = "")
  }

  if (duplicate_target_tiles) {
    cat(
      "At least one target-date tile has multiple production copies. ",
      "This supports the reprocessed-duplicate explanation.\n",
      sep = ""
    )
  } else {
    cat("No repeated target-date tile ID was found.\n")
  }

  if (!is.null(midday_data)) {
    cat(
      "Whole-day versus midday unique URL counts: ",
      nrow(whole_day_data), " versus ", nrow(midday_data), ".\n",
      sep = ""
    )
    if (nrow(midday_data) < nrow(whole_day_data)) {
      cat(
        "The reduced midday count indicates that some whole-day hits only ",
        "overlap the temporal query near a day boundary.\n",
        sep = ""
      )
    }
  }
}

cat("\nPlease share the complete console output from this diagnostic.\n")

