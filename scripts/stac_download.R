library(tidyverse)
library(sf)
library(terra)
library(httr2)
library(jsonlite)

sif_sf <- readRDS("data/ns_sif_mgrs_crop_composition.rds")
sif_sf <- sif_sf %>% select(-geometry)

write.csv(sif_sf, 'data/ns_sif_mgrs_crop_composition.csv',  row.names = FALSE)

sif_sf <- sif_sf[sif_sf$mgrs_tile == '32UNC',]

unique(format(sif_sf$Delta_Date, '%Y-%m'))


stac_url <- "https://geoservice.dlr.de/eoc/ogc/stac/v1"
collection <- "S2_L3A_WASP"
band_asset <- "FRC_B8"

target_tiles <- sort(unique(sif_sf$mgrs_tile))
target_months <- sort(unique(format(as.Date(sif_sf$Delta_Date), "%Y-%m")))

# Use the SIF footprint as the search area. STAC returns all monthly WASP
# products intersecting this bbox; we then filter exactly to the MGRS tiles.
sif_bbox <- sf::st_bbox(sif_sf)
query_bbox <- c(
  unname(sif_bbox["xmin"]) - 0.05,
  unname(sif_bbox["ymin"]) - 0.05,
  unname(sif_bbox["xmax"]) + 0.05,
  unname(sif_bbox["ymax"]) + 0.05
)

month_window <- function(year_month) {
  start <- as.Date(paste0(year_month, "-01"))
  end <- seq(start, by = "1 month", length.out = 2)[2]

  paste0(format(start, "%Y-%m-%d"), "/", format(end, "%Y-%m-%d"))
}

item_tile <- function(item_id) {
  stringr::str_match(item_id, "_T([0-9]{2}[A-Z]{3})_")[, 2]
}

target_item_date <- function(year_month) {
  paste0(stringr::str_replace(year_month, "-", ""), "15")
}

asset_href <- function(feature, asset_name) {
  href <- feature$assets[[asset_name]]$href

  if (is.null(href)) {
    NA_character_
  } else {
    href
  }
}

url_exists <- function(url) {
  if (is.na(url)) {
    return(NA)
  }

  response <- try(
    request(url) |>
      req_method("HEAD") |>
      req_error(is_error = function(resp) FALSE) |>
      req_perform(),
    silent = TRUE
  )

  !inherits(response, "try-error") && resp_status(response) < 400
}

stac_items_request <- function(year_month) {
  request(paste0(stac_url, "/collections/", collection, "/items")) |>
    req_url_query(
      limit = 100,
      bbox = paste(query_bbox, collapse = ","),
      datetime = month_window(year_month),
      f = "json"
    ) |>
    req_error(is_error = function(resp) FALSE)
}

query_month <- function(year_month) {
  message("Querying ", year_month)

  response <- stac_items_request(year_month) |>
    req_perform()

  content_type <- resp_header(response, "content-type")

  if (resp_status(response) >= 400 || !stringr::str_detect(content_type, "json|geo\\+json")) {
    stop(
      "DLR STAC returned status ", resp_status(response),
      " with content-type ", content_type,
      " for ", year_month,
      call. = FALSE
    )
  }

  items <- resp_body_string(response) |>
    jsonlite::fromJSON(simplifyVector = FALSE)

  features <- items$features

  if (length(features) == 0) {
    return(tibble(
      year_month = year_month,
      mgrs_tile = target_tiles,
      item_id = NA_character_,
      band = band_asset,
      asset_href = NA_character_,
      stac_has_item = FALSE,
      stac_has_band = FALSE,
      asset_reachable = NA
    ))
  }

  found <- map_dfr(features, function(feature) {
    id <- feature$id
    tile <- item_tile(id)
    href <- asset_href(feature, band_asset)

    tibble(
      year_month = year_month,
      mgrs_tile = tile,
      item_id = id,
      band = band_asset,
      asset_href = href,
      stac_has_item = TRUE,
      stac_has_band = !is.na(href)
    )
  }) |>
    filter(
      mgrs_tile %in% target_tiles,
      stringr::str_detect(item_id, target_item_date(year_month))
    )

  missing_tiles <- setdiff(target_tiles, found$mgrs_tile)

  missing_rows <- if (length(missing_tiles) > 0) {
    tibble(
      year_month = year_month,
      mgrs_tile = missing_tiles,
      item_id = NA_character_,
      band = band_asset,
      asset_href = NA_character_,
      stac_has_item = FALSE,
      stac_has_band = FALSE
    )
  } else {
    tibble(
      year_month = character(),
      mgrs_tile = character(),
      item_id = character(),
      band = character(),
      asset_href = character(),
      stac_has_item = logical(),
      stac_has_band = logical()
    )
  }

  bind_rows(found, missing_rows) |>
    mutate(asset_reachable = map_lgl(asset_href, url_exists))
}

band8_availability <- map_dfr(target_months, query_month) |>
  arrange(year_month, mgrs_tile, item_id)

availability_summary <- band8_availability |>
  group_by(year_month, mgrs_tile) |>
  summarise(
    stac_items = sum(stac_has_item),
    band8_assets = sum(stac_has_band),
    reachable_band8_assets = sum(asset_reachable, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(year_month, mgrs_tile)

print(band8_availability, n = Inf)
print(availability_summary, n = Inf)

readr::write_csv(
  band8_availability,
  "data/s2_l3a_wasp_band8_stac_availability.csv"
)

readr::write_csv(
  availability_summary,
  "data/s2_l3a_wasp_band8_stac_availability_summary.csv"
)





