library(tidyverse)
library(sf)
library(terra)


#-------------------------------------------------------------------------------
# Get temporal ranges for each tile-year-month combo
# read only 32UQV and 32UMC geodes folders
geodes_root <- "data/geodes_wasp_zips"
target_mgrs_tiles <- c("32UMC", "32UQV")
target_geodes_roots <- file.path(geodes_root, target_mgrs_tiles)
sentinel_product_pattern <- "^SENTINEL2[A-Z]_\\d{8}-000000-000_L3A_T\\d{2}[A-Z]{3}_C_V[0-9]-[0-9]$"
dts_workers <- 2

dts_product_index <- target_geodes_roots %>%
  map(~ list.dirs(.x, recursive = TRUE, full.names = TRUE)) %>%
  unlist(use.names = FALSE) %>%
  keep(~ stringr::str_detect(basename(.x), sentinel_product_pattern)) %>%
  tibble(product_path = .) %>%
  mutate(product_id = basename(product_path),
         mgrs_tile = stringr::str_remove(stringr::str_extract(product_id, "T\\d{2}[A-Z]{3}"), "^T"),
         product_date = as.Date(stringr::str_extract(product_id, "\\d{8}"), format = "%Y%m%d"),
         product_year = lubridate::year(product_date),
         product_month = lubridate::month(product_date),
         product_version = stringr::str_extract(product_id, "V[0-9]-[0-9]$"),
         dts_r1_path = map_chr(product_path, ~ list.files(file.path(.x, "MASKS"), pattern = "DTS_R1\\.tif$", full.names = TRUE)[1])) %>%
  filter(!is.na(dts_r1_path))

summarise_dts_r1 <- function(dts_r1_path) {
  r <- terra::rast(dts_r1_path)
  terra::global(r, fun = function(x, ...) {
    c(
      min = min(x, na.rm = TRUE),
      q10 = as.numeric(stats::quantile(x, 0.1, na.rm = TRUE, names = FALSE)),
      median = median(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }) %>%
    as_tibble(rownames = "raster_layer")
}

dts_old_plan <- future::plan()
future::plan(future::multisession, workers = dts_workers)

progressr::handlers(global = TRUE)

dts_temporal_ranges <- progressr::with_progress({
  p <- progressr::progressor(along = dts_product_index$dts_r1_path)
  dts_product_index %>%
    mutate(dts_summary = furrr::future_map2(dts_r1_path, product_id, function(path, id) {
      result <- summarise_dts_r1(path)
      p(sprintf("done %s", id))
      result
    }, .options = furrr::furrr_options(seed = FALSE))) %>%
    unnest(dts_summary) %>%
    arrange(mgrs_tile, product_year, product_month)
})

future::plan(dts_old_plan)



dts_temporal_ranges[dts_temporal_ranges$min == -10000,]$min <- 0
dts_temporal_ranges[dts_temporal_ranges$q10 == -10000,]$q10 <- 0
dts_temporal_ranges[dts_temporal_ranges$median == -10000,]$median <- 0

#saveRDS(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges.rds")
#dts_temporal_ranges <- readRDS("data/geodes_wasp_zips/dts_r1_temporal_ranges.rds")

dts_temporal_ranges <- dts_temporal_ranges %>%
  mutate(min_days = case_when(min != 0 ~ min, q10 != 0 ~ q10, median != 0 ~ median, TRUE ~ NA_real_),
         max_days = max,
         year_start = as.Date(paste0(product_year, "-01-01")),
         temporal_low = year_start + min_days,
         temporal_high = year_start + max_days) %>%
  dplyr::select(-year_start)

dts_temporal_ranges %>% select(c(mgrs_tile, product_date, min_days, max_days, temporal_low, temporal_high))


colSums(is.na(dts_temporal_ranges))

saveRDS(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges_dates_UMC_UQV.rds")


