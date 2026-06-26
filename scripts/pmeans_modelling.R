library(tidyverse)
library(sf)
library(lubridate)

set.seed(42)

model_data_path <- "data/model_data_wpar.csv"
target_col <- "Daily_SIF_740nm"
train_fraction <- 0.8
area_crs <- 3035
pixel_area_m2 <- 10 * 10

df <- read_csv(
  model_data_path,
  col_types = cols(
    Delta_Time = col_character(),
    Delta_Date = col_date(),
    sif_id = col_character(),
    mgrs_tile = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

required_cols <- c(
  target_col,
  "Delta_Time", "Delta_Date",
  "Latitude", "Longitude",
  corner_cols,
  "crop_pixel_count",
  "ww_pct",
  "mean_fapar", "mean_par", "apar"
)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

parse_delta_hour <- function(delta_time) {
  if (inherits(delta_time, "POSIXt")) {
    return(hour(delta_time))
  }

  if (is.numeric(delta_time)) {
    return(if_else(
      !is.na(delta_time) & delta_time >= 0 & delta_time < 24,
      as.integer(floor(delta_time)),
      hour(as.POSIXct(delta_time, origin = "1990-01-01", tz = "UTC")),
      missing = NA_integer_
    ))
  }

  delta_time_chr <- as.character(delta_time)
  parsed_time <- suppressWarnings(parse_date_time(
    delta_time_chr,
    orders = c(
      "ymd HMS z", "ymd HM z",
      "ymd HMS", "ymd HM",
      "ymd IMS p", "ymd IM p",
      "HMS", "HM", "IMS p", "IM p"
    ),
    tz = "UTC",
    quiet = TRUE
  ))
  parsed_hour <- hour(parsed_time)

  numeric_time <- suppressWarnings(as.numeric(delta_time_chr))
  numeric_hour <- if_else(
    !is.na(numeric_time) & numeric_time >= 0 & numeric_time < 24,
    as.integer(floor(numeric_time)),
    hour(as.POSIXct(numeric_time, origin = "1990-01-01", tz = "UTC")),
    missing = NA_integer_
  )

  coalesce(parsed_hour, numeric_hour)
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

active_growth_months <- list(
  crop_count_winter_wheat = c(2:7, 10:11),
  crop_count_winter_barley = c(2:6, 10:11),
  crop_count_winter_rye = c(2:7, 10:11),
  crop_count_other_winter_cereals = c(2:7, 10:11),
  crop_count_spring_wheat = 3:8,
  crop_count_spring_barley = 3:8,
  crop_count_spring_oat = 3:8,
  crop_count_maize = 5:10,
  crop_count_legumes = 4:9,
  crop_count_potato = 4:9,
  crop_count_sugar_beet = 4:10,
  crop_count_rapeseed = c(2:7, 9:11),
  crop_count_clover_alfalfa = 3:10,
  crop_count_arable_grass = 3:11,
  crop_count_permanent_grassland = 3:11,
  crop_count_vineyard = 4:10,
  crop_count_fruit_trees_and_other_woody_vegetation = 3:10,
  crop_count_hops = 4:9,
  crop_count_other_agricultural_use = 3:10
)

calculate_active_growth_crop_pixels <- function(data, active_growth_months) {
  active_growth_cols <- intersect(names(active_growth_months), names(data))

  if (length(active_growth_cols) == 0) {
    return(rep(NA_real_, nrow(data)))
  }

  reduce(
    active_growth_cols,
    .init = rep(0, nrow(data)),
    .f = function(total, crop_col) {
      crop_count <- replace_na(as.numeric(data[[crop_col]]), 0)
      is_active_month <- data$month %in% active_growth_months[[crop_col]]

      total + if_else(is_active_month, crop_count, 0)
    }
  )
}

sif_sf <- df %>%
  mutate(
    across(all_of(corner_cols), as.numeric),
    crop_pixel_count = as.numeric(crop_pixel_count)
  ) %>%
  filter(if_all(all_of(corner_cols), ~ !is.na(.x))) %>%
  mutate(
    geometry = pmap(
      list(
        Lon_corner1, Lat_corner1,
        Lon_corner2, Lat_corner2,
        Lon_corner3, Lat_corner3,
        Lon_corner4, Lat_corner4
      ),
      make_sif_polygon
    )
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

sif_sf$sif_area_m2 <- as.numeric(st_area(st_transform(sif_sf, area_crs)))

sif_sf <- sif_sf %>%
  mutate(
    pct_vegetation = if_else(
      is.na(sif_area_m2) | sif_area_m2 <= 0,
      NA_real_,
      100 * crop_pixel_count * pixel_area_m2 / sif_area_m2
    )
  )

drop_cols <- c(
  corner_cols,
  "Quality_Flag",
  "Meteo.specific_humidity",
  "Meteo.surface_pressure",
  "Meteo.temperature_skin",
  "Meteo.temperature_two_meter",
  "Meteo.vapor_pressure_deficit",
  "Metadata.MeasurementMode",
  "file_id",
  "source_file",
  "sif_id",
  "sif_year",
  "mgrs_tile",
  "sif_extract_id",
  "wasp_product_id",
  "wasp_year_month",
  "non_crop_low_ndvi_threshold",
  "contaminated_20m_fraction_threshold",
  "sif_row_id",
  "sif_doy",
  "closest_doy",
  "closest_glass_date",
  "glass_day_diff"
)

crop_count_cols <- names(sif_sf) %>%
  str_subset("^crop_count_")

missing_active_growth_cols <- setdiff(names(active_growth_months), crop_count_cols)
if (length(missing_active_growth_cols) > 0) {
  warning("Missing crop-count columns for active-growth feature: ", paste(missing_active_growth_cols, collapse = ", "))
}

sif_sf[is.na(sif_sf$ww_pct) == TRUE,]$ww_pct <- 0

model_df <- sif_sf %>%
  st_drop_geometry() %>%
  mutate(
    sif = .data[[target_col]],
    delta_hour = parse_delta_hour(Delta_Time),
    month = month(Delta_Date),
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude),
    ww_pct = as.numeric(ww_pct),
    across(starts_with("mean_"), as.numeric),
    across(all_of(crop_count_cols), as.numeric),
    apar = as.numeric(apar)
  )

model_df$active_growth_crop_pixel_count <- calculate_active_growth_crop_pixels(
  model_df,
  active_growth_months
)

model_df <- model_df %>%
  mutate(
    active_growth_crop_pct = if_else(
      is.na(crop_pixel_count) | crop_pixel_count <= 0,
      0,
      active_growth_crop_pixel_count / crop_pixel_count
    )
  ) %>%
  select(-any_of(c(drop_cols, crop_count_cols))) %>%
  select(
    sif,
    delta_hour,
    month,
    Latitude,
    Longitude,
    starts_with("mean_"),
    apar,
    ww_pct,
    active_growth_crop_pct,
    pct_vegetation
  ) %>%
  drop_na()

if (nrow(model_df) < 2) {
  stop("Need at least two complete modeling rows after preprocessing.")
}

model_df <- model_df %>% select(-c(delta_hour))

write.csv(model_df, 'data/final_model.csv', row.names = FALSE)





