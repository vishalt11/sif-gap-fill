library(tidyverse)
library(sf)
library(terra)
library(purrr)

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

#target_col <- "Daily_SIF_740nm"
centroid_crs <- 3035

make_sif_polygon <- function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
  st_polygon(list(rbind(
    c(lon1, lat1),
    c(lon2, lat2),
    c(lon3, lat3),
    c(lon4, lat4),
    c(lon1, lat1)
  )))
}

# Main sif data
final_df <- readRDS('data/extracted_modis_data/modis_1_12.rds')

#757 correction
final_df <- final_df %>%
  mutate(
    sigma_757_daily = Science.SIF_Uncertainty_757nm * Science.daily_correction_factor,
  
    neg_status_757 = case_when(
      Daily_SIF_757nm + 2 * sigma_757_daily >= 0 ~ "accept",
      Daily_SIF_757nm + 2 * sigma_757_daily < 0 &
        Daily_SIF_757nm + 3 * sigma_757_daily >= 0 ~ "questionable",
      Daily_SIF_757nm + 3 * sigma_757_daily < 0 ~ "reject"
    )
  )

#771 correction
final_df <- final_df %>%
  mutate(
    sigma_771_daily = Science.SIF_Uncertainty_771nm * Science.daily_correction_factor,
    
    neg_status_771 = case_when(
      Daily_SIF_771nm + 2 * sigma_771_daily >= 0 ~ "accept",
      Daily_SIF_771nm + 2 * sigma_771_daily < 0 &
        Daily_SIF_771nm + 3 * sigma_771_daily >= 0 ~ "questionable",
      Daily_SIF_771nm + 3 * sigma_771_daily < 0 ~ "reject"
    )
  )

# modis sif (Daily_SIF_757nm + 1.5 * Daily_SIF_771nm) / 2 correction
final_df <- final_df %>%
  mutate(
    sigma_757_daily = Science.SIF_Uncertainty_757nm * Science.daily_correction_factor,
    sigma_771_daily = Science.SIF_Uncertainty_771nm * Science.daily_correction_factor,
    
    target_modis_sif = (Daily_SIF_757nm + 1.5 * Daily_SIF_771nm) / 2,
    
    sigma_modis_sif = 0.5 * sqrt(sigma_757_daily^2 + (1.5 * sigma_771_daily)^2),
    
    neg_status_modis_sif = case_when(
      target_modis_sif + 2 * sigma_modis_sif >= 0 ~ "accept",
      target_modis_sif + 2 * sigma_modis_sif < 0 &
        target_modis_sif + 3 * sigma_modis_sif >= 0 ~ "questionable",
      target_modis_sif + 3 * sigma_modis_sif < 0 ~ "reject"
    )
  )

temp_df <- final_df[final_df$neg_status_modis_sif == 'accept' & final_df$Quality_Flag == 0,]
sif_breaks <- seq(floor(min(temp_df$target_modis_sif, na.rm = TRUE) / 0.25) * 0.25,
                  ceiling(max(temp_df$target_modis_sif, na.rm = TRUE) / 0.25) * 0.25,
                  by = 0.25)
df_binned <- temp_df %>%
  mutate(sif_bin = cut(target_modis_sif, breaks = sif_breaks, include.lowest = TRUE, right = FALSE)) %>%
  select(sif_bin)
table(df_binned$sif_bin)


final_df <- final_df %>%
  mutate(
    bin_check_757 = if_else(Daily_SIF_757nm >= -0.5 & Daily_SIF_757nm < 1.5, "accept", "reject"),
    bin_check_771 = if_else(Daily_SIF_771nm >= -0.5 & Daily_SIF_771nm < 2, "accept", "reject"),
    bin_check_modis_sif = if_else(target_modis_sif >= -0.5 & target_modis_sif < 2, "accept", "reject"),
    final_check_757 = if_else(bin_check_757 == "accept" & neg_status_757 == "accept", "accept", "reject"),
    final_check_771 = if_else(bin_check_771 == "accept" & neg_status_771 == "accept", "accept", "reject"),
    final_check_modis_sif = if_else(bin_check_modis_sif == "accept" & neg_status_modis_sif == "accept", "accept", "reject")
  )

final_df <- final_df %>% select(Daily_SIF_757nm, Daily_SIF_771nm, target_modis_sif,
                                final_check_757, final_check_771, final_check_modis_sif,
                                SZA, VZA, VAz, SAz,
                                Latitude, Longitude, Delta_Time, Delta_Date, sif_doy,
                                Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
                                Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4,
                                Metadata.MeasurementMode, Quality_Flag, 
                                sif_area_km2_evi, state, hzs)

final_df <- final_df[lubridate::month(final_df$Delta_Date) %in% 2:7,]

mgrs_tif_paths <- list.files("data/temp_data/mgrs_tifs", pattern = "\\.tif$", full.names = TRUE)

mgrs_tile_bboxes <- mgrs_tif_paths %>%
  map(function(tif_path) {
    tif_rast <- terra::rast(tif_path)
    tif_ext <- unname(as.vector(terra::ext(tif_rast)))
    tif_bbox <- st_bbox(c(xmin = tif_ext[1], ymin = tif_ext[3], xmax = tif_ext[2], ymax = tif_ext[4]), crs = st_crs(terra::crs(tif_rast)))
    st_sf(mgrs_tile = stringr::str_extract(basename(tif_path), "T[0-9]{2}[A-Z]{3}"), tif_file = basename(tif_path), geometry = st_as_sfc(tif_bbox)) %>%
      st_transform(4326)
  }) %>%
  bind_rows() %>%
  st_make_valid()

mgrs_tile_labels <- mgrs_tile_bboxes %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326)

sif_mgrs_polygons <- final_df %>%
  mutate(source_row = row_number(), target_modis_sif = as.numeric(target_modis_sif), across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric)) %>%
  filter(!is.na(target_modis_sif), if_all(all_of(corner_cols), ~ !is.na(.x))) %>%
  mutate(geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

sif_mgrs_centroids <- sif_mgrs_polygons %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326)

sif_mgrs_centroids_in_tiles <- sif_mgrs_centroids %>%
  st_join(mgrs_tile_bboxes %>% dplyr::select(mgrs_tile), join = st_within, left = FALSE) %>%
  arrange(source_row, mgrs_tile) %>%
  distinct(source_row, .keep_all = TRUE)

# Verification check: should be no rows
# sif_mgrs_centroids_in_tiles %>%
#   st_drop_geometry() %>%
#   count(source_row, name = "n_tiles") %>%
#   filter(n_tiles > 1)
saveRDS(sif_mgrs_centroids_in_tiles, 'data/main_sif_data/9tiles_2_7_M01_QF01.rds')
#-------------------------------------------------------------------------------
# Add sentinel 2 data range compatibility with sif date
final_df <- readRDS('data/main_sif_data/9tiles_2_7_M01_QF01.rds')
final_df <- final_df %>% st_drop_geometry()
dts_temporal_ranges <- readRDS("data/geodes_wasp_zips/dts_r1_temporal_ranges_dates.rds")
dts_temporal_ranges <- dts_temporal_ranges %>% mutate(mgrs_tile = paste0("T", mgrs_tile))

dts_temporal_ranges[117,]$min_days <- 37
dts_temporal_ranges[117,]$temporal_low  <- as.Date('2020-02-07')

sif_date_aligned <- final_df %>%
  mutate(sif_date = as.Date(Delta_Date), sif_year = lubridate::year(sif_date), sif_month = lubridate::month(sif_date)) %>%
  left_join(dts_temporal_ranges %>% dplyr::select(product_path, mgrs_tile, product_year, product_month, product_date, temporal_low, temporal_high, min_days, max_days), by = c("mgrs_tile" = "mgrs_tile", "sif_year" = "product_year", "sif_month" = "product_month")) %>%
  mutate(date_align = case_when(is.na(product_date) ~ "no_product", sif_date >= temporal_low & sif_date <= temporal_high ~ "inrange", TRUE ~ "outrange"))

table(sif_date_aligned$date_align)

sif_date_aligned <- sif_date_aligned %>% filter(date_align %in% c('inrange', 'outrange'))

table(sif_date_aligned$Metadata.MeasurementMode)
sort(colSums(is.na(sif_date_aligned)))

df <- sif_date_aligned %>% select(c(target_modis_sif, Delta_Date))

write_csv(df, 'data/main_sif_data/sif_dates.csv')

saveRDS(sif_date_aligned, 'data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange.rds')

#-------------------------------------------------------------------------------
# Add Phase angle

final_df <- readRDS('data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange.rds')

final_df[71958,]$hzs <- '8b'
final_df[19899,]$hzs <- '8b'

deg_to_rad <- pi / 180
rad_to_deg <- 180 / pi

final_df <- final_df %>%
  mutate(
    # Relative azimuth wrapped to [-180, 180] degrees
    relative_azimuth = ((SAz - VAz + 180) %% 360) - 180,
    
    # Cosine of the angle between surface-to-Sun and
    # surface-to-sensor vectors
    phase_cos = (
      cos(SZA * deg_to_rad) * cos(VZA * deg_to_rad) +
        sin(SZA * deg_to_rad) * sin(VZA * deg_to_rad) *
        cos(relative_azimuth * deg_to_rad)
    ),
    
    # Clamp for floating-point safety, then convert to degrees
    phase_angle = acos(pmax(-1, pmin(1, phase_cos))) * rad_to_deg
  )

final_df <- final_df %>%
  mutate(
    signed_phase_angle = if_else(
      relative_azimuth < 0,
      -phase_angle,
      phase_angle
    )
  )

# updated for new sif upperbounds, 2-7 months M01 QF01
#2022-07-29
#2022-07-31
#2024-05-31
#2024-06-02
#2024-07-13
#2024-07-15
#2024-07-29
#remove sif rows that dont have PAR data
dates_to_remove <- as.Date(c(
  "2022-07-29",
  "2022-07-31",
  "2024-05-31",
  "2024-06-02",
  "2024-07-13",
  "2024-07-15",
  "2024-07-29"
))

final_df <- final_df |>
  dplyr::filter(!Delta_Date %in% dates_to_remove)


final_df <- final_df %>% select(Daily_SIF_757nm, Daily_SIF_771nm, target_modis_sif,
                                final_check_757, final_check_771, final_check_modis_sif,
                                Latitude, Longitude, Delta_Time, sif_doy, 
                                Delta_Date, temporal_low, temporal_high, date_align,
                                Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
                                Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4,
                                sif_area_km2_evi, Metadata.MeasurementMode, Quality_Flag, 
                                state, hzs, mgrs_tile, product_path, 
                                SZA, VZA, VAz, SAz, phase_angle, signed_phase_angle)


colSums(is.na(final_df))

saveRDS(final_df, 'data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.rds')
write_csv(final_df, 'data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.csv')

summary(final_df[final_df$date_align == 'inrange' & final_df$final_check_modis_sif == 'accept',]$target_modis_sif)
summary(final_df[final_df$final_check_modis_sif == 'accept',]$target_modis_sif)


inrange_df <- final_df %>% filter(date_align == 'inrange')
saveRDS(inrange_df, 'data/main_sif_data/9tiles_2_7_M01_QF01_inrange_PARrm.rds')
write_csv(inrange_df, 'data/main_sif_data/9tiles_2_7_M01_QF01_inrange_PARrm.csv')

temp_df <- final_df %>% select(c(target_modis_sif, Delta_Date, temporal_low, temporal_high, date_align))
head(temp_df)

temp_df <- temp_df %>%
  mutate(
    temporal_range_days = as.numeric(temporal_high - temporal_low),
    temporal_mid_date = temporal_low + floor(temporal_range_days / 2),
    abs_delta_from_temporal_mid_days = abs(as.numeric(Delta_Date - temporal_mid_date))
  )

summary(temp_df[temp_df$date_align == 'outrange',]$temporal_range_days)
summary(temp_df[temp_df$date_align == 'outrange',]$abs_delta_from_temporal_mid_days)

summary(temp_df[temp_df$date_align == 'inrange',]$temporal_range_days)
summary(temp_df[temp_df$date_align == 'inrange',]$abs_delta_from_temporal_mid_days)

temp_df[temp_df$temporal_range_days == 0,]



dts_temporal_ranges[dts_temporal_ranges$temporal_low == as.Date('2020-02-07'),]

