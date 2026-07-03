library(tidyverse)
library(sf)
library(purrr)

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

#required_cols <- c("Latitude", "Longitude", corner_cols)
target_col <- "Daily_SIF_740nm"
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



# flag extreme sif rows depending on certain criteria
final_df <- readRDS('data/extracted_modis_data/modis_1_12.rds')

#high_sif <- final_df[final_df$Daily_SIF_740nm > 0,]
#high_sif %>% select(Daily_SIF_757nm, mean_evi, mean_ndvi, mean_fapar)
#plot(high_sif$mean_evi, high_sif$Daily_SIF_740nm)


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

table(final_df$neg_status_757)
summary(final_df[final_df$neg_status_757 == 'accept',]$Daily_SIF_757nm)

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

table(final_df$neg_status_771)
summary(final_df[final_df$neg_status_771 == 'accept',]$Daily_SIF_771nm)

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

table(final_df$neg_status_modis_sif)
summary(final_df[final_df$neg_status_modis_sif == 'accept',]$target_modis_sif)

sif_breaks <- seq(floor(min(final_df$target_modis_sif, na.rm = TRUE) / 0.25) * 0.25,
                  ceiling(max(final_df$target_modis_sif, na.rm = TRUE) / 0.25) * 0.25,
                  by = 0.25)

df_binned <- final_df %>%
  mutate(sif_bin = cut(target_modis_sif, breaks = sif_breaks, include.lowest = TRUE, right = FALSE)) %>%
  select(sif_bin)

table(df_binned$sif_bin)

final_df <- final_df %>%
  mutate(
    bin_check_757 = if_else(Daily_SIF_757nm >= -0.5 & Daily_SIF_757nm < 1.25, "accept", "reject"),
    bin_check_771 = if_else(Daily_SIF_771nm >= -0.5 & Daily_SIF_771nm < 1, "accept", "reject"),
    bin_check_modis_sif = if_else(target_modis_sif >= -0.5 & target_modis_sif < 1.25, "accept", "reject"),
    final_check_757 = if_else(bin_check_757 == "accept" & neg_status_757 == "accept", "accept", "reject"),
    final_check_771 = if_else(bin_check_771 == "accept" & neg_status_771 == "accept", "accept", "reject"),
    final_check_modis_sif = if_else(bin_check_modis_sif == "accept" & neg_status_modis_sif == "accept", "accept", "reject")
  )

summary(final_df[final_df$final_check_757 == 'accept',]$Daily_SIF_757nm)
nrow(final_df[final_df$final_check_757 == 'accept',])

saveRDS(final_df, 'data/extracted_modis_data/modis_1_12_bin_uncertainity_corrected_raw.rds')

final_df <- readRDS('data/extracted_modis_data/modis_1_12_bin_uncertainity_corrected_raw.rds')
colnames(final_df)

final_df <- final_df %>% select(Daily_SIF_757nm, Daily_SIF_771nm, target_modis_sif,
                                Latitude, Longitude, Delta_Time,
                                Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
                                Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4,
                                sif_area_km2_evi, Delta_Date, Metadata.MeasurementMode, Quality_Flag, state,
                                hzs, final_check_757, final_check_771, final_check_modis_sif)

final_df <- final_df[final_df$Quality_Flag == 0 & final_df$Metadata.MeasurementMode == 0,]

saveRDS(final_df, 'data/extracted_modis_data/modis_1_12_bin_uncertainity_corrected_M0QF0_cnn.rds')

final_df <- final_df[lubridate::month(final_df$Delta_Date) %in% 2:7,]

saveRDS(final_df, 'data/extracted_modis_data/modis_2_7_bin_uncertainity_corrected_M0QF0_cnn.rds')
write_csv(final_df, 'data/extracted_modis_data/modis_2_7_bin_uncertainity_corrected_M0QF0_cnn.csv')
#IQR rule

# q <- quantile(final_df$Daily_SIF_757nm, probs = c(0.25, 0.75), na.rm = TRUE)
# iqr <- q[2] - q[1]
# lower <- q[1] - 1.5 * iqr
# upper <- q[2] + 1.5 * iqr
# c(upper, lower)
# upper <- 1.5
# final_df <- final_df %>% filter(Daily_SIF_740nm >= lower, Daily_SIF_740nm <= upper)

sort(colSums(is.na(final_df)))

#-------------------------------------------------------------------------------
# checking quality of extreme sif values based on neighbour values

final_sif_polygons <- final_df %>%
  mutate(sif_row_id = row_number(), sif_date = as.Date(Delta_Date), across(all_of(c("Latitude", "Longitude", corner_cols, target_col)), as.numeric)) %>%
  mutate(geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

final_sif_centroids_3035 <- final_sif_polygons %>%
  st_transform(centroid_crs) %>%
  st_centroid()

extreme_samples <- final_sif_polygons %>%
  st_drop_geometry() %>%
  mutate(
    extreme_bin = case_when(
      Daily_SIF_740nm >= -0.75 & Daily_SIF_740nm < -0.5 ~ "[-0.75,-0.5)",
      Daily_SIF_740nm >= -0.5 & Daily_SIF_740nm < -0.25 ~ "[-0.5,-0.25)",
      #Daily_SIF_740nm >= 1 & Daily_SIF_740nm < 1.25 ~ "[1,1.25)",
      Daily_SIF_740nm >= 1.25 & Daily_SIF_740nm < 1.5 ~ "[1.25,1.5)",
      Daily_SIF_740nm >= 1.5 & Daily_SIF_740nm <= 1.75 ~ "[1.5,1.75]",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(extreme_bin)) %>%
  mutate(sample_id = row_number(), sample_label = paste0("sample ", sample_id, " | ", extreme_bin, " | SIF=", round(Daily_SIF_740nm, 3), " | ", sif_date))

nearest_neighbor_rows <- map_dfr(seq_len(nrow(extreme_samples)), function(i) {
  sample_row <- extreme_samples[i, ]
  sample_centroid <- final_sif_centroids_3035 %>% filter(sif_row_id == sample_row$sif_row_id)
  same_date_centroids <- final_sif_centroids_3035 %>% filter(sif_date == sample_row$sif_date, sif_row_id != sample_row$sif_row_id)
  if (nrow(same_date_centroids) == 0) return(tibble())
  same_date_centroids %>%
    st_drop_geometry() %>%
    mutate(sample_id = sample_row$sample_id, sample_sif_row_id = sample_row$sif_row_id, extreme_bin = sample_row$extreme_bin, neighbor_distance_m = as.numeric(st_distance(sample_centroid, same_date_centroids))) %>%
    filter(neighbor_distance_m <= 3000) %>%
    arrange(neighbor_distance_m) %>%
    slice_head(n = 14) %>%
    transmute(sample_id, sample_sif_row_id, extreme_bin, neighbor_sif_row_id = sif_row_id, neighbor_distance_m)
})

extreme_sample_polygons <- final_sif_polygons %>%
  inner_join(extreme_samples %>% dplyr::select(sif_row_id, sample_id, extreme_bin, sample_label), by = "sif_row_id") %>%
  mutate(leaflet_group = paste0("Extreme samples ", extreme_bin), map_label = sample_label)

neighbor_polygons <- final_sif_polygons %>%
  inner_join(nearest_neighbor_rows, by = c("sif_row_id" = "neighbor_sif_row_id")) %>%
  mutate(leaflet_group = paste0("Neighbors ", extreme_bin), map_label = paste0("neighbor of sample ", sample_id, " | SIF=", round(Daily_SIF_740nm, 3), " | dist=", round(neighbor_distance_m, 1), " m | ", sif_date))

neighbor_summary <- neighbor_polygons %>%
  st_drop_geometry() %>%
  group_by(sample_id, sample_sif_row_id, extreme_bin) %>%
  summarise(
    n_neighbors = n(),
    neighbor_sif_min = min(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_mean = mean(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_median = median(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_max = max(Daily_SIF_740nm, na.rm = TRUE),
    local_sif_mad = median(abs(Daily_SIF_740nm - neighbor_sif_median), na.rm = TRUE),
    max_d = max(neighbor_distance_m, na.rm = TRUE),
    neighbor_active_growth_pct_mean = mean(active_growth_pct, na.rm = TRUE),
    neighbor_crop_pixel_count_mean = mean(crop_pixel_count, na.rm = TRUE),
    neighbor_mean_fapar_mean = mean(mean_fapar, na.rm = TRUE),
    neighbor_mean_evi_mean = mean(mean_evi, na.rm = TRUE),
    neighbor_mean_ndvi_mean = mean(mean_ndvi, na.rm = TRUE),
    .groups = "drop"
  )

min_neighbors_for_reliability <- 5
local_outlier_score_threshold <- 4
local_mad_floor <- 0.01

extreme_sample_neighbor_summary <- extreme_samples %>%
  dplyr::select(sample_id, sample_sif_row_id = sif_row_id, sample_date = sif_date, extreme_bin, sample_sif = Daily_SIF_740nm, sample_active_growth_pct = active_growth_pct, sample_crop_pixel_count = crop_pixel_count, sample_mean_fapar = mean_fapar, sample_mean_evi = mean_evi, sample_mean_ndvi = mean_ndvi) %>%
  left_join(neighbor_summary, by = c("sample_id", "sample_sif_row_id", "extreme_bin")) %>%
  mutate(
    n_neighbors = replace_na(n_neighbors, 0L),
    local_outlier_score = abs(sample_sif - neighbor_sif_median) / pmax(local_sif_mad, local_mad_floor, na.rm = TRUE),
    sif_reliability = case_when(
      n_neighbors < min_neighbors_for_reliability | is.na(local_outlier_score) ~ "keep",
      local_outlier_score >= local_outlier_score_threshold & sample_sif < neighbor_sif_median ~ "suspicious_low",
      local_outlier_score >= local_outlier_score_threshold & sample_sif > neighbor_sif_median ~ "suspicious_high",
      TRUE ~ "keep"
    )
  ) %>%
  relocate(sample_id, sample_sif_row_id, sample_date, extreme_bin, sif_reliability, sample_sif, sample_active_growth_pct, sample_crop_pixel_count, sample_mean_fapar, sample_mean_evi, sample_mean_ndvi, n_neighbors, max_d, neighbor_sif_min, neighbor_sif_mean, neighbor_sif_median, neighbor_sif_max, local_sif_mad, local_outlier_score, neighbor_active_growth_pct_mean, neighbor_crop_pixel_count_mean, neighbor_mean_fapar_mean, neighbor_mean_evi_mean, neighbor_mean_ndvi_mean) %>%
  arrange(extreme_bin, sample_id)

print(extreme_sample_neighbor_summary)




# ex_sum <- extreme_sample_neighbor_summary %>% 
#   select(c(sample_id, sif_reliability, local_outlier_score, local_sif_mad, n_neighbors, max_d, sample_sif, neighbor_sif_mean,  
#            sample_active_growth_pct,  neighbor_active_growth_pct_mean, 
#            sample_crop_pixel_count, neighbor_crop_pixel_count_mean,
#            sample_mean_fapar, neighbor_mean_fapar_mean,
#            sample_mean_evi, neighbor_mean_evi_mean, 
#            sample_mean_ndvi, neighbor_mean_ndvi_mean)) %>% 
#   mutate(diff_sif = round(abs(neighbor_sif_mean - sample_sif),3),
#          diff_active = round(abs(neighbor_active_growth_pct_mean - sample_active_growth_pct),3),
#          diff_crop = round(abs(neighbor_crop_pixel_count_mean - sample_crop_pixel_count),3),
#          diff_fapar = round(abs(neighbor_mean_fapar_mean - sample_mean_fapar),3),
#          diff_evi = round(abs(neighbor_mean_evi_mean - sample_mean_evi),3),
#          diff_ndvi = round(abs(neighbor_mean_ndvi_mean - sample_mean_ndvi),3))
# 
# ex_sum %>% select(sample_id, sif_reliability, local_outlier_score, local_sif_mad, max_d, n_neighbors, sample_sif, neighbor_sif_mean, diff_sif, diff_active, diff_crop, diff_fapar, diff_evi, diff_ndvi) %>% print(n=500)
# ex_sum %>% filter(sample_id == 59) %>% select(sample_id, sif_reliability, local_outlier_score, local_sif_mad, max_d, n_neighbors, sample_sif, neighbor_sif_mean, diff_sif, diff_active, diff_crop, diff_fapar, diff_evi, diff_ndvi)
