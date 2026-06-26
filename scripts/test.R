library(tidyverse)
library(sf)
library(terra)

ba <- readRDS("data/spectral_indices_means/ba_sif_32UNA_wasp_spectral_indices.rds") %>% st_drop_geometry()
head(ba, 1)

summary(ba$crop_pixel_count)
plot(density(ba$crop_pixel_count))

ba_band_long <- ba %>%
  select(sif = Daily_SIF_740nm, mean_ndvi, 
         mean_ndwi, mean_ndre8a, mean_nirv,
         mean_osavi, mean_psri, mean_ndre) %>%
  pivot_longer(cols = c(mean_ndvi, mean_ndwi, mean_ndre8a, mean_nirv, mean_osavi, mean_psri, mean_ndre), names_to = "band",values_to = "reflectance") %>%
  filter(!is.na(sif), !is.na(reflectance))

ggplot(ba_band_long, aes(x = sif, y = reflectance, color = band)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  facet_wrap(~ band, scales = "free_y") +
  scale_color_manual(
    values = c(
      mean_ndvi = "#0072B2",
      mean_ndwi = "#009E73",
      mean_ndre8a = "#D55E00",
      mean_nirv = "#CC79A7",
      mean_osavi = "red", 
      mean_psri = "blue", 
      mean_ndre = "green"
    ),
    labels = c(
      mean_ndvi = "ndvi",
      mean_ndwi = "ndwi",
      mean_ndre8a = "ndre8a",
      mean_nirv = "nirv",
      mean_osavi = "osavi", 
      mean_psri = "psri", 
      mean_ndre = "ndre"
    )
  ) +
  labs(
    x = "SIF 740 nm",
    y = "Mean Sentinel-2 reflectance",
    color = "Band"
  ) +
  theme_minimal()

ba_ndvi_2022 <- ba %>%
  filter(ww_pct >= 0.3) %>%
  mutate(
    Delta_Date = as.Date(Delta_Date),
    month = lubridate::month(Delta_Date, label = TRUE, abbr = TRUE)
  ) %>%
  filter(lubridate::year(Delta_Date) == 2022) %>%
  select(
    sif = Daily_SIF_740nm,
    mean_ndvi,
    crop_pixel_count,
    month
  ) %>%
  filter(!is.na(sif), !is.na(mean_ndvi), !is.na(crop_pixel_count), !is.na(month))

ggplot(ba_ndvi_2022, aes(x = mean_ndvi, y = sif)) +
  geom_point(aes(color = crop_pixel_count, shape = month), alpha = 0.85, size = 3) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "black") +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "Mean NDVI",
    y = "SIF 740 nm",
    color = "Crop pixel count",
    shape = "Month"
  ) +
  theme_minimal()

ggplot(ba_ndvi_2022, aes(x = mean_ndvi, y = sif)) +
  geom_point(aes(color = crop_pixel_count), alpha = 0.85, size = 3) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "black") +
  facet_wrap(~ month, scales = "free") +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "Mean NDVI",
    y = "SIF 740 nm",
    color = "Crop pixel count"
  ) +
  theme_minimal()


#-------------------------------------------------------------------------------
# area to point sif row model data inspection
library(tidyverse)
library(arrow)

polygon_targets_dir <- file.path("data", "backup", "test")

polygon_targets_files <- list.files(
  polygon_targets_dir,
  pattern = "\\.parquet$",
  full.names = TRUE,
  recursive = FALSE
)



df <- polygon_targets_files |>
  set_names(tools::file_path_sans_ext(basename(polygon_targets_files))) |>
  map_dfr(
    ~ read_parquet(.x) |>
      mutate(polygon_targets_file = .x),
    .id = "polygon_targets_source"
  ) |>
  mutate(
    Delta_Date = as.Date(Delta_Date),
    sif = Daily_SIF_740nm,
    #ww_pct = replace_na(ww_pct, 0),
    #crop_pixel_count = replace_na(crop_pixel_count, 0)
  )

colnames(df)


llif (!"month" %in% names(df)) {
  df <- df |> mutate(month = as.integer(format(Delta_Date, "%m")))
}

if (!"year" %in% names(df)) {
  df <- df |> mutate(year = as.integer(format(Delta_Date, "%Y")))
}

crop_count_cols <- names(df)[str_detect(names(df), "^crop_count_")]

eda_dir <- file.path("eda_images", "area_to_point_sif_row_filters")
dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)

q_fun <- function(x, p) {
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

sif_q1 <- q_fun(df$sif, 0.25)
sif_q3 <- q_fun(df$sif, 0.75)
sif_iqr <- sif_q3 - sif_q1
sif_iqr_lower <- sif_q1 - 1.5 * sif_iqr
sif_iqr_upper <- sif_q3 + 1.5 * sif_iqr
sif_median <- median(df$sif, na.rm = TRUE)
sif_mad <- mad(df$sif, constant = 1.4826, na.rm = TRUE)
sif_mad_lower <- sif_median - 3 * sif_mad
sif_mad_upper <- sif_median + 3 * sif_mad
sif_p001 <- q_fun(df$sif, 0.001)
sif_p999 <- q_fun(df$sif, 0.999)

df_eda <- df |>
  mutate(
    sif_negative = sif < 0,
    sif_iqr_outlier = sif < sif_iqr_lower | sif > sif_iqr_upper,
    sif_mad_outlier = sif < sif_mad_lower | sif > sif_mad_upper,
    sif_extreme_quantile = sif < sif_p001 | sif > sif_p999,
    crop_count_bin = cut(
      crop_pixel_count,
      breaks = c(-Inf, 0, 10, 100, 500, 1000, 5000, 10000, Inf),
      labels = c("0", "1-10", "11-100", "101-500", "501-1000", "1001-5000", "5001-10000", ">10000")
    ),
    ww_pct_bin = cut(
      ww_pct,
      breaks = c(-Inf, 0, 0.05, 0.2, 0.5, 0.7, Inf),
      labels = c("0", "(0,0.05]", "(0.05,0.2]", "(0.2,0.5]", "(0.5,0.7]", ">0.7")
    )
  )

cat("\n# 1. Overall row counts\n")
overall_counts <- df_eda |>
  summarise(
    n_rows = n(),
    n_missing_sif = sum(is.na(sif)),
    n_negative_sif = sum(sif_negative, na.rm = TRUE),
    pct_negative_sif = mean(sif_negative, na.rm = TRUE),
    n_crop_pixel_count_zero = sum(crop_pixel_count == 0, na.rm = TRUE),
    pct_crop_pixel_count_zero = mean(crop_pixel_count == 0, na.rm = TRUE),
    n_ww_pct_zero = sum(ww_pct == 0, na.rm = TRUE),
    pct_ww_pct_zero = mean(ww_pct == 0, na.rm = TRUE),
    n_ww_pct_gt_0 = sum(ww_pct > 0, na.rm = TRUE),
    n_ww_pct_ge_0_2 = sum(ww_pct >= 0.2, na.rm = TRUE),
    n_ww_pct_ge_0_5 = sum(ww_pct >= 0.5, na.rm = TRUE)
  )
print(overall_counts)

cat("\n# 2. SIF distribution and outlier thresholds\n")
sif_distribution <- df_eda |>
  summarise(
    n = sum(!is.na(sif)),
    min = min(sif, na.rm = TRUE),
    p001 = q_fun(sif, 0.001),
    p01 = q_fun(sif, 0.01),
    p05 = q_fun(sif, 0.05),
    p25 = q_fun(sif, 0.25),
    median = median(sif, na.rm = TRUE),
    mean = mean(sif, na.rm = TRUE),
    p75 = q_fun(sif, 0.75),
    p95 = q_fun(sif, 0.95),
    p99 = q_fun(sif, 0.99),
    p999 = q_fun(sif, 0.999),
    max = max(sif, na.rm = TRUE),
    sd = sd(sif, na.rm = TRUE),
    mad = mad(sif, constant = 1.4826, na.rm = TRUE),
    iqr_lower = sif_iqr_lower,
    iqr_upper = sif_iqr_upper,
    mad_lower = sif_mad_lower,
    mad_upper = sif_mad_upper,
    n_iqr_outliers = sum(sif_iqr_outlier, na.rm = TRUE),
    n_mad_outliers = sum(sif_mad_outlier, na.rm = TRUE),
    n_extreme_quantile = sum(sif_extreme_quantile, na.rm = TRUE)
  )
print(sif_distribution)

cat("\n# 3. Most negative SIF rows\n")
lowest_sif_rows <- df_eda |>
  arrange(sif) |>
  select(
    sif_extract_id, sif_id, sif, Delta_Date, year_month, mgrs_tile,
    Latitude, Longitude, Quality_Flag, Metadata.MeasurementMode,
    crop_pixel_count, ww_pct, all_of(crop_count_cols)
  ) |>
  head(30)
print(lowest_sif_rows, n = 30, width = Inf)

cat("\n# 4. Highest SIF rows\n")
highest_sif_rows <- df_eda |>
  arrange(desc(sif)) |>
  select(
    sif_extract_id, sif_id, sif, Delta_Date, year_month, mgrs_tile,
    Latitude, Longitude, Quality_Flag, Metadata.MeasurementMode,
    crop_pixel_count, ww_pct, all_of(crop_count_cols)
  ) |>
  head(30)
print(highest_sif_rows, n = 30, width = Inf)

cat("\n# 5. SIF summary by crop_pixel_count bin\n")
sif_by_crop_count_bin <- df_eda |>
  group_by(crop_count_bin) |>
  summarise(
    n = n(),
    pct_rows = n / nrow(df_eda),
    n_negative = sum(sif < 0, na.rm = TRUE),
    pct_negative = mean(sif < 0, na.rm = TRUE),
    sif_min = min(sif, na.rm = TRUE),
    sif_p05 = q_fun(sif, 0.05),
    sif_median = median(sif, na.rm = TRUE),
    sif_mean = mean(sif, na.rm = TRUE),
    sif_p95 = q_fun(sif, 0.95),
    sif_max = max(sif, na.rm = TRUE),
    ww_pct_median = median(ww_pct, na.rm = TRUE),
    ww_pct_mean = mean(ww_pct, na.rm = TRUE),
    .groups = "drop"
  )
print(sif_by_crop_count_bin, n = Inf, width = Inf)

cat("\n# 6. SIF summary by winter-wheat fraction bin\n")
sif_by_ww_pct_bin <- df_eda |>
  group_by(ww_pct_bin) |>
  summarise(
    n = n(),
    pct_rows = n / nrow(df_eda),
    n_negative = sum(sif < 0, na.rm = TRUE),
    pct_negative = mean(sif < 0, na.rm = TRUE),
    crop_pixel_count_median = median(crop_pixel_count, na.rm = TRUE),
    sif_min = min(sif, na.rm = TRUE),
    sif_p05 = q_fun(sif, 0.05),
    sif_median = median(sif, na.rm = TRUE),
    sif_mean = mean(sif, na.rm = TRUE),
    sif_p95 = q_fun(sif, 0.95),
    sif_max = max(sif, na.rm = TRUE),
    .groups = "drop"
  )
print(sif_by_ww_pct_bin, n = Inf, width = Inf)

cat("\n# 7. SIF summary for candidate filter sets\n")
candidate_filter_summary <- tibble(
  filter_name = c(
    "all_rows",
    "quality_flag_0",
    "measurement_mode_0_1",
    "non_negative_sif",
    "crop_pixel_count_gt_0",
    "ww_pct_gt_0",
    "ww_pct_ge_0_2",
    "ww_pct_ge_0_5",
    "quality0_crop_gt0",
    "quality0_crop_gt0_ww_gt0",
    "quality0_crop_gt0_ww_ge_0_2",
    "quality0_crop_gt0_ww_ge_0_5",
    "quality0_crop_gt0_keep_negative_non_extreme"
  ),
  keep = list(
    rep(TRUE, nrow(df_eda)),
    df_eda$Quality_Flag == 0,
    df_eda$Metadata.MeasurementMode %in% c(0, 1),
    df_eda$sif >= 0,
    df_eda$crop_pixel_count > 0,
    df_eda$ww_pct > 0,
    df_eda$ww_pct >= 0.2,
    df_eda$ww_pct >= 0.5,
    df_eda$Quality_Flag == 0 & df_eda$crop_pixel_count > 0,
    df_eda$Quality_Flag == 0 & df_eda$crop_pixel_count > 0 & df_eda$ww_pct > 0,
    df_eda$Quality_Flag == 0 & df_eda$crop_pixel_count > 0 & df_eda$ww_pct >= 0.2,
    df_eda$Quality_Flag == 0 & df_eda$crop_pixel_count > 0 & df_eda$ww_pct >= 0.5,
    df_eda$Quality_Flag == 0 & df_eda$crop_pixel_count > 0 & !df_eda$sif_extreme_quantile
  )
) |>
  mutate(
    n = map_int(keep, ~ sum(.x, na.rm = TRUE)),
    pct_rows = n / nrow(df_eda),
    sif_mean = map_dbl(keep, ~ mean(df_eda$sif[.x], na.rm = TRUE)),
    sif_median = map_dbl(keep, ~ median(df_eda$sif[.x], na.rm = TRUE)),
    sif_min = map_dbl(keep, ~ min(df_eda$sif[.x], na.rm = TRUE)),
    sif_p05 = map_dbl(keep, ~ q_fun(df_eda$sif[.x], 0.05)),
    sif_p95 = map_dbl(keep, ~ q_fun(df_eda$sif[.x], 0.95)),
    sif_max = map_dbl(keep, ~ max(df_eda$sif[.x], na.rm = TRUE)),
    pct_negative = map_dbl(keep, ~ mean(df_eda$sif[.x] < 0, na.rm = TRUE)),
    crop_pixel_count_median = map_dbl(keep, ~ median(df_eda$crop_pixel_count[.x], na.rm = TRUE)),
    ww_pct_median = map_dbl(keep, ~ median(df_eda$ww_pct[.x], na.rm = TRUE))
  ) |>
  select(-keep)
print(candidate_filter_summary, n = Inf, width = Inf)

cat("\n# 8. SIF summary by quality flag and measurement mode\n")
sif_by_quality_mode <- df_eda |>
  group_by(Quality_Flag, Metadata.MeasurementMode) |>
  summarise(
    n = n(),
    pct_rows = n / nrow(df_eda),
    pct_negative = mean(sif < 0, na.rm = TRUE),
    sif_median = median(sif, na.rm = TRUE),
    sif_mean = mean(sif, na.rm = TRUE),
    sif_p05 = q_fun(sif, 0.05),
    sif_p95 = q_fun(sif, 0.95),
    crop_pixel_count_median = median(crop_pixel_count, na.rm = TRUE),
    ww_pct_median = median(ww_pct, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Quality_Flag, Metadata.MeasurementMode)
print(sif_by_quality_mode, n = Inf, width = Inf)

cat("\n# 9. SIF summary by month\n")
sif_by_month <- df_eda |>
  group_by(month) |>
  summarise(
    n = n(),
    pct_negative = mean(sif < 0, na.rm = TRUE),
    sif_median = median(sif, na.rm = TRUE),
    sif_mean = mean(sif, na.rm = TRUE),
    sif_p05 = q_fun(sif, 0.05),
    sif_p95 = q_fun(sif, 0.95),
    crop_pixel_count_median = median(crop_pixel_count, na.rm = TRUE),
    ww_pct_median = median(ww_pct, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(month)
print(sif_by_month, n = Inf, width = Inf)

cat("\n# 10. Dominant crop-count class by row\n")
dominant_crop_summary <- df_eda |>
  select(sif_extract_id, sif, all_of(crop_count_cols)) |>
  pivot_longer(
    cols = all_of(crop_count_cols),
    names_to = "crop_class",
    values_to = "crop_count"
  ) |>
  group_by(sif_extract_id) |>
  slice_max(crop_count, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(crop_class = str_remove(crop_class, "^crop_count_")) |>
  group_by(crop_class) |>
  summarise(
    n = n(),
    sif_median = median(sif, na.rm = TRUE),
    sif_mean = mean(sif, na.rm = TRUE),
    pct_negative = mean(sif < 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n))
print(dominant_crop_summary, n = Inf, width = Inf)

p_sif_hist <- ggplot(df_eda, aes(x = sif)) +
  geom_histogram(bins = 120) +
  geom_vline(xintercept = c(sif_iqr_lower, sif_iqr_upper), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(sif_mad_lower, sif_mad_upper), color = "blue", linetype = "dotted") +
  labs(
    x = "Daily SIF 740 nm",
    y = "Number of SIF rows",
    title = "SIF distribution with IQR and MAD outlier thresholds"
  ) +
  theme_minimal()

p_crop_sif <- ggplot(df_eda, aes(x = crop_pixel_count, y = sif)) +
  geom_point(aes(color = ww_pct), alpha = 0.35, size = 0.8) +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "Crop pixel count",
    y = "Daily SIF 740 nm",
    color = "Winter wheat fraction",
    title = "SIF versus crop pixel count"
  ) +
  theme_minimal()

p_ww_sif <- ggplot(df_eda, aes(x = ww_pct, y = sif)) +
  geom_point(aes(color = crop_pixel_count), alpha = 0.35, size = 0.8) +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "Winter wheat fraction",
    y = "Daily SIF 740 nm",
    color = "Crop pixel count",
    title = "SIF versus winter wheat fraction"
  ) +
  theme_minimal()

p_month_sif <- ggplot(df_eda, aes(x = factor(month), y = sif)) +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(
    x = "Month",
    y = "Daily SIF 740 nm",
    title = "SIF distribution by month"
  ) +
  theme_minimal()

ggsave(file.path(eda_dir, "sif_hist_outlier_thresholds.png"), p_sif_hist, width = 9, height = 5, dpi = 200)
ggsave(file.path(eda_dir, "sif_vs_crop_pixel_count.png"), p_crop_sif, width = 8, height = 6, dpi = 200)
ggsave(file.path(eda_dir, "sif_vs_ww_pct.png"), p_ww_sif, width = 8, height = 6, dpi = 200)
ggsave(file.path(eda_dir, "sif_by_month_boxplot.png"), p_month_sif, width = 8, height = 5, dpi = 200)

cat("\nSaved EDA plots to: ", eda_dir, "\n", sep = "")

#-------------------------------------------------------------------------------
library(terra)
library(tidyverse)
library(sf)

df <- readRDS('data/model_data_wpar.rds')
mgrs <- readRDS('data/mgrs_de.rds')

df_tiles <- unique(df$mgrs_tile)

mgrs_df_tiles <- mgrs %>%
  filter(mgrs_tile %in% df_tiles) %>%
  st_transform(4326)

germany_states <- giscoR::gisco_get_nuts(
  country = "DE",
  nuts_level = 1,
  resolution = "03",
  epsg = 4326
) %>%
  st_make_valid()
germany_bbox <- st_bbox(germany_states)

ggplot() +
  geom_sf(
    data = germany_states,
    fill = "grey95",
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_sf(
    data = mgrs_df_tiles,
    aes(fill = mgrs_tile),
    color = "black",
    linewidth = 0.45,
    alpha = 0.28
  ) +
  geom_sf(
    data = df,
    fill = NA,
    color = "red",
    linewidth = 0.08,
    alpha = 0.25
  ) +
  geom_sf_text(
    data = mgrs_df_tiles,
    aes(label = mgrs_tile),
    size = 3,
    fontface = "bold"
  ) +
  coord_sf(
    xlim = germany_bbox[c("xmin", "xmax")],
    ylim = germany_bbox[c("ymin", "ymax")],
    expand = FALSE
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "SIF Polygons and MGRS Tiles",
    fill = "MGRS tile",
    x = NULL,
    y = NULL
  )



