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
