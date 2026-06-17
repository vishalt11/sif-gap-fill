library(tidyverse)
library(sf)
library(terra)

wasp_root <- "data/geodes_wasp_zips/32UQV/2019"
product_id <- "SENTINEL2X_20190715-000000-000_L3A_T32UQV_C_V4-0"
product_dir <- file.path(wasp_root, product_id, product_id)

wasp_file <- function(product_dir, product_id, type, band_or_res = NULL) {
  filename <- if (is.null(band_or_res)) {
    paste0(product_id, "_", type, ".tif")
  } else {
    paste0(product_id, "_", type, "_", band_or_res, ".tif")
  }
  
  subdir <- if (type %in% c("FLG", "WGT", "DTS")) "MASKS" else ""
  file.path(product_dir, subdir, filename)
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

make_sif_polygons <- function(sif_df) {
  sif_polygon_geometry <- pmap(
    list(
      sif_df$Lon_corner1, sif_df$Lat_corner1,
      sif_df$Lon_corner2, sif_df$Lat_corner2,
      sif_df$Lon_corner3, sif_df$Lat_corner3,
      sif_df$Lon_corner4, sif_df$Lat_corner4
    ),
    make_sif_polygon
  )
  
  st_sf(
    sif_df,
    geometry = do.call(st_sfc, c(sif_polygon_geometry, list(crs = 4326)))
  ) %>%
    st_make_valid()
}
reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4

flg_r1 <- rast(wasp_file(product_dir, product_id, "FLG", "R1"))
b4 <- rast(wasp_file(product_dir, product_id, "FRC", "B4"))
#b8 <- rast(wasp_file(product_dir, product_id, "FRC", "B8"))

plot(flg_r1)

ba <- readRDS("data/ba_sif_mgrs_crop_composition.rds") %>% st_drop_geometry()
ba <- ba %>% filter(mgrs_tile == '32UQV')
ba <- ba %>% filter(format(Delta_Date, '%Y-%m') == '2019-07')
nrow(ba)

target_crs <- crs(flg_r1)

sif_polygons <- make_sif_polygons(ba)

sif_polygons_utm <- sif_polygons %>% st_transform(target_crs)

b4 <- crop(b4, vect(sif_polygons_utm))

flg_crop <- crop(flg_r1, b4)

b4 <- ifel(
  flg_crop == land_flag_value & b4 != reflectance_nodata,
  b4 / reflectance_quantification_value,
  NA
)

b4_plot <- mask(b4, vect(sif_polygons_utm))
plot(b4)
lines(sif_polygons_utm)

#mgrs <- readRDS('data/mgrs_de.rds')











# band mean plotting with sif values

ba <- readRDS("data/ba_sif_32UQV_2019_wasp_10m_band_means.rds") %>% st_drop_geometry()

head(ba, 1)

# How the mean b2-b8 values change depending on the sif value
# SIF on the x-axis and each Sentinel-2 band mean on the y-axis, then facet by band. That shows how B2/B3/B4/B8 change as SIF changes.

ba_band_long <- ba %>%
  select(sif = Daily_SIF_740nm, mean_b2, mean_b3, mean_b4, mean_b8) %>%
  pivot_longer(cols = c(mean_b2, mean_b3, mean_b4, mean_b8), names_to = "band",values_to = "reflectance") %>%
  filter(!is.na(sif), !is.na(reflectance))

ggplot(ba_band_long, aes(x = sif, y = reflectance, color = band)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  facet_wrap(~ band, scales = "free_y") +
  scale_color_manual(
    values = c(
      mean_b2 = "#0072B2",
      mean_b3 = "#009E73",
      mean_b4 = "#D55E00",
      mean_b8 = "#CC79A7"
    ),
    labels = c(
      mean_b2 = "B2 blue",
      mean_b3 = "B3 green",
      mean_b4 = "B4 red",
      mean_b8 = "B8 NIR"
    )
  ) +
  labs(
    x = "SIF 740 nm",
    y = "Mean Sentinel-2 reflectance",
    color = "Band"
  ) +
  theme_minimal()








# EDA plotting
# Sanity check mgrs_de and sen2 band extent lineup

germany_states <- giscoR::gisco_get_nuts(
  country = "DE",
  nuts_level = 1,
  resolution = "01",
  epsg = 4326
) %>%
  select(state = NUTS_NAME, geometry)

# Raster CRS
b4_crs <- crs(b4)

# Transform Germany states to B4 CRS
germany_states_utm <- st_transform(germany_states, b4_crs)

# If needed, make/confirm B4 extent polygon
b4_extent_sf <- as.polygons(ext(b4), crs = crs(b4)) %>%
  st_as_sf()

# Make sure MGRS tile polygon is also in B4 CRS
mgrs_32uqv_utm <- st_transform(mgrs_32uqv, b4_crs)

plot(
  st_geometry(germany_states_utm),
  border = "grey45",
  col = "grey95",
  main = "Germany States with B4 Extent and MGRS 32UQV"
)

plot(
  st_geometry(b4_extent_sf),
  add = TRUE,
  border = "red",
  col = NA,
  lwd = 2
)

plot(
  st_geometry(mgrs_32uqv_utm),
  add = TRUE,
  border = "cyan",
  col = NA,
  lwd = 2
)

legend(
  "bottomleft",
  legend = c("Germany states", "B4 raster extent", "mgrs_de 32UQV"),
  col = c("grey45", "red", "cyan"),
  lwd = c(1, 2, 2),
  bg = "white"
)

