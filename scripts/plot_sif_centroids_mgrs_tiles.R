library(tidyverse)
library(sf)
library(terra)

sif_path <- "data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.csv"
mgrs_tif_dir <- "data/temp_data/mgrs_tifs"
germany_path <- "data/germany_boundaries.gpkg"
output_path <- "eda_images/pdf/sif_centroids_mgrs_tile_boundaries.pdf"
satellite_svg_path <- "eda_images/svg/sif_footprints_satellite_2022-05-07.svg"

sif_corner_columns <- c(
  "Lon_corner1", "Lat_corner1",
  "Lon_corner2", "Lat_corner2",
  "Lon_corner3", "Lat_corner3",
  "Lon_corner4", "Lat_corner4"
)

sif <- readr::read_csv(sif_path, show_col_types = FALSE) %>%
  mutate(
    Delta_Date = as.Date(Delta_Date),
    across(all_of(sif_corner_columns), as.numeric)
  ) %>%
  filter(
    !is.na(Delta_Date),
    if_all(all_of(sif_corner_columns), is.finite)
  )
sif_mgrs_tiles <- paste0("T", str_remove(unique(sif$mgrs_tile), "^T"))
germany_states <- st_read(germany_path, quiet = TRUE) %>% st_transform(4326)

b5_paths <- list.files(mgrs_tif_dir, pattern = "_B5\\.tif$", full.names = TRUE)

mgrs_tiles <- map(b5_paths, \(path) { raster <- rast(path); st_as_sf(as.polygons(ext(raster), crs = crs(raster))) %>% mutate(mgrs_tile = str_extract(basename(path), "T\\d{2}[A-Z]{3}")) }) %>%
  bind_rows() %>%
  distinct(mgrs_tile, .keep_all = TRUE) %>%
  filter(mgrs_tile %in% sif_mgrs_tiles) %>%
  st_transform(4326)

sif_polygons <- lapply(
  seq_len(nrow(sif)),
  \(i) st_polygon(list(matrix(
    c(
      sif$Lon_corner1[i], sif$Lat_corner1[i],
      sif$Lon_corner2[i], sif$Lat_corner2[i],
      sif$Lon_corner3[i], sif$Lat_corner3[i],
      sif$Lon_corner4[i], sif$Lat_corner4[i],
      sif$Lon_corner1[i], sif$Lat_corner1[i]
    ),
    ncol = 2,
    byrow = TRUE
  )))
)

sif_footprints <- st_sf(
  sif,
  geometry = st_sfc(sif_polygons, crs = 4326)
) %>%
  st_make_valid()

#-------------------------------------------------------------------------------
# Annotation-free satellite view of the complete footprint track for one date.

selected_date <- as.Date("2022-05-07")

date_footprints_3035 <- sif_footprints %>%
  filter(
    Delta_Date == selected_date,
    is.finite(target_modis_sif)
  ) %>%
  st_transform(3035)

if (nrow(date_footprints_3035) == 0) {
  stop("No finite SIF footprints were found for 2022-05-07.")
}

map_extent_3035 <- date_footprints_3035 %>%
  st_union() %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(5000)

map_extent_4326 <- st_transform(map_extent_3035, 4326)
satellite_tiles_track <- maptiles::get_tiles(
  map_extent_4326,
  provider = "Esri.WorldImagery",
  zoom = 11,
  crop = TRUE,
  project = FALSE
)

satellite_crs <- st_crs(terra::crs(satellite_tiles_track))
map_extent_satellite <- st_transform(map_extent_3035, satellite_crs)
date_footprints_satellite <- st_transform(
  date_footprints_3035,
  satellite_crs
)
map_bbox_satellite <- st_bbox(map_extent_satellite)
sif_limits <- range(
  date_footprints_satellite$target_modis_sif,
  na.rm = TRUE
)

sif_satellite_plot <- ggplot() +
  tidyterra::geom_spatraster_rgb(data = satellite_tiles_track) +
  geom_sf(
    data = date_footprints_satellite,
    aes(fill = target_modis_sif),
    color = NA,
    alpha = 0.70
  ) +
  geom_sf(
    data = date_footprints_satellite,
    aes(color = target_modis_sif),
    fill = NA,
    linewidth = 0.3
  ) +
  scale_fill_viridis_c(limits = sif_limits, guide = "none") +
  scale_color_viridis_c(limits = sif_limits, guide = "none") +
  coord_sf(
    xlim = c(map_bbox_satellite["xmin"], map_bbox_satellite["xmax"]),
    ylim = c(map_bbox_satellite["ymin"], map_bbox_satellite["ymax"]),
    crs = satellite_crs,
    datum = NA,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

dir.create(dirname(satellite_svg_path), recursive = TRUE, showWarnings = FALSE)
ggsave(
  satellite_svg_path,
  plot = sif_satellite_plot,
  device = "svg",
  width = 7,
  height = 9,
  units = "in",
  bg = "transparent"
)

message("Selected SIF date: ", selected_date)
message("Footprints plotted: ", nrow(date_footprints_3035))
message("Saved SVG: ", satellite_svg_path)

sif_centroids <- sif_footprints %>%
  st_transform(3035) %>%
  st_centroid() %>%
  st_transform(4326)

sif_centroids <- st_filter(sif_centroids, st_union(mgrs_tiles), .predicate = st_within)

sif_centroids <- sif_centroids %>%
  mutate(
    sif_year = lubridate::year(as.Date(Delta_Date)),
    sample_order = runif(n())
  ) %>%
  group_by(sif_year) %>%
  arrange(sample_order, .by_group = TRUE) %>%
  slice_head(n = 2000) %>%
  ungroup() %>%
  select(-sample_order)

table(sif_centroids$sif_year)

mgrs_labels <- mgrs_tiles %>%
  st_transform(3035) %>%
  st_centroid() %>%
  st_transform(4326) %>%
  mutate(mgrs_tile = str_remove(mgrs_tile, "^T"))

sif_plot <- ggplot() +
  geom_sf(data = germany_states, fill = "grey96", color = "grey60", linewidth = 0.25) +
  geom_sf(data = sif_centroids, aes(color = target_modis_sif), size = 0.25, alpha = 0.8) +
  geom_sf(data = mgrs_tiles, fill = NA, color = "grey20", linewidth = 0.35) +
  geom_sf_text(data = mgrs_labels, aes(label = mgrs_tile), size = 2.5, color = "grey20") +
  scale_color_viridis_c(name = "target_modis_sif") +
  coord_sf(expand = FALSE) +
  labs(title = "SIF centroids inside selected MGRS tile bounding boxes", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3), legend.position = "right", plot.title = element_text(face = "bold"))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, plot = sif_plot, width = 14, height = 10, units = "in")

#-------------------------------------------------------------------------------

ctr <- rast("data/crop_type_tif/croptypes_2024.tif")
crop_classes <- readr::read_delim("data/crop_type_tif/LEGEND_CropTypes.txt", delim = "\t", show_col_types = FALSE)
colnames(crop_classes) <- c("code", "label")
levels(ctr) <- data.frame(value = crop_classes$code, crop = crop_classes$label)
germany_states_crop <- st_transform(germany_states, st_crs(crs(ctr)))

crop_1km <- aggregate(ctr > 0, fact = 100, fun = "max", na.rm = TRUE)
wheat_1km <- aggregate(ctr == 11, fact = 100, fun = "mean", na.rm = TRUE)
crop_map_1km <- ifel(wheat_1km >= 0.4, 2, ifel(crop_1km == 1, 1, NA)) %>% mask(vect(st_union(germany_states_crop)))
names(crop_map_1km) <- "crop_group"

high_wheat <- wheat_1km >= 0.4
plot(high_wheat)

crop_map_df <- as.data.frame(crop_map_1km, xy = TRUE, na.rm = TRUE) %>%
  mutate(crop_group = factor(crop_group, levels = c(1, 2), labels = c("All crop types", "Winter wheat >= 40%")))

crop_map_plot <- ggplot() +
  geom_sf(data = germany_states_crop, fill = "black", color = "grey60", linewidth = 0.3) +
  geom_raster(data = crop_map_df, aes(x = x, y = y, fill = crop_group)) +
  geom_sf(data = germany_states_crop, fill = NA, color = "grey60", linewidth = 0.3) +
  scale_fill_manual(values = c("All crop types" = "forestgreen", "Winter wheat >= 40%" = "gold"), name = NULL) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = "Crop types and concentrated winter wheat in Germany, 2024", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave("eda_images/pdf/crop_types_winter_wheat_2024_1km.pdf", plot = crop_map_plot, width = 11, height = 10, units = "in")

#-------------------------------------------------------------------------------

wheat_fraction_1km <- mask(wheat_1km, vect(st_union(germany_states_crop)))
names(wheat_fraction_1km) <- "wheat_fraction"
wheat_fraction_df <- as.data.frame(wheat_fraction_1km, xy = TRUE, na.rm = TRUE)

wheat_fraction_plot <- ggplot() +
  geom_raster(data = wheat_fraction_df, aes(x = x, y = y, fill = wheat_fraction)) +
  geom_sf(data = germany_states_crop, fill = NA, color = "grey90", linewidth = 0.3) +
  scale_fill_viridis_c(name = "Winter wheat fraction", limits = c(0, 1)) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(title = "Winter wheat fraction at 1 km resolution, 2024", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave("eda_images/pdf/winter_wheat_fraction_2024_1km.pdf", plot = wheat_fraction_plot, width = 11, height = 10, units = "in")

#-------------------------------------------------------------------------------

phase_data <- sif %>% filter(Delta_Date == as.Date("2022-05-07"))
phase_polygons <- lapply(seq_len(nrow(phase_data)), \(i) st_polygon(list(matrix(c(phase_data$Lon_corner1[i], phase_data$Lat_corner1[i], phase_data$Lon_corner2[i], phase_data$Lat_corner2[i], phase_data$Lon_corner3[i], phase_data$Lat_corner3[i], phase_data$Lon_corner4[i], phase_data$Lat_corner4[i], phase_data$Lon_corner1[i], phase_data$Lat_corner1[i]), ncol = 2, byrow = TRUE))))
phase_polygons <- st_sf(phase_data, geometry = st_sfc(phase_polygons, crs = 4326))

track_extent <- phase_polygons %>%
  st_union() %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(3035) %>%
  st_buffer(5000) %>%
  st_transform(4326)

satellite_tiles <- maptiles::get_tiles(track_extent, provider = "Esri.WorldImagery", zoom = 11, crop = TRUE, project = FALSE)
phase_polygons_3857 <- phase_polygons %>% st_transform(3035) %>% st_buffer(-20) %>% st_transform(3857)
track_bbox_3857 <- st_bbox(st_transform(track_extent, 3857))

phase_angle_plot <- ggplot() +
  tidyterra::geom_spatraster_rgb(data = satellite_tiles) +
  geom_sf(data = phase_polygons_3857, aes(fill = phase_angle, color = phase_angle), alpha = 0.9, linewidth = 0.2) +
  scale_fill_viridis_c(name = "Phase angle") +
  scale_color_viridis_c(guide = "none") +
  coord_sf(xlim = c(track_bbox_3857["xmin"], track_bbox_3857["xmax"]), ylim = c(track_bbox_3857["ymin"], track_bbox_3857["ymax"]), crs = st_crs(3857), datum = NA, expand = FALSE) +
  labs(title = "SIF footprint phase angle", x = NULL, y = NULL, caption = "Basemap: Esri World Imagery") +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

sif_value_plot <- ggplot() +
  tidyterra::geom_spatraster_rgb(data = satellite_tiles) +
  geom_sf(data = phase_polygons_3857, aes(fill = target_modis_sif, color = target_modis_sif), alpha = 0.9, linewidth = 0.2) +
  scale_fill_viridis_c(name = "SIF") +
  scale_color_viridis_c(guide = "none") +
  coord_sf(xlim = c(track_bbox_3857["xmin"], track_bbox_3857["xmax"]), ylim = c(track_bbox_3857["ymin"], track_bbox_3857["ymax"]), crs = st_crs(3857), datum = NA, expand = FALSE) +
  labs(title = "SIF footprint value", x = NULL, y = NULL, caption = "Basemap: Esri World Imagery") +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

phase_sif_plot <- patchwork::wrap_plots(phase_angle_plot, sif_value_plot, ncol = 2)
ggsave("eda_images/pdf/sif_phase_angle_and_value_satellite_2022-05-07.pdf", plot = phase_sif_plot, width = 16, height = 8, units = "in")

