library(tidyverse)
library(sf)
library(giscoR)
library(terra)
library(raster)
library(leaflet)
library(htmlwidgets)

set.seed(42)
corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

required_cols <- c("Latitude", "Longitude", corner_cols)
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

sif_path <- "data/sif_sf_months2_7_cleaned.rds"

sif_raw <- readRDS(sif_path)

sif_df <- if (inherits(sif_raw, "sf")) {
  st_drop_geometry(sif_raw)
} else {
  as_tibble(sif_raw)
}

has_target_col <- target_col %in% names(sif_df)


missing_cols <- setdiff(required_cols, names(sif_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}



sif_sample <- sif_df %>%
  mutate(
    across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric),
    source_row = row_number()
  ) %>%
  filter(
    !is.na(Latitude),
    !is.na(Longitude),
    if_all(all_of(corner_cols), ~ !is.na(.x))
  ) %>%
  slice_sample(n = 5) %>%
  mutate(
    plot_id = row_number(),
    sif_value = if (has_target_col) as.numeric(.data[[target_col]]) else NA_real_,
    facet_label = if_else(
      is.na(sif_value),
      paste0("sample ", plot_id, " | source row ", source_row),
      paste0("sample ", plot_id, " | SIF ", round(sif_value, 4))
    ),
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

centroid_points <- sif_sample %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326) %>%
  mutate(point_type = "Polygon centroid")

latlon_points <- sif_sample %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  mutate(point_type = "Latitude/Longitude")

corner_mean_points <- sif_sample %>%
  st_drop_geometry() %>%
  mutate(
    corner_mean_lon = rowMeans(across(starts_with("Lon_corner"))),
    corner_mean_lat = rowMeans(across(starts_with("Lat_corner")))
  ) %>%
  st_as_sf(coords = c("corner_mean_lon", "corner_mean_lat"), crs = 4326, remove = FALSE) %>%
  mutate(point_type = "Mean of 4 corners")

comparison_points <- bind_rows(
  centroid_points,
  latlon_points,
  corner_mean_points
)

polygon_plot_df <- sif_sample %>%
  st_drop_geometry() %>%
  transmute(
    plot_id,
    facet_label,
    lon = pmap(
      list(Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4),
      function(lon1, lon2, lon3, lon4) c(lon1, lon2, lon3, lon4, lon1)
    ),
    lat = pmap(
      list(Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4),
      function(lat1, lat2, lat3, lat4) c(lat1, lat2, lat3, lat4, lat1)
    )
  ) %>%
  unnest(c(lon, lat))

point_plot_df <- comparison_points %>%
  mutate(
    lon = st_coordinates(.)[, "X"],
    lat = st_coordinates(.)[, "Y"]
  ) %>%
  st_drop_geometry()

centroid_offset <- tibble(
  plot_id = centroid_points$plot_id,
  facet_label = centroid_points$facet_label,
  source_row = centroid_points$source_row,
  sif_value = centroid_points$sif_value,
  centroid_lon = st_coordinates(centroid_points)[, "X"],
  centroid_lat = st_coordinates(centroid_points)[, "Y"],
  original_lon = latlon_points$Longitude,
  original_lat = latlon_points$Latitude,
  corner_mean_lon = corner_mean_points$corner_mean_lon,
  corner_mean_lat = corner_mean_points$corner_mean_lat,
  centroid_to_original_m = as.numeric(st_distance(
    st_transform(centroid_points, centroid_crs),
    st_transform(latlon_points, centroid_crs),
    by_element = TRUE
  )),
  original_to_corner_mean_m = as.numeric(st_distance(
    st_transform(latlon_points, centroid_crs),
    st_transform(corner_mean_points, centroid_crs),
    by_element = TRUE
  )),
  centroid_to_corner_mean_m = as.numeric(st_distance(
    st_transform(centroid_points, centroid_crs),
    st_transform(corner_mean_points, centroid_crs),
    by_element = TRUE
  )),
  original_in_polygon = diag(st_covers(
    sif_sample,
    latlon_points,
    sparse = FALSE
  ))
)

offset_segment_df <- centroid_offset %>%
  transmute(
    plot_id,
    facet_label,
    x = original_lon,
    y = original_lat,
    xend = centroid_lon,
    yend = centroid_lat
  )

centroid_check_plot <- ggplot() +
  geom_polygon(
    data = polygon_plot_df,
    aes(x = lon, y = lat, group = plot_id),
    fill = "grey90",
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_segment(
    data = offset_segment_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey35",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_point(
    data = point_plot_df,
    aes(x = lon, y = lat, color = point_type),
    size = 2.8,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      "Polygon centroid" = "blue",
      "Mean of 4 corners" = "purple",
      "Latitude/Longitude" = "red"
    )
  ) +
  facet_wrap(~ facet_label, scales = "free") +
  labs(
    title = "SIF footprint centroid vs Latitude/Longitude point",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    aspect.ratio = 1,
    panel.grid.major = element_line(linewidth = 0.2, color = "grey80"),
    strip.text = element_text(face = "bold")
  )

print(centroid_offset,width = Inf)
print(centroid_check_plot)


#-------------------------------------------------------------------------------
# oco2 revisits every 16 days

germany_states <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 1, resolution = "01", epsg = 4326) %>%
  dplyr::select(state = NUTS_NAME, geometry) %>%
  st_make_valid()

sif_polygons_2019 <- sif_df %>%
  mutate(
    sif_date = as.Date(Delta_Date),
    measurement_mode = as.integer(Metadata.MeasurementMode),
    across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric)
  ) %>%
  filter(
    lubridate::year(sif_date) == 2019,
    measurement_mode %in% c(0L, 1L)
  ) %>%
  mutate(
    geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

sif_centroids_2019 <- sif_polygons_2019 %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326) %>%
  mutate(
    sif_date = as.Date(sif_date),
    measurement_mode = factor(measurement_mode, levels = c(0, 1), labels = c("0 Nadir", "1 Glint"))
  )

sif_centroids_april_2019 <- sif_centroids_2019 %>%
  filter(lubridate::month(sif_date) == 4) %>%
  mutate(sif_date = factor(sif_date))

mode_dates_2019 <- sif_centroids_2019 %>%
  st_drop_geometry() %>%
  count(sif_date, measurement_mode, name = "n_centroids") %>%
  arrange(sif_date, measurement_mode)

oco2_mode_centroids_2019_plot <- ggplot() +
  geom_sf(data = germany_states, fill = "grey95", color = "grey45", linewidth = 0.25) +
  geom_sf(data = sif_centroids_april_2019, aes(color = sif_date, shape = measurement_mode), size = 3.5, alpha = 0.55) +
  scale_shape_manual(values = c("0 Nadir" = 16, "1 Glint" = 17)) +
  coord_sf(xlim = c(5.8, 15.2), ylim = c(47.2, 55.1), expand = FALSE) +
  labs(
    title = "April 2019 OCO-2 SIF polygon centroids over Germany",
    subtitle = "Color = April observation date, shape = measurement mode",
    color = "Date",
    shape = "Mode"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.major = element_line(linewidth = 0.2, color = "grey82"))

print(mode_dates_2019)
print(oco2_mode_centroids_2019_plot)

ggsave(filename = "eda_images/oco2_mode_centroids_2019_april_plot.png", plot = oco2_mode_centroids_2019_plot, width = 4400, height = 5000, units = "px", dpi = 300)


#-------------------------------------------------------------------------------
# plot top 100 sif polys w crop raster

crop_dir <- "data/crop_type_tif"
legend_path <- file.path(crop_dir, "LEGEND_CropTypes.txt")
top_sif_crop_leaflet_html <- "eda_images/top100_sif_crop_type_leaflet.html"

crop_classes <- read.delim(legend_path, sep = "\t", stringsAsFactors = FALSE)
names(crop_classes) <- c("code", "label")
crop_classes$code <- as.integer(crop_classes$code)

legend_classes <- crop_classes[crop_classes$code != 0, ]

fallback_palette <- grDevices::hcl.colors(nrow(legend_classes), "Dark 3")
names(fallback_palette) <- as.character(legend_classes$code)

semantic_palette <- c(
  "11" = "blue",
  "12" = "#c9a227",
  "13" = "#b08935",
  "14" = "#d7b56d",
  "21" = "#ffe08a",
  "22" = "#bf8f2f",
  "23" = "#e6cc98",
  "30" = "#2f9e44",
  "40" = "#20c997",
  "50" = "#8d6e63",
  "60" = "#e64980",
  "71" = "#ffd43b",
  "81" = "#51cf66",
  "82" = "#69db7c",
  "83" = "#2b8a3e",
  "90" = "#9c36b5",
  "100" = "#e03131",
  "110" = "#087f5b",
  "111" = "#adb5bd"
)

crop_palette <- fallback_palette
crop_palette[intersect(names(crop_palette), names(semantic_palette))] <- semantic_palette[intersect(names(crop_palette), names(semantic_palette))]

crop_pal <- leaflet::colorFactor(palette = unname(crop_palette), domain = as.integer(names(crop_palette)), na.color = "transparent")

top_sif_polygons <- sif_df %>%
  mutate(sif_date = as.Date(Delta_Date), crop_year = lubridate::year(sif_date), across(all_of(c("Latitude", "Longitude", corner_cols, target_col)), as.numeric)) %>%
  filter(crop_year >= 2019, crop_year <= 2024) %>%
  slice_max(order_by = .data[[target_col]], n = 100, with_ties = FALSE, na_rm = TRUE) %>%
  arrange(desc(.data[[target_col]])) %>%
  mutate(
    sif_rank = row_number(),
    layer_label = paste0("rank ", stringr::str_pad(sif_rank, 3, pad = "0"), " | SIF ", round(.data[[target_col]], 3), " | ", sif_date),
    geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

crop_raster_cache <- new.env(parent = emptyenv())

get_crop_raster <- function(year) {
  key <- as.character(year)
  if (!exists(key, envir = crop_raster_cache)) {
    crop_raster <- terra::rast(file.path(crop_dir, paste0("croptypes_", year, ".tif")))
    levels(crop_raster) <- data.frame(value = crop_classes$code, crop = crop_classes$label)
    assign(key, crop_raster, envir = crop_raster_cache)
  }
  get(key, envir = crop_raster_cache)
}

top_sif_crop_to_leaflet_raster <- function(i) {
  sif_poly <- top_sif_polygons[i, ]
  crop_raster <- get_crop_raster(sif_poly$crop_year)
  sif_poly_crop_crs <- st_transform(sif_poly, terra::crs(crop_raster))
  crop_sif <- terra::crop(crop_raster, terra::vect(sif_poly_crop_crs), snap = "out")
  crop_sif <- terra::mask(crop_sif, terra::vect(sif_poly_crop_crs))
  crop_sif[crop_sif == 0] <- NA
  crop_sif_wgs84 <- terra::project(crop_sif, "EPSG:4326", method = "near")
  names(crop_sif_wgs84) <- paste0("crop_type_rank_", stringr::str_pad(sif_poly$sif_rank, 3, pad = "0"))
  raster::raster(crop_sif_wgs84)
}

top_sif_crop_rasters <- stats::setNames(lapply(seq_len(nrow(top_sif_polygons)), top_sif_crop_to_leaflet_raster), top_sif_polygons$layer_label)

top_sif_crop_map <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery")

for (group_name in names(top_sif_crop_rasters)) {
  top_sif_crop_map <- leaflet::addRasterImage(top_sif_crop_map, top_sif_crop_rasters[[group_name]], colors = crop_pal, opacity = 0.72, group = group_name, project = FALSE, maxBytes = 64 * 1024 * 1024)
}

top_sif_crop_map <- top_sif_crop_map %>%
  leaflet::addPolygons(data = top_sif_polygons, color = "#00d5ff", weight = 1.2, opacity = 1, fill = FALSE, group = "Top 100 SIF footprints", label = ~layer_label) %>%
  leaflet::addCircleMarkers(data = st_centroid(st_transform(top_sif_polygons, centroid_crs)) %>% st_transform(4326), radius = 3, color = "white", fillColor = "#ff2d55", fillOpacity = 0.9, weight = 1, group = "Top 100 SIF centroids", label = ~layer_label) %>%
  leaflet::addLegend(position = "bottomright", colors = unname(crop_palette[as.character(legend_classes$code)]), labels = gsub("_", " ", legend_classes$label), opacity = 0.85, title = "Crop type") %>%
  leaflet::addLayersControl(baseGroups = "Esri World Imagery", overlayGroups = c(names(top_sif_crop_rasters), "Top 100 SIF footprints", "Top 100 SIF centroids"), options = leaflet::layersControlOptions(collapsed = TRUE))

top_sif_bbox <- sf::st_bbox(top_sif_polygons)
top_sif_crop_map <- leaflet::fitBounds(top_sif_crop_map, lng1 = top_sif_bbox[["xmin"]], lat1 = top_sif_bbox[["ymin"]], lng2 = top_sif_bbox[["xmax"]], lat2 = top_sif_bbox[["ymax"]])

#htmlwidgets::saveWidget(top_sif_crop_map, file = top_sif_crop_leaflet_html, selfcontained = FALSE)

top_sif_crop_map

#-------------------------------------------------------------------------------
# daily OCO-2 scan-line approximation from centroid extremes

sif_daily_centroids <- sif_df %>%
  mutate(sif_date = as.Date(Delta_Date), measurement_mode = as.integer(Metadata.MeasurementMode), across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric)) %>%
  mutate(geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid() %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326)

sif_daily_centroid_df <- sif_daily_centroids %>%
  mutate(lon = st_coordinates(.)[, "X"], lat = st_coordinates(.)[, "Y"]) %>%
  st_drop_geometry()

scan_line_endpoints <- sif_daily_centroid_df %>%
  filter(measurement_mode %in% c(0L, 1L)) %>%
  group_by(sif_date, measurement_mode) %>%
  summarise(
    south_lon = lon[which.min(lat)],
    south_lat = lat[which.min(lat)],
    north_lon = lon[which.max(lat)],
    north_lat = lat[which.max(lat)],
    n_soundings = n(),
    .groups = "drop"
  ) %>%
  mutate(
    track_year = factor(lubridate::year(sif_date)),
    measurement_mode = factor(measurement_mode, levels = c(0, 1), labels = c("0 Nadir", "1 Glint")),
    track_label = paste0(sif_date, " | ", measurement_mode, " | n=", n_soundings),
    geometry = pmap(list(south_lon, south_lat, north_lon, north_lat), function(x1, y1, x2, y2) st_linestring(rbind(c(x1, y1), c(x2, y2))))
  ) %>%
  st_as_sf(crs = 4326)

scan_lines_nadir <- scan_line_endpoints %>% filter(measurement_mode == "0 Nadir")
scan_lines_glint <- scan_line_endpoints %>% filter(measurement_mode == "1 Glint")

track_year_pal <- leaflet::colorFactor(palette = grDevices::hcl.colors(length(unique(scan_line_endpoints$track_year)), "Dark 3"), domain = scan_line_endpoints$track_year)

oco2_daily_scan_lines_leaflet <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  leaflet::addPolylines(data = scan_lines_nadir, color = ~track_year_pal(track_year), weight = 4, opacity = 0.8, label = ~track_label, group = "Nadir scan lines") %>%
  leaflet::addPolylines(data = scan_lines_glint, color = ~track_year_pal(track_year), weight = 4, opacity = 0.8, label = ~track_label, group = "Glint scan lines") %>%
  leaflet::addLegend(position = "bottomright", pal = track_year_pal, values = scan_line_endpoints$track_year, title = "Year", opacity = 0.85) %>%
  leaflet::addLayersControl(baseGroups = "Esri World Imagery", overlayGroups = c("Nadir scan lines", "Glint scan lines"), options = leaflet::layersControlOptions(collapsed = FALSE))

scan_lines_bbox <- sf::st_bbox(scan_line_endpoints)
oco2_daily_scan_lines_leaflet <- leaflet::fitBounds(oco2_daily_scan_lines_leaflet, lng1 = scan_lines_bbox[["xmin"]], lat1 = scan_lines_bbox[["ymin"]], lng2 = scan_lines_bbox[["xmax"]], lat2 = scan_lines_bbox[["ymax"]])

#htmlwidgets::saveWidget(oco2_daily_scan_lines_leaflet, file = "eda_images/oco2_daily_scan_lines_leaflet.html", selfcontained = FALSE)

oco2_daily_scan_lines_leaflet

#-------------------------------------------------------------------------------
# sif polygons and gosif comparison, how close they match?

r <- rast('temp/GOSIF_2024.M04.tif')
r[r >= 32766] <- NA
r <- r * 0.0001
terra::plot(r)

global(r, fun = "min", na.rm = TRUE)
global(r, fun = "max", na.rm = TRUE)
nlyr(r)
names(r)
is.factor(r)
levels(r)
res(r)
crs(r)

sif_df <- readRDS('data/sif_sf_months2_7_cleaned.rds') 
sif_df <- sif_df %>% st_drop_geometry()
nrow(sif_df[format(as.Date(sif_df$Delta_Date), '%Y-%m') == '2024-04',])

sif_april_2024_nadir_polygons <- sif_df %>%
  mutate(sif_date = as.Date(Delta_Date), measurement_mode = as.integer(Metadata.MeasurementMode), across(all_of(c("Latitude", "Longitude", corner_cols, target_col)), as.numeric)) %>%
  filter(format(sif_date, "%Y-%m") == "2024-04", measurement_mode == 0L) %>%
  mutate(geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

sif_april_2024_nadir_plot <- ggplot() +
  geom_sf(data = germany_states, fill = "grey95", color = "grey45", linewidth = 0.25) +
  geom_sf(data = sif_april_2024_nadir_polygons, aes(fill = .data[[target_col]]), color = "red", linewidth = 0.08, alpha = 0.55) +
  scale_fill_viridis_c(option = "viridis", na.value = "transparent") +
  coord_sf(xlim = c(5.8, 15.2), ylim = c(47.2, 55.1), expand = FALSE) +
  labs(title = "April 2024 OCO-2 nadir SIF footprints over Germany", fill = "SIF") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.major = element_line(linewidth = 0.2, color = "grey82"))

print(sif_april_2024_nadir_plot)

gosif_april_crop <- terra::crop(r, terra::vect(sif_april_2024_nadir_polygons), snap = "out")

gosif_april_cells <- terra::as.polygons(gosif_april_crop, values = TRUE, na.rm = TRUE, dissolve = FALSE) %>%
  st_as_sf() %>%
  rename(gosif_value = all_of(names(gosif_april_crop)[1])) %>%
  mutate(gosif_cell = row_number()) %>%
  dplyr::select(gosif_cell, gosif_value, geometry)

gosif_cell_bounds <- do.call(rbind, lapply(st_geometry(gosif_april_cells), st_bbox)) %>%
  as_tibble() %>%
  rename(gosif_xmin = xmin, gosif_ymin = ymin, gosif_xmax = xmax, gosif_ymax = ymax) %>%
  mutate(gosif_centroid_lon = (gosif_xmin + gosif_xmax) / 2, gosif_centroid_lat = (gosif_ymin + gosif_ymax) / 2) %>%
  bind_cols(gosif_april_cells %>% st_drop_geometry() %>% dplyr::select(gosif_cell, gosif_value), .)

gosif_april_cells_3035 <- gosif_april_cells %>%
  st_transform(centroid_crs) %>%
  mutate(gosif_cell_area_m2 = as.numeric(st_area(geometry)))

sif_april_2024_nadir_for_overlap <- sif_april_2024_nadir_polygons %>%
  mutate(sif_poly_id = row_number(), oco2_sif = .data[[target_col]]) %>%
  dplyr::select(sif_poly_id, sif_date, oco2_sif, geometry) %>%
  st_transform(centroid_crs)

gosif_oco2_april_2024_intersections <- st_intersection(gosif_april_cells_3035, sif_april_2024_nadir_for_overlap) %>%
  mutate(overlap_area_m2 = as.numeric(st_area(geometry))) %>%
  st_drop_geometry()

gosif_oco2_april_2024_compare <- gosif_oco2_april_2024_intersections %>%
  group_by(gosif_cell, gosif_value, gosif_cell_area_m2) %>%
  summarise(
    n_oco2_polygons = n_distinct(sif_poly_id),
    oco2_overlap_area_m2 = sum(overlap_area_m2, na.rm = TRUE),
    oco2_overlap_fraction = oco2_overlap_area_m2 / first(gosif_cell_area_m2),
    oco2_sif_mean = mean(oco2_sif, na.rm = TRUE),
    oco2_sif_area_weighted_mean = weighted.mean(oco2_sif, overlap_area_m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gosif_cell_bounds, by = c("gosif_cell", "gosif_value")) %>%
  relocate(gosif_cell, gosif_xmin, gosif_ymin, gosif_xmax, gosif_ymax, gosif_centroid_lon, gosif_centroid_lat, gosif_value, gosif_cell_area_m2, oco2_overlap_area_m2, oco2_overlap_fraction, n_oco2_polygons, oco2_sif_mean, oco2_sif_area_weighted_mean) %>%
  arrange(desc(oco2_overlap_fraction), gosif_cell)

print(gosif_oco2_april_2024_compare)

gosif_oco2_april_2024_compare %>% 
  dplyr::select(gosif_cell, gosif_value, oco2_overlap_fraction, oco2_sif_mean, oco2_sif_area_weighted_mean) %>%
  mutate(diff = oco2_sif_area_weighted_mean - gosif_value) %>%
  summarise(mean_diff = mean(diff, na.rm = TRUE))

gosif_cell_920 <- gosif_april_cells %>% filter(gosif_cell == 920)

gosif_cell_920_overlap <- gosif_oco2_april_2024_intersections %>%
  filter(gosif_cell == 920) %>%
  dplyr::select(sif_poly_id, overlap_area_m2)

sif_cell_920_polygons <- sif_april_2024_nadir_polygons %>%
  mutate(sif_poly_id = row_number(), oco2_sif = .data[[target_col]]) %>%
  inner_join(gosif_cell_920_overlap, by = "sif_poly_id") %>%
  mutate(sif_label = paste0("poly ", sif_poly_id, " | OCO-2 SIF=", round(oco2_sif, 3), " | overlap m2=", round(overlap_area_m2, 1)))

gosif_cell_920_summary <- gosif_oco2_april_2024_compare %>%
  filter(gosif_cell == 920) %>%
  mutate(cell_label = paste0("GOSIF cell 920 | GOSIF=", round(gosif_value, 3), " | OCO-2 weighted mean=", round(oco2_sif_area_weighted_mean, 3), " | covered=", round(100 * oco2_overlap_fraction, 1), "%"))

gosif_cell_920 <- gosif_cell_920 %>%
  left_join(gosif_cell_920_summary %>% st_drop_geometry() %>% dplyr::select(gosif_cell, cell_label), by = "gosif_cell")

sif_cell_920_pal <- leaflet::colorNumeric(palette = "viridis", domain = sif_cell_920_polygons$oco2_sif, na.color = "transparent")

gosif_cell_920_map <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  leaflet::addPolygons(data = gosif_cell_920, color = "#ffd43b", weight = 4, opacity = 1, fill = FALSE, label = ~cell_label, group = "GOSIF cell 920") %>%
  leaflet::addPolygons(data = sif_cell_920_polygons, color = ~sif_cell_920_pal(oco2_sif), fillColor = ~sif_cell_920_pal(oco2_sif), weight = 1.2, opacity = 1, fillOpacity = 0.55, label = ~sif_label, group = "Overlapping OCO-2 SIF polygons") %>%
  leaflet::addLegend(position = "bottomright", pal = sif_cell_920_pal, values = sif_cell_920_polygons$oco2_sif, title = "OCO-2 SIF", opacity = 0.85) %>%
  leaflet::addLayersControl(baseGroups = "Esri World Imagery", overlayGroups = c("GOSIF cell 920", "Overlapping OCO-2 SIF polygons"), options = leaflet::layersControlOptions(collapsed = FALSE))

gosif_cell_920_bbox <- sf::st_bbox(st_transform(st_buffer(st_transform(gosif_cell_920, centroid_crs), 3000), 4326))
gosif_cell_920_map <- leaflet::fitBounds(gosif_cell_920_map, lng1 = gosif_cell_920_bbox[["xmin"]], lat1 = gosif_cell_920_bbox[["ymin"]], lng2 = gosif_cell_920_bbox[["xmax"]], lat2 = gosif_cell_920_bbox[["ymax"]])

gosif_cell_920_map


#-------------------------------------------------------------------------------

df <- readRDS('data/sif_sf_months2_7_cleaned.rds')
df <- df[lubridate::year(df$Delta_Date) %in% c('2019', '2020', '2021', '2022', '2023', '2024'),]
df <- df %>% sf::st_drop_geometry()

df %>%
  write_csv("338k_cnn_sif.csv")

#-------------------------------------------------------------------------------
df <- read.csv('data/extracted_modis_data/df_wpar_338k_crop_counts.csv')
colnames(df)

temp_df <- df[df$Quality_Flag == 0 & df$Metadata.MeasurementMode == 0,]

summary(temp_df[temp_df$Daily_SIF_740nm >= 0,]$active_growth_pct)
summary(temp_df[temp_df$Daily_SIF_740nm < 0,]$active_growth_pct)

sif_breaks <- seq(
  floor(min(temp_df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  ceiling(max(temp_df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  by = 0.25
)

temp_df_binned <- temp_df %>%
  mutate(
    sif_bin = cut(Daily_SIF_740nm, breaks = sif_breaks, include.lowest = TRUE, right = FALSE)
  )

ggplot(temp_df_binned, aes(x = Daily_SIF_740nm)) +
  geom_histogram(binwidth = 0.25, boundary = 0, fill = "steelblue", color = "white") +
  labs(
    title = "Daily SIF 740 nm Histogram",
    x = "Daily_SIF_740nm",
    y = "Count"
  ) +
  theme_minimal()


summary(df[df$Metadata.MeasurementMode == 0,]$Daily_SIF_740nm)
summary(df[df$Metadata.MeasurementMode == 1,]$Daily_SIF_740nm)



#-------------------------------------------
df <- read.csv('data/cnn_modis_chips/modis_8day_250m_16x16/chip_metadata.csv')
colnames(df)

summary(df$fapar_valid_fraction)
summary(df$evi_valid_fraction)
summary(df$ndvi_valid_fraction)
summary(df$par_valid_fraction)

sum(df$par_valid_fraction == 0)
sum(df$par_valid_fraction < 1)





#-------------------------------------------------------------------------------
# 338k sif data readjustments
library(tidyverse)
library(sf)
df <- read.csv('data/extracted_modis_data/df_wpar_338k_crop_counts.csv')
# keep only nadir and quality flag = 0
df <- tdf[tdf$Quality_Flag == 0 & tdf$Metadata.MeasurementMode == 0, ]

sif_breaks <- seq(
  floor(min(df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  ceiling(max(df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  by = 0.25
)

df_binned <- df %>%
  mutate(sif_bin = cut(Daily_SIF_740nm, breaks = sif_breaks, include.lowest = TRUE, right = FALSE))

table(df_binned$sif_bin)

# remove SIF outliers
#df <- df %>% filter(Daily_SIF_740nm >= -0.5, Daily_SIF_740nm < 1.75)

#IQR rule
q <- quantile(df$Daily_SIF_740nm, probs = c(0.25, 0.75), na.rm = TRUE)
iqr <- q[2] - q[1]

lower <- q[1] - 1.5 * iqr
upper <- q[2] + 1.5 * iqr

df <- df %>%
  filter(Daily_SIF_740nm >= lower, Daily_SIF_740nm <= upper)

# give each sif row a aggro climatic zone depending on which of these geometries they fall in
root_dir <- 'temp'
files_zones <- c('DE_HZ_5b.geojson', 'DE_HZ_6a.geojson', 'DE_HZ_6b.geojson', 'DE_HZ_7a.geojson', 
                 'DE_HZ_7b.geojson', 'DE_HZ_8a.geojson', 'DE_HZ_8b.geojson', 'DE_HZ_9a.geojson')

agro_climatic_zones <- files_zones %>%
  map_dfr(~ sf::read_sf(file.path(root_dir, .x), quiet = TRUE) %>% st_make_valid() %>% dplyr::select(hzs, geometry))

df_sif_polygons <- df %>%
  mutate(sif_row_id = row_number(), across(all_of(c("Latitude", "Longitude", corner_cols, target_col)), as.numeric)) %>%
  mutate(geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

df_sif_centroids <- df_sif_polygons %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326)

df_zone_lookup <- df_sif_centroids %>%
  dplyr::select(sif_row_id) %>%
  st_join(agro_climatic_zones %>% dplyr::select(hzs), join = st_within, left = TRUE) %>%
  st_drop_geometry() %>%
  dplyr::select(sif_row_id, hzs)

df <- df_sif_polygons %>%
  st_drop_geometry() %>%
  left_join(df_zone_lookup, by = "sif_row_id") %>%
  arrange(sif_row_id) %>%
  dplyr::select(-sif_row_id)

table(df$hzs, useNA = "ifany")

colSums(is.na(df))

#summary(df[is.na(df$active_growth_pct),]$crop_pixel_count)

df[is.na(df$ww_pct),]$ww_pct <- 0
df[is.na(df$active_growth_pct),]$active_growth_pct <- 0

final_df <- df %>% select(c(Daily_SIF_740nm, Delta_Time, Latitude, Longitude, Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
                            Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4, Delta_Date, state, hzs, active_growth_pct))
colSums(is.na(final_df))
final_df <- final_df %>% drop_na()
final_df %>%
  write_csv("data/148k_cnn_sif.csv")

# df <- df %>% drop_na()
# df %>%
#   write_csv("data/148k_nneighbour.csv")


tdf <- df %>% select(-c(sif_doy, closest_doy, closest_glass_date, glass_day_diff, mean_fapar, mean_par, apar))
tdf %>%
  write_csv("data/338k_base_crop_hzs.csv")
#-------------------------------------------------------------------------------
# take extreme  sif values and see their neighbourhood (how is the distribution)

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

library(tidyverse)
library(sf)
library(purrr)
#final_df <- read_csv('data/148k_nneighbour.csv')

evi_df <- read_csv('data/extracted_modis_data/338k_crop_hzs_evi.csv')
ndvi_df <- read_csv('data/extracted_modis_data/338k_crop_hzs_ndvi.csv')
wpar_df <- read_csv('data/extracted_modis_data/338k_crop_hzs_wpar.csv')

key_cols <- Reduce(intersect,list(names(evi_df), names(ndvi_df), names(wpar_df)))

final_df <- evi_df %>%
  full_join(ndvi_df, by = key_cols, suffix = c("_evi", "_ndvi")) %>%
  full_join(wpar_df, by = key_cols)

colnames(final_df)

#IQR rule

hist(final_df[final_df$Daily_SIF_740nm >= 1.5 & final_df$Daily_SIF_740nm < 1.75,]$Daily_SIF_740nm)

q <- quantile(final_df$Daily_SIF_740nm, probs = c(0.25, 0.75), na.rm = TRUE)
iqr <- q[2] - q[1]

lower <- q[1] - 1.5 * iqr
#upper <- q[2] + 1.5 * iqr
upper <- 1.75

final_df <- final_df %>%
  filter(Daily_SIF_740nm >= lower, Daily_SIF_740nm <= upper)


sif_breaks <- seq(
  floor(min(final_df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  ceiling(max(final_df$Daily_SIF_740nm, na.rm = TRUE) / 0.25) * 0.25,
  by = 0.25
)

df_binned <- final_df %>%
  mutate(sif_bin = cut(Daily_SIF_740nm, breaks = sif_breaks, include.lowest = TRUE, right = FALSE))

table(df_binned$sif_bin)

nrow(evi_df) - nrow(df_binned)
rm(evi_df, ndvi_df, wpar_df)
sort(colSums(is.na(final_df)))

# ggplot(df_binned, aes(x = Daily_SIF_740nm)) +
#   geom_histogram(binwidth = 0.25, boundary = 0, fill = "steelblue", color = "white") +
#   labs(title = "Daily SIF 740 nm Histogram",x = "Daily_SIF_740nm",y = "Count") +
#   theme_minimal()


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
      Daily_SIF_740nm >= 1 & Daily_SIF_740nm < 1.25 ~ "[1,1.25)",
      Daily_SIF_740nm >= 1.25 & Daily_SIF_740nm < 1.5 ~ "[1.25,1.5)",
      Daily_SIF_740nm >= 1.5 & Daily_SIF_740nm <= 1.75 ~ "[1.5,1.75]",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(extreme_bin)) %>%
  group_by(extreme_bin) %>%
  group_modify(~ slice_sample(.x, n = min(100, nrow(.x)))) %>%
  ungroup() %>%
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

extreme_sample_neighbor_summary <- neighbor_polygons %>%
  st_drop_geometry() %>%
  group_by(sample_id, sample_sif_row_id, extreme_bin) %>%
  summarise(
    n_neighbors = n(),
    neighbor_sif_min = min(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_mean = mean(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_median = median(Daily_SIF_740nm, na.rm = TRUE),
    neighbor_sif_max = max(Daily_SIF_740nm, na.rm = TRUE),
    max_d = max(neighbor_distance_m, na.rm = TRUE),
    neighbor_active_growth_pct_mean = mean(active_growth_pct, na.rm = TRUE),
    neighbor_crop_pixel_count_mean = mean(crop_pixel_count, na.rm = TRUE),
    neighbor_mean_fapar_mean = mean(mean_fapar, na.rm = TRUE),
    neighbor_mean_evi_mean = mean(mean_evi, na.rm = TRUE),
    neighbor_mean_ndvi_mean = mean(mean_ndvi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(extreme_samples %>% dplyr::select(sample_id, sample_sif_row_id = sif_row_id, sample_sif = Daily_SIF_740nm, sample_date = sif_date, sample_active_growth_pct = active_growth_pct, sample_crop_pixel_count = crop_pixel_count, sample_mean_fapar = mean_fapar, sample_mean_evi = mean_evi, sample_mean_ndvi = mean_ndvi), by = c("sample_id", "sample_sif_row_id")) %>%
  relocate(sample_id, sample_sif_row_id, sample_date, extreme_bin, sample_sif, sample_active_growth_pct, sample_crop_pixel_count, sample_mean_fapar, sample_mean_evi, sample_mean_ndvi, n_neighbors, max_d, neighbor_sif_min, neighbor_sif_mean, neighbor_sif_median, neighbor_sif_max, neighbor_active_growth_pct_mean, neighbor_crop_pixel_count_mean, neighbor_mean_fapar_mean, neighbor_mean_evi_mean, neighbor_mean_ndvi_mean) %>%
  arrange(extreme_bin, sample_id)

print(extreme_sample_neighbor_summary)

ex_sum <- extreme_sample_neighbor_summary %>% 
      select(c(sample_id, n_neighbors, max_d, sample_sif, neighbor_sif_mean,  
               sample_active_growth_pct,  neighbor_active_growth_pct_mean, 
               sample_crop_pixel_count, neighbor_crop_pixel_count_mean,
               sample_mean_fapar, neighbor_mean_fapar_mean,
               sample_mean_evi, neighbor_mean_evi_mean, 
               sample_mean_ndvi, neighbor_mean_ndvi_mean)) %>% 
      mutate(diff_sif = round(abs(neighbor_sif_mean - sample_sif),3),
             diff_active = round(abs(neighbor_active_growth_pct_mean - sample_active_growth_pct),3),
             diff_crop = round(abs(neighbor_crop_pixel_count_mean - sample_crop_pixel_count),3),
             diff_fapar = round(abs(neighbor_mean_fapar_mean - sample_mean_fapar),3),
             diff_evi = round(abs(neighbor_mean_evi_mean - sample_mean_evi),3),
             diff_ndvi = round(abs(neighbor_mean_ndvi_mean - sample_mean_ndvi),3))

ex_sum %>% select(sample_id,max_d,n_neighbors, sample_sif, neighbor_sif_mean, diff_sif, diff_active, diff_crop, diff_fapar, diff_evi, diff_ndvi) %>% print(n=500)
ex_sum %>% filter(sample_id == 59) %>% select(sample_id,max_d,n_neighbors, sample_sif, neighbor_sif_mean, diff_sif, diff_active, diff_crop, diff_fapar, diff_evi, diff_ndvi)

                                              
extreme_neighborhood_pal <- leaflet::colorNumeric(palette = "viridis", domain = c(extreme_sample_polygons$Daily_SIF_740nm, neighbor_polygons$Daily_SIF_740nm), na.color = "transparent")

extreme_neighborhood_groups <- c("Extreme samples [-0.75,-0.5)", "Neighbors [-0.75,-0.5)", "Extreme samples [-0.5,-0.25)", "Neighbors [-0.5,-0.25)", "Extreme samples [1,1.25)", "Neighbors [1,1.25)", "Extreme samples [1.25,1.5)", "Neighbors [1.25,1.5)", "Extreme samples [1.5,1.75]", "Neighbors [1.5,1.75]")

extreme_neighborhood_map <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  #leaflet::addProviderTiles(leaflet::providers$OpenStreetMap, group = "OpenStreetMap") %>%
  leaflet::addPolygons(data = extreme_sample_polygons %>% filter(extreme_bin == "[-0.75,-0.5)"), color = "black", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 4, opacity = 1, fillOpacity = 0.85, label = ~map_label, group = "Extreme samples [-0.75,-0.5)") %>%
  leaflet::addPolygons(data = neighbor_polygons %>% filter(extreme_bin == "[-0.75,-0.5)"), color = "#555555", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 0.7, opacity = 0.8, fillOpacity = 0.85, label = ~map_label, group = "Neighbors [-0.75,-0.5)") %>%
  leaflet::addPolygons(data = extreme_sample_polygons %>% filter(extreme_bin == "[-0.5,-0.25)"), color = "black", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 4, opacity = 1, fillOpacity = 0.85, label = ~map_label, group = "Extreme samples [-0.5,-0.25)") %>%
  leaflet::addPolygons(data = neighbor_polygons %>% filter(extreme_bin == "[-0.5,-0.25)"), color = "#555555", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 0.7, opacity = 0.8, fillOpacity = 0.85, label = ~map_label, group = "Neighbors [-0.5,-0.25)") %>%
  leaflet::addPolygons(data = extreme_sample_polygons %>% filter(extreme_bin == "[1,1.25)"), color = "black", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 4, opacity = 1, fillOpacity = 0.85, label = ~map_label, group = "Extreme samples [1,1.25)") %>%
  leaflet::addPolygons(data = neighbor_polygons %>% filter(extreme_bin == "[1,1.25)"), color = "#555555", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 0.7, opacity = 0.8, fillOpacity = 0.85, label = ~map_label, group = "Neighbors [1,1.25)") %>%
  leaflet::addPolygons(data = extreme_sample_polygons %>% filter(extreme_bin == "[1.25,1.5)"), color = "black", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 4, opacity = 1, fillOpacity = 0.85, label = ~map_label, group = "Extreme samples [1.25,1.5)") %>%
  leaflet::addPolygons(data = neighbor_polygons %>% filter(extreme_bin == "[1.25,1.5)"), color = "#555555", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 0.7, opacity = 0.8, fillOpacity = 0.85, label = ~map_label, group = "Neighbors [1.25,1.5)") %>%
  leaflet::addPolygons(data = extreme_sample_polygons %>% filter(extreme_bin == "[1.5,1.75]"), color = "black", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 4, opacity = 1, fillOpacity = 0.85, label = ~map_label, group = "Extreme samples [1.5,1.75]") %>%
  leaflet::addPolygons(data = neighbor_polygons %>% filter(extreme_bin == "[1.5,1.75]"), color = "#555555", fillColor = ~extreme_neighborhood_pal(Daily_SIF_740nm), weight = 0.7, opacity = 0.8, fillOpacity = 0.85, label = ~map_label, group = "Neighbors [1.5,1.75]") %>%
  leaflet::addLegend(position = "bottomright", pal = extreme_neighborhood_pal, values = c(extreme_sample_polygons$Daily_SIF_740nm, neighbor_polygons$Daily_SIF_740nm), title = "Daily SIF 740 nm", opacity = 0.85) %>%
  leaflet::addLayersControl(baseGroups = "Esri World Imagery", overlayGroups = extreme_neighborhood_groups, options = leaflet::layersControlOptions(collapsed = FALSE))

extreme_neighborhood_bbox <- sf::st_bbox(bind_rows(extreme_sample_polygons %>% dplyr::select(geometry), neighbor_polygons %>% dplyr::select(geometry)))
extreme_neighborhood_map <- leaflet::fitBounds(extreme_neighborhood_map, lng1 = extreme_neighborhood_bbox[["xmin"]], lat1 = extreme_neighborhood_bbox[["ymin"]], lng2 = extreme_neighborhood_bbox[["xmax"]], lat2 = extreme_neighborhood_bbox[["ymax"]])

extreme_neighborhood_map


