library(tidyverse)
library(sf)


library(giscoR)
library(terra)
library(leaflet)
library(htmlwidgets)


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

set.seed(42)
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

#-------------------------------------------------------------------------------
# Check df
df <- readRDS('data/base_sif/sif_sf_1_12_cleaned.rds')
df <- df %>% st_drop_geometry()
nrow(df)
colnames(df)

table(df$Metadata.MeasurementMode)
table(df$Quality_Flag)

unique(lubridate::year(df$Delta_Date))
df <- df[lubridate::year(df$Delta_Date) %in% c('2018', '2019', '2020', '2021', '2022', '2023', '2024'),]

colSums(is.na(df))

saveRDS(df, 'data/base_sif/sif_sf_1_12_cleaned.rds')

#-------------------------------------------------------------------------------
# combine modis data files and add hzs col
evi_df <- readRDS('data/extracted_modis_data/evi_1_12.rds')
evi_df <- evi_df %>% st_drop_geometry()
ndvi_df <- readRDS('data/extracted_modis_data/ndvi_1_12.rds')
ndvi_df <- ndvi_df %>% st_drop_geometry()
wpar_df <- readRDS('data/extracted_modis_data/wpar_1_12.rds')
wpar_df <- wpar_df %>% st_drop_geometry()

key_cols <- Reduce(intersect,list(names(evi_df), names(ndvi_df), names(wpar_df)))

final_df <- evi_df %>%
  full_join(ndvi_df, by = key_cols, suffix = c("_evi", "_ndvi")) %>%
  full_join(wpar_df, by = key_cols)

colSums(is.na(final_df))
colnames(final_df)

final_df <- final_df %>% select(-c(evi_doy, evi_start_date, evi_day_offset, 
                                   ndvi_doy, ndvi_start_date, ndvi_day_offset,
                                   fapar_doy, fapar_start_date, fapar_day_offset,
                                   par_doy, par_date, par_day_diff, sif_area_km2_ndvi))

rm(evi_df, ndvi_df, wpar_df)
saveRDS(final_df, 'data/extracted_modis_data/modis_1_12.rds')

final_df <- readRDS('data/extracted_modis_data/modis_1_12.rds')
# give each sif row a aggro climatic zone depending on which of these geometries they fall in
root_dir <- 'temp'
files_zones <- c('DE_HZ_5b.geojson', 'DE_HZ_6a.geojson', 'DE_HZ_6b.geojson', 'DE_HZ_7a.geojson', 
                 'DE_HZ_7b.geojson', 'DE_HZ_8a.geojson', 'DE_HZ_8b.geojson', 'DE_HZ_9a.geojson')

agro_climatic_zones <- files_zones %>%
  map_dfr(~ sf::read_sf(file.path(root_dir, .x), quiet = TRUE) %>% st_make_valid() %>% dplyr::select(hzs, geometry))

df_sif_polygons <- final_df %>%
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

final_df <- df_sif_polygons %>%
  st_drop_geometry() %>%
  left_join(df_zone_lookup, by = "sif_row_id") %>%
  arrange(sif_row_id) %>%
  dplyr::select(-sif_row_id)

final_df <- final_df %>%
  filter(!is.na(mean_par), !is.na(apar))

saveRDS(final_df, 'data/extracted_modis_data/modis_1_12.rds')

#-------------------------------------------------------------------------------

df <- read_csv('data/extracted_modis_data/modis_2_7_bin_uncertainity_corrected_M0QF0_gt0_6area_cnn.csv')

final_df <- df[df$final_check_modis_sif == 'accept' & df$Quality_Flag == 0,]

#final_df <- df[df$Metadata.MeasurementMode == 0 & df$final_check_modis_sif == 'accept' & df$Quality_Flag == 0,]
final_df <- final_df[lubridate::month(final_df$Delta_Date) %in% 2:7,]

inc <- 0.25
sif_breaks <- seq(
  floor(min(final_df$target_modis_sif, na.rm = TRUE) / inc) * inc,
  ceiling(max(final_df$target_modis_sif, na.rm = TRUE) / inc) * inc,
  by = inc
)
df_binned <- final_df %>%
  mutate(area_bin = cut(target_modis_sif, breaks = sif_breaks, include.lowest = TRUE, right = FALSE)) %>%
  select(area_bin)
table(df_binned$area_bin)


#-------------------------------------------------------------------------------
# Plot small areas and big areas
df <- readRDS('data/extracted_modis_data/modis_1_12_bin_uncertainity_corrected_raw.rds')
final_df <- df[df$Metadata.MeasurementMode == 0 & df$final_check_modis_sif == 'accept' & df$Quality_Flag == 0,]
#final_df <- df[df$final_check_modis_sif == 'accept' & df$Quality_Flag == 0,]
final_df <- final_df[lubridate::month(final_df$Delta_Date) %in% 2:7,]

#final_df <- final_df[final_df$sif_area_km2_evi >= 0.6,]

final_df <- final_df %>%
  mutate(crop_pixel_area_km2 = crop_pixel_count * 0.0001,
         crop_pixel_area_pct = (crop_pixel_area_km2 / sif_area_km2_evi) * 100)

final_df <- final_df[final_df$crop_pixel_area_pct < 0.1,]
final_df <- final_df[lubridate::month(final_df$Delta_Date) == 5,]

# final_df <- final_df %>% select(Daily_SIF_757nm, Daily_SIF_771nm, target_modis_sif,
#                                 Latitude, Longitude, Delta_Time,
#                                 Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
#                                 Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4,
#                                 sif_area_km2_evi, Delta_Date, Metadata.MeasurementMode, Quality_Flag, state,
#                                 hzs, final_check_757, final_check_771, final_check_modis_sif)
# 
# write_csv(final_df, 'data/extracted_modis_data/modis_2_7_bin_uncertainity_corrected_M0QF0_gt0_6area_cnn.csv')

low_crop_pct_polygons <- final_df %>%
  mutate(source_row = row_number(), crop_pixel_area_pct = as.numeric(crop_pixel_area_pct), target_modis_sif = as.numeric(target_modis_sif), across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric)) %>%
  filter(!is.na(crop_pixel_area_pct), !is.na(target_modis_sif), if_all(all_of(corner_cols), ~ !is.na(.x))) %>%
  slice_min(order_by = crop_pixel_area_pct, n = 1000, with_ties = FALSE) %>%
  arrange(crop_pixel_area_pct) %>%
  mutate(crop_pct_rank = row_number(), crop_pct_label = paste0("rank=", crop_pct_rank, " | crop pct=", round(crop_pixel_area_pct, 4), "% | target SIF=", round(target_modis_sif, 3), " | area=", round(sif_area_km2_evi, 4), " km2 | date=", Delta_Date), geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

low_crop_pct_pal <- leaflet::colorNumeric(palette = "viridis", domain = low_crop_pct_polygons$target_modis_sif, na.color = "transparent")
low_crop_pct_bbox <- st_bbox(low_crop_pct_polygons)

low_crop_pct_map <- leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  leaflet::addPolygons(data = low_crop_pct_polygons, color = ~low_crop_pct_pal(target_modis_sif), fillColor = ~low_crop_pct_pal(target_modis_sif), weight = 1.2, opacity = 1, fillOpacity = 0.90, label = ~crop_pct_label, group = "Lowest 100 crop pixel pct SIF polygons") %>%
  leaflet::addLegend(position = "bottomright", pal = low_crop_pct_pal, values = low_crop_pct_polygons$target_modis_sif, title = "target_modis_sif", opacity = 0.85) %>%
  leaflet::addLayersControl(baseGroups = "Esri World Imagery", overlayGroups = "Lowest 100 crop pixel pct SIF polygons", options = leaflet::layersControlOptions(collapsed = FALSE))

low_crop_pct_map <- leaflet::fitBounds(low_crop_pct_map, lng1 = low_crop_pct_bbox[["xmin"]], lat1 = low_crop_pct_bbox[["ymin"]], lng2 = low_crop_pct_bbox[["xmax"]], lat2 = low_crop_pct_bbox[["ymax"]])

low_crop_pct_map

#-------------------------------------------------------------------------------
# Plot sif polygons only in the respective mgrs tiles

df <- readRDS('data/extracted_modis_data/modis_1_12_bin_uncertainity_corrected_raw.rds')
#final_df <- df[df$Metadata.MeasurementMode == 0,]
final_df <- df[df$final_check_modis_sif == 'accept' & df$Quality_Flag == 0,]
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

# sif_mgrs_centroids_in_tiles %>%
#   st_drop_geometry() %>%
#   count(source_row, name = "n_tiles") %>%
#   filter(n_tiles > 1)

saveRDS(sif_mgrs_centroids_in_tiles, 'data/11tiles_2_7.rds')

#sif_mgrs_centroids_in_tiles <- sif_mgrs_centroids_in_tiles[sif_mgrs_centroids_in_tiles$Metadata.MeasurementMode == 0,]

#saveRDS(sif_mgrs_centroids_in_tiles, 'data/9tiles_M0_2_7.rds')

# germany_states_mgrs <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 1, resolution = "01", epsg = 4326) %>%
#   dplyr::select(state = NUTS_NAME, geometry) %>%
#   st_make_valid()
# 
# mgrs_sif_centroids_plot <- ggplot() +
#   geom_sf(data = germany_states_mgrs, fill = "grey95", color = "grey55", linewidth = 0.25) +
#   geom_sf(data = mgrs_tile_bboxes, fill = NA, color = "black", linewidth = 0.7) +
#   geom_sf_text(data = mgrs_tile_labels, aes(label = mgrs_tile), size = 3, color = "black") +
#   geom_sf(data = sif_mgrs_centroids_in_tiles, aes(color = target_modis_sif), size = 0.35, alpha = 0.7) +
#   scale_color_viridis_c(name = "target_modis_sif") +
#   labs(title = "SIF centroids inside selected MGRS tile bounding boxes", x = NULL, y = NULL) +
#   theme_minimal() +
#   theme(panel.grid.major = element_line(linewidth = 0.15, color = "grey85"), legend.position = "right")
# 
# print(mgrs_sif_centroids_plot)


#-------------------------------------------------------------------------------
# Get temporal ranges for each tile-year-month combo

geodes_root <- "data/geodes_wasp_zips"
sentinel_product_pattern <- "^SENTINEL2[A-Z]_\\d{8}-000000-000_L3A_T\\d{2}[A-Z]{3}_C_V[0-9]-[0-9]$"
dts_workers <- 2

dts_product_index <- list.dirs(geodes_root, recursive = TRUE, full.names = TRUE) %>%
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

write_csv(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges.csv")
saveRDS(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges.rds")

#-------------------------------------------------------------------------------

r <- rast('data/geodes_wasp_zips/32UNC/2020/SENTINEL2A_20200215-000000-000_L3A_T32UNC_C_V4-0/MASKS/SENTINEL2A_20200215-000000-000_L3A_T32UNC_C_V4-0_DTS_R1.tif')

#b4 <- rast('data/geodes_wasp_zips/32ULA/2023/SENTINEL2X_20230415-000000-000_L3A_T32ULA_C_V4-0/SENTINEL2X_20230415-000000-000_L3A_T32ULA_C_V4-0_FRC_B4.tif')
#b8 <- rast('data/geodes_wasp_zips/32ULA/2023/SENTINEL2X_20230415-000000-000_L3A_T32ULA_C_V4-0/SENTINEL2X_20230415-000000-000_L3A_T32ULA_C_V4-0_FRC_B8.tif')

#ndvi <- (b8-b4)/(b8+b4)

plot(ndvi)

r[r == -10000] <- 0
global(r, fun = function(x, ...) {
  c(
    min = min(x, na.rm = TRUE),
    q35 = quantile(x, 0.35, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
})

plot(r)

#-------------------------------------------------------------------------------
# convert days since start of year to date

dts_temporal_ranges <- readRDS("data/geodes_wasp_zips/dts_r1_temporal_ranges.rds")

dts_temporal_ranges <- dts_temporal_ranges %>%
  mutate(min_days = case_when(min != 0 ~ min, q10 != 0 ~ q10, median != 0 ~ median, TRUE ~ NA_real_),
         max_days = max,
         year_start = as.Date(paste0(product_year, "-01-01")),
         temporal_low = year_start + min_days,
         temporal_high = year_start + max_days) %>%
  dplyr::select(-year_start)

dts_temporal_ranges %>% select(c(mgrs_tile, product_date, min_days, max_days, temporal_low, temporal_high))

write_csv(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges_dates.csv")
saveRDS(dts_temporal_ranges, "data/geodes_wasp_zips/dts_r1_temporal_ranges_dates.rds")

#-------------------------------------------------------------------------------
# assign inrange/outrange to sif row based year-month match with product_date and 
# whether it falls within 

sif <- readRDS('data/9tiles_2_7.rds')
dts_temporal_ranges <- readRDS("data/geodes_wasp_zips/dts_r1_temporal_ranges_dates.rds")
dts_temporal_ranges[117,]$min_days <- 37
dts_temporal_ranges[117,]$temporal_low  <- as.Date('2020-02-07')

sif_date_aligned <- sif %>%
  mutate(sif_date = as.Date(Delta_Date), sif_year = lubridate::year(sif_date), sif_month = lubridate::month(sif_date)) %>%
  left_join(dts_temporal_ranges %>% dplyr::select(product_path, mgrs_tile, product_year, product_month, product_date, temporal_low, temporal_high, min_days, max_days), by = c("mgrs_tile" = "mgrs_tile", "sif_year" = "product_year", "sif_month" = "product_month")) %>%
  mutate(date_align = case_when(is.na(product_date) ~ "no_product", sif_date >= temporal_low & sif_date <= temporal_high ~ "inrange", TRUE ~ "outrange"))

sif_date_aligned <- sif_date_aligned %>% filter(date_align == 'inrange')

table(sif_date_aligned$Metadata.MeasurementMode)

saveRDS(sif_date_aligned, 'data/9tiles_2_7_inrange_M01.rds')

sif_date_aligned <- sif_date_aligned[sif_date_aligned$Metadata.MeasurementMode == 0,]

saveRDS(sif_date_aligned, 'data/9tiles_2_7_inrange_M0.rds')

sort(colSums(is.na(sif_date_aligned)))


#-------------------------------------------------------------------------------

df <- readRDS('data/9tiles_2_7_inrange_M01.rds')
table(df$Metadata.MeasurementMode)
df %>%
  st_drop_geometry() %>%
  distinct(Delta_Date, Metadata.MeasurementMode) %>%
  count(Delta_Date, name = "n_modes") %>%
  filter(n_modes > 1)

df %>%
  st_drop_geometry() %>%
  distinct(mgrs_tile, Delta_Date, Metadata.MeasurementMode) %>%
  count(mgrs_tile, Delta_Date, name = "n_modes") %>%
  filter(n_modes > 1)

summary(df[df$Metadata.MeasurementMode == 0,]$sif_area_km2_evi)
summary(df[df$Metadata.MeasurementMode == 1,]$sif_area_km2_evi)
summary(df[df$ww_pct >= 0.5 & lubridate::month(df$Delta_Date) == 2,]$target_modis_sif)

df %>%
  st_drop_geometry() %>%
  filter(ww_pct >= 0.5) %>%
  mutate(doy = lubridate::yday(Delta_Date),
         doy_bin = 8 * floor((doy - 1) / 8) + 1) %>%
  group_by(doy_bin) %>%
  summarise(n = n(),mean_sif = mean(target_modis_sif, na.rm = TRUE),median_sif = median(target_modis_sif, na.rm = TRUE),.groups = "drop") %>%
  print(n=100)

df <- df %>% select(c(target_modis_sif, Delta_Date))
df <- df %>% st_drop_geometry()

write_csv(df, 'data/sif_dates.csv')


#-------------------------------------------------------------------------------

df <- readRDS('data/9tiles_2_7_inrange_M01.rds')

colnames(df)

df <- df %>% select(Daily_SIF_757nm, Daily_SIF_771nm, target_modis_sif,
                    final_check_757, final_check_771, final_check_modis_sif,
                    Latitude, Longitude, Delta_Time, Delta_Date, sif_doy,
                    Lat_corner1, Lat_corner2, Lat_corner3, Lat_corner4,
                    Lon_corner1, Lon_corner2, Lon_corner3, Lon_corner4,
                    sif_area_km2_evi, Metadata.MeasurementMode, 
                    state, hzs, mgrs_tile, product_path, 
                    SZA, VZA, VAz, SAz)


df[is.na(df$hzs),]$hzs <- '8b'

deg_to_rad <- pi / 180
rad_to_deg <- 180 / pi

df <- df %>%
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

df <- df %>%
  mutate(
    signed_phase_angle = if_else(
      relative_azimuth < 0,
      -phase_angle,
      phase_angle
    )
  )

head(df,1)

#2022-07-29
#2022-07-31
#2024-07-13
#2024-07-15
#2024-07-29

#remove sif rows that dont have PAR data
dates_to_remove <- as.Date(c(
  "2022-07-29",
  "2022-07-31",
  "2024-07-13",
  "2024-07-15",
  "2024-07-29"
))

df <- df |>
  dplyr::filter(!Delta_Date %in% dates_to_remove)

saveRDS(df, 'data/9tiles_2_7_inrange_M01_cnn.rds')

df %>%
  st_drop_geometry() %>%
  write_csv('data/9tiles_2_7_inrange_M01_cnn.csv')

# sensor angle checks
# nadir_phase_check <- df %>%
#   filter(Metadata.MeasurementMode == 0) %>%
#   mutate(
#     relative_azimuth = ((SAz - VAz + 180) %% 360) - 180,
#     
#     phase_cos =
#       cos(SZA * pi / 180) * cos(VZA * pi / 180) +
#       sin(SZA * pi / 180) * sin(VZA * pi / 180) *
#       cos(relative_azimuth * pi / 180),
#     
#     phase_angle =
#       acos(pmax(-1, pmin(1, phase_cos))) * 180 / pi,
#     
#     phase_minus_sza = phase_angle - SZA,
#     absolute_difference = abs(phase_minus_sza)
#   )
# 
# summary(nadir_phase_check$phase_angle)
# summary(nadir_phase_check$phase_minus_sza)
# summary(nadir_phase_check$absolute_difference)
# 
# max(nadir_phase_check$absolute_difference, na.rm = TRUE)
# max(nadir_phase_check$VZA, na.rm = TRUE)

# Your values are also all well above the paper’s \(20^\circ\) region where the strongest hotspot/geometry effects occur. 
# Therefore, you do not have low-phase-angle hotspot observations in this nadir subset.

#-------------------------------------------------------------------------------
# plot 2022-05-07 phase angle and sif (compare angle shifts vs sif changes for same land composition)

df_1track <- df %>% filter(Delta_Date == as.Date('2022-05-07'))

df_1track_polygons <- df_1track %>%
  st_drop_geometry() %>%
  mutate(across(all_of(c("Latitude", "Longitude", corner_cols)), as.numeric),
         track_label = paste0("date=", Delta_Date, " | SIF=", round(target_modis_sif, 3), " | phase=", round(phase_angle, 2), " deg"),
         geometry = pmap(list(Lon_corner1, Lat_corner1, Lon_corner2, Lat_corner2, Lon_corner3, Lat_corner3, Lon_corner4, Lat_corner4), make_sif_polygon)) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

track_bbox <- st_bbox(df_1track_polygons)

track_sif_pal <- leaflet::colorNumeric(palette = "viridis", domain = df_1track_polygons$target_modis_sif, na.color = "transparent")

track_sif_leaflet <- leaflet::leaflet(data = df_1track_polygons, options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  leaflet::addPolygons(color = ~track_sif_pal(target_modis_sif), fillColor = ~track_sif_pal(target_modis_sif), fillOpacity = 0.75, opacity = 1, weight = 0.8, label = ~track_label, group = "SIF polygons") %>%
  leaflet::addLegend(position = "bottomright", pal = track_sif_pal, values = df_1track_polygons$target_modis_sif, title = "target_modis_sif", opacity = 0.85)

track_sif_leaflet <- leaflet::fitBounds(track_sif_leaflet, lng1 = track_bbox[["xmin"]], lat1 = track_bbox[["ymin"]], lng2 = track_bbox[["xmax"]], lat2 = track_bbox[["ymax"]])

track_phase_pal <- leaflet::colorNumeric(palette = "viridis", domain = df_1track_polygons$phase_angle, na.color = "transparent")

track_phase_leaflet <- leaflet::leaflet(data = df_1track_polygons, options = leaflet::leafletOptions(preferCanvas = TRUE)) %>%
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
  leaflet::addPolygons(color = ~track_phase_pal(phase_angle), fillColor = ~track_phase_pal(phase_angle), fillOpacity = 0.75, opacity = 1, weight = 0.8, label = ~track_label, group = "Phase angle polygons") %>%
  leaflet::addLegend(position = "bottomright", pal = track_phase_pal, values = df_1track_polygons$phase_angle, title = "phase_angle", opacity = 0.85)

track_phase_leaflet <- leaflet::fitBounds(track_phase_leaflet, lng1 = track_bbox[["xmin"]], lat1 = track_bbox[["ymin"]], lng2 = track_bbox[["xmax"]], lat2 = track_bbox[["ymax"]])

htmlwidgets::saveWidget(track_sif_leaflet, file = "eda_images/track_2022_05_07_target_modis_sif.html", selfcontained = FALSE)
htmlwidgets::saveWidget(track_phase_leaflet, file = "eda_images/track_2022_05_07_phase_angle.html", selfcontained = FALSE)

track_sif_leaflet
track_phase_leaflet

#-------------------------------------------------------------------------------

r <- rast('data/viirs_vnp18a2_daily_mean_par_germany_native/2023/VNP18A2.002_2023-04-10_Daily_Mean_PAR_VIIRS_Sinusoidal_native.tif')
r_check <- rast('data/viirs_vnp18a2_daily_mean_par_germany_native/2023/VNP18A2.002_2023-04-10_PAR_Quality_VIIRS_Sinusoidal_native.tif')

plot(r)

res(r)
crs(r)
res(r_check)
crs(r_check)

unique(r_check)

global(r_check, fun = function(x, ...) {
  c(
    min = min(x, na.rm = TRUE),
    q35 = quantile(x, 0.35, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
})

plot(r_check)
res(r)
crs(r)
res(r_check)
crs(r_check)

#-------------------------------------------------------------------------------
df <- readRDS("data/main_sif_data/2tiles_2_7_M01_QF01_inoutrange.rds")

library(sf)
library(tidyverse)
library(leaflet)
library(purrr)

df <- df %>%
  mutate(row_id = row_number(), Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

df_3035 <- df %>%
  st_transform(3035)

missing_hzs_points <- df %>%
  filter(is.na(hzs))

nearest_5_ids <- map_dfr(missing_hzs_points$row_id, function(id) {
  target_point <- df_3035 %>%
    filter(row_id == id)
  
  candidate_points <- df_3035 %>%
    filter(row_id != id)
  
  tibble(
    missing_row_id = id,
    row_id = candidate_points$row_id,
    distance_m = as.numeric(st_distance(candidate_points, target_point))
  ) %>%
    arrange(distance_m) %>%
    slice_head(n = 5)
})

nearest_5 <- df %>%
  inner_join(nearest_5_ids, by = "row_id")

nearest_5 %>%
  st_drop_geometry() %>%
  select(missing_row_id, row_id, distance_m, Daily_SIF_757nm, Daily_SIF_771nm,
         target_modis_sif, Delta_Date, Latitude, Longitude, state, hzs, mgrs_tile) %>%
  arrange(missing_row_id, distance_m)

missing_hzs_map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addCircleMarkers(data = nearest_5, radius = 5, color = "blue", fillColor = "blue",
                   fillOpacity = 0.8, group = "nearest_5",
                   label = ~paste0("nearest | missing_row_id=", missing_row_id,
                                   " | row_id=", row_id,
                                   " | d=", round(distance_m, 1), " m",
                                   " | hzs=", hzs)) %>%
  addCircleMarkers(data = missing_hzs_points, radius = 7, color = "red", fillColor = "red",
                   fillOpacity = 1, group = "missing_hzs",
                   label = ~paste0("missing hzs | row_id=", row_id,
                                   " | SIF=", round(target_modis_sif, 3))) %>%
  addLayersControl(overlayGroups = c("nearest_5", "missing_hzs"), options = layersControlOptions(collapsed = FALSE))

missing_hzs_map

#-------------------------------------------------------------------------------

df <- read_csv('data/sentinel2_spatial_aggregation_4000m_original_tiles_inout/fixed_grid_4000m_aggregate_manifest.csv')
sort(df$aggregated_target_modis_sif, decreasing = TRUE)

temp <- df %>% filter(aggregated_target_modis_sif > 0.9)
#-------------------------------------------------------------------------------

b2 <- rast('data/geodes_wasp_zips/32ULA/2019/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0_FRC_B2.tif')
b4 <- rast('data/geodes_wasp_zips/32ULA/2019/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0_FRC_B4.tif')
b8 <- rast('data/geodes_wasp_zips/32ULA/2019/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0/SENTINEL2B_20190215-000000-000_L3A_T32ULA_C_V4-0_FRC_B8.tif')

#r[r == -10000] <- 0
global(b8, fun = function(x, ...) {
  c(
    min = min(x, na.rm = TRUE),
    q35 = quantile(x, 0.35, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
})

plot(b8)


#-------------------------------------------------------------------------------
df <- read_csv('data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.csv')

inc <- 0.25
sif_breaks <- seq(
  floor(min(df$target_modis_sif, na.rm = TRUE) / inc) * inc,
  ceiling(max(df$target_modis_sif, na.rm = TRUE) / inc) * inc,
  by = inc
)
df_binned <- df %>%
  mutate(area_bin = cut(target_modis_sif, breaks = sif_breaks, include.lowest = TRUE, right = FALSE)) %>%
  select(area_bin)
table(df_binned$area_bin)

#-------------------------------------------------------------------------------

df1 <- read_csv(
  "data/pmeans_model_data/model_data_pmeans_6tiles.csv",
  show_col_types = FALSE
)

df2 <- read_csv(
  "data/main_sif_data/9tiles_2_7_M01_QF01_inoutrange_PARrm.csv",
  show_col_types = FALSE
)

coordinate_cols <- c(
  "Latitude", "Longitude",
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

# Six decimal places is approximately decimetre-level coordinate precision.
make_join_keys <- function(data) {
  data %>%
    mutate(
      join_time = as.numeric(as.POSIXct(Delta_Time, tz = "UTC")),
      across(
        all_of(coordinate_cols),
        ~ round(as.numeric(.x), 6),
        .names = "join_{.col}"
      )
    )
}

join_columns <- c(
  "join_time",
  paste0("join_", coordinate_cols)
)

df1_keyed <- make_join_keys(df1)
df2_keyed <- make_join_keys(df2)

# Only append columns that are not already present in df1.
columns_to_append <- setdiff(
  names(df2),
  c(names(df1), "Delta_Time", coordinate_cols)
)

# Make sure one geometry/time key does not match multiple df2 rows.
duplicate_df2_keys <- df2_keyed %>%
  group_by(across(all_of(join_columns))) %>%
  summarise(n_matches = n(), .groups = "drop") %>%
  filter(n_matches > 1)

if (nrow(duplicate_df2_keys) > 0) {
  stop(
    "df2 contains duplicated coordinate/time keys. ",
    "Inspect duplicate_df2_keys before joining."
  )
}

df2_for_join <- df2_keyed %>%
  select(
    all_of(join_columns),
    all_of(columns_to_append)
  )

joined_df <- df1_keyed %>%
  left_join(
    df2_for_join,
    by = join_columns,
    relationship = "many-to-one"
  ) %>%
  select(-all_of(join_columns))

# Confirm that the join did not duplicate df1 rows.
stopifnot(nrow(joined_df) == nrow(df1))

message("Original df1 rows: ", nrow(df1))
message(
  "Rows matched to target_modis_sif: ",
  sum(!is.na(joined_df$target_modis_sif))
)
message(
  "Unmatched rows: ",
  sum(is.na(joined_df$target_modis_sif))
)

colSums(is.na(joined_df))

joined_df <- joined_df %>%
  filter(!is.na(target_modis_sif))

joined_df[is.na(joined_df$ww_pct),]$ww_pct <- 0

joined_df <- joined_df %>%
  drop_na()

message("Complete rows retained: ", nrow(joined_df))

# Inspect the newly acquired target and metadata.
joined_df %>%
  select(
    Daily_SIF_740nm,
    Daily_SIF_757nm,
    Daily_SIF_771nm,
    target_modis_sif,
    final_check_modis_sif,
    hzs,
    date_align,
    product_path
  ) %>%
  head()

write_csv(joined_df,"data/pmeans_model_data/model_data_pmeans_6tiles_with_modis_sif.csv")

#-------------------------------------------------------------------------------

df <- read_csv('data/pmeans_model_data/model_data_pmeans_6tiles_with_modis_sif_fapar_par.csv')

df %>%
  filter(is.na(active_crop_pct)) %>%
  count(crop_pixel_count)

df <- df %>%
  mutate(
    active_crop_pct = if_else(
      is.na(active_crop_pct) & crop_pixel_count == 0,
      0,
      active_crop_pct
    )
  )

sort(colSums(is.na(df)))

df <- df %>% drop_na()

write_csv(df,"data/pmeans_model_data/model_data_pmeans_6tiles_with_modis_sif_fapar_par_cleaned.csv")

colnames(df)


#-------------------------------------------------------------------------------

target_mgrs_tiles <- c("32UNA", "32UQV", "32UPU")

target_mgrs_tif_paths <- list.files("data/temp_data/mgrs_tifs", pattern = "\\.tif$", full.names = TRUE) %>%
  keep(~ stringr::str_extract(basename(.x), "T[0-9]{2}[A-Z]{3}") %>% stringr::str_remove("^T") %in% target_mgrs_tiles)

target_mgrs_bboxes <- target_mgrs_tif_paths %>%
  map(function(tif_path) {
    tif_rast <- terra::rast(tif_path)
    tif_ext <- unname(as.vector(terra::ext(tif_rast)))
    tif_bbox <- st_bbox(c(xmin = tif_ext[1], ymin = tif_ext[3], xmax = tif_ext[2], ymax = tif_ext[4]), crs = st_crs(terra::crs(tif_rast)))
    st_sf(mgrs_tile = stringr::str_extract(basename(tif_path), "T[0-9]{2}[A-Z]{3}") %>% stringr::str_remove("^T"), tif_file = basename(tif_path), geometry = st_as_sfc(tif_bbox)) %>%
      st_transform(4326)
  }) %>%
  bind_rows() %>%
  st_make_valid()

target_mgrs_labels <- target_mgrs_bboxes %>%
  st_transform(centroid_crs) %>%
  st_centroid() %>%
  st_transform(4326)

germany_states_mgrs_bavaria <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 1, resolution = "01", epsg = 4326) %>%
  dplyr::select(nuts_id = NUTS_ID, state = NUTS_NAME, geometry) %>%
  st_make_valid()

bavaria_nuts1 <- germany_states_mgrs_bavaria %>%
  filter(nuts_id == "DE2")

bavaria_nuts3 <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 3, resolution = "01", epsg = 4326) %>%
  filter(stringr::str_starts(NUTS_ID, "DE2")) %>%
  dplyr::select(nuts_id = NUTS_ID, nuts3 = NUTS_NAME, geometry) %>%
  st_make_valid()

mgrs_bavaria_bbox <- st_bbox(st_union(st_geometry(bind_rows(bavaria_nuts1 %>% dplyr::select(geometry), target_mgrs_bboxes %>% dplyr::select(geometry)))))
mgrs_bavaria_xpad <- (mgrs_bavaria_bbox[["xmax"]] - mgrs_bavaria_bbox[["xmin"]]) * 0.08
mgrs_bavaria_ypad <- (mgrs_bavaria_bbox[["ymax"]] - mgrs_bavaria_bbox[["ymin"]]) * 0.08

mgrs_bavaria_plot <- ggplot() +
  geom_sf(data = germany_states_mgrs_bavaria, fill = "grey96", color = "grey65", linewidth = 0.25) +
  geom_sf(data = bavaria_nuts1, fill = "#e8f3ff", color = "#1f78b4", linewidth = 0.7, alpha = 0.55) +
  geom_sf(data = bavaria_nuts3, fill = NA, color = "#1f78b4", linewidth = 0.25) +
  geom_sf(data = target_mgrs_bboxes, fill = NA, color = "#e31a1c", linewidth = 0.9) +
  geom_sf_text(data = target_mgrs_labels, aes(label = mgrs_tile), color = "#e31a1c", size = 3.5, fontface = "bold") +
  coord_sf(ylim = c(47, 51), xlim = c(8, 14), expand = FALSE) +
  labs(title = "Selected MGRS Tile Borders With Bavaria NUTS 3", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.15, color = "grey85"), legend.position = "none")

print(mgrs_bavaria_plot)

nuts3_tile_overlap <- st_intersection(
  bavaria_nuts3 %>%
    st_transform(centroid_crs) %>%
    mutate(nuts3_area_m2 = as.numeric(st_area(geometry))),
  target_mgrs_bboxes %>%
    st_transform(centroid_crs) %>%
    dplyr::select(mgrs_tile)
) %>%
  mutate(overlap_area_m2 = as.numeric(st_area(geometry)),
         nuts3_area_in_mgrs_pct = (overlap_area_m2 / nuts3_area_m2) * 100) %>%
  st_drop_geometry() %>%
  dplyr::select(nuts_id, nuts3, mgrs_tile, nuts3_area_m2, overlap_area_m2, nuts3_area_in_mgrs_pct)

nuts3_regions_80pct_in_mgrs <- bavaria_nuts3 %>%
  left_join(nuts3_tile_overlap %>% filter(nuts3_area_in_mgrs_pct >= 80), by = c("nuts_id", "nuts3")) %>%
  filter(!is.na(mgrs_tile)) %>%
  st_make_valid()


saveRDS(nuts3_regions_80pct_in_mgrs, 'nuts3_regions_80pct_in_mgrs.rds')

mgrs_nuts3_80pct_plot <- ggplot() +
  geom_sf(data = bavaria_nuts1, fill = "grey96", color = "grey55", linewidth = 0.4) +
  geom_sf(data = bavaria_nuts3, fill = NA, color = "grey78", linewidth = 0.18) +
  geom_sf(data = nuts3_regions_80pct_in_mgrs, aes(fill = mgrs_tile), color = "#1f78b4", linewidth = 0.25, alpha = 0.65) +
  geom_sf(data = target_mgrs_bboxes, fill = NA, color = "#e31a1c", linewidth = 0.9) +
  geom_sf_text(data = target_mgrs_labels, aes(label = mgrs_tile), color = "#e31a1c", size = 3.5, fontface = "bold") +
  coord_sf(ylim = c(47, 51), xlim = c(8, 14), expand = FALSE) +
  labs(title = "Bavaria NUTS 3 Regions With At Least 80 Percent Area Inside Selected MGRS Tiles", fill = "MGRS tile", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.15, color = "grey85"), legend.position = "right")

print(mgrs_nuts3_80pct_plot)

#-------------------------------------------------------------------------------

availability_mgrs_tiles <- c("32UPU", "32UNA", "32UQV")
availability_product_pattern <- "^SENTINEL2[A-Z]_\\d{8}-000000-000_L3A_T\\d{2}[A-Z]{3}_C_V[0-9]-[0-9]$"

mgrs_product_folder_index <- availability_mgrs_tiles %>%
  map(~ list.dirs(file.path("data/geodes_wasp_zips", .x), recursive = TRUE, full.names = TRUE)) %>%
  unlist(use.names = FALSE) %>%
  keep(~ stringr::str_detect(basename(.x), availability_product_pattern)) %>%
  tibble(product_path = .) %>%
  mutate(product_id = basename(product_path),
         mgrs_tile = stringr::str_remove(stringr::str_extract(product_id, "T\\d{2}[A-Z]{3}"), "^T"),
         product_date = as.Date(stringr::str_extract(product_id, "\\d{8}"), format = "%Y%m%d"),
         product_year = lubridate::year(product_date),
         product_month = lubridate::month(product_date),
         year_month = format(product_date, "%Y-%m"),
         product_version = stringr::str_extract(product_id, "V[0-9]-[0-9]$")) %>%
  arrange(mgrs_tile, product_year, product_month)

mgrs_year_month_availability <- mgrs_product_folder_index %>%
  group_by(mgrs_tile, product_year, product_month, year_month) %>%
  summarise(n_products = n(), product_ids = paste(product_id, collapse = "; "), product_paths = paste(product_path, collapse = "; "), .groups = "drop") %>%
  mutate(is_available = TRUE) %>%
  complete(mgrs_tile = availability_mgrs_tiles, product_year = min(product_year):max(product_year), product_month = 1:12, fill = list(n_products = 0, product_ids = NA_character_, product_paths = NA_character_, is_available = FALSE)) %>%
  mutate(year_month = if_else(is.na(year_month), sprintf("%s-%02d", product_year, product_month), year_month),
         month_label = factor(month.abb[product_month], levels = month.abb)) %>%
  arrange(mgrs_tile, product_year, product_month)

mgrs_year_month_availability %>%
  filter(is_available) %>%
  dplyr::select(mgrs_tile, product_year, product_month, year_month, n_products, product_ids, product_paths)

mgrs_year_month_availability_plot <- ggplot(mgrs_year_month_availability, aes(x = month_label, y = factor(product_year), fill = is_available)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(is_available, "yes", "")), size = 2.6, color = "white") +
  facet_wrap(~ mgrs_tile, ncol = 1) +
  scale_fill_manual(values = c("TRUE" = "#1f78b4", "FALSE" = "grey88"), labels = c("FALSE" = "missing", "TRUE" = "available"), name = NULL) +
  labs(title = "Available Sentinel-2 Monthly Composite Folders", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))

print(mgrs_year_month_availability_plot)


#-------------------------------------------------------------------------------

df <- read_csv('data/pmeans_model_data/model_data_pmeans_6tiles_with_modis_sif_fapar_par_cleaned.csv')
colnames(df)

colSums(is.na(df))

#-------------------------------------------------------------------------------

germany_state_borders <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 1, resolution = "01", epsg = 4326) %>%
  dplyr::select(state = NUTS_NAME, nuts_id = NUTS_ID, geometry) %>%
  st_make_valid()

germany_boundaries_gpkg <- "data/germany_boundaries.gpkg"

st_write(germany_state_borders, germany_boundaries_gpkg, delete_dsn = TRUE)

#-------------------------------------------------------------------------------

hzs_geojson_paths <- list.files("temp/hzs", pattern = "\\.geojson$", full.names = TRUE)

hzs_zones <- hzs_geojson_paths %>%
  map(sf::read_sf) %>%
  bind_rows() %>%
  st_make_valid()

st_write(hzs_zones, "data/hzs_zones.gpkg", delete_dsn = TRUE)

#-------------------------------------------------------------------------------



