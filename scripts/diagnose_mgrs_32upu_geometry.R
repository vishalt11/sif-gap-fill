library(tidyverse)
library(sf)
library(giscoR)
library(patchwork)

# Diagnose whether rows originally assigned to Sentinel/MGRS tile 32UPU are
# spatially consistent with the 32UPU geometry stored in mgrs_de.rds.

sf::sf_use_s2(FALSE)

input_csv <- "data/9tiles_2_7_inrange_M01_cnn.csv"
mgrs_rds <- "data/mgrs_de.rds"
target_tile <- "32UPU"
utm_crs <- 32632
output_dir <- "eda_images/mgrs_32upu_geometry_diagnostic"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

normalize_mgrs_tile <- function(x) {
  stringr::str_remove(
    stringr::str_to_upper(stringr::str_trim(as.character(x))),
    "^T"
  )
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

df <- read_csv(
  input_csv,
  col_types = cols(
    Delta_Date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
  mutate(
    source_csv_row = row_number(),
    mgrs_tile_normalized = normalize_mgrs_tile(mgrs_tile)
  )

required_cols <- c("mgrs_tile", corner_cols)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_32upu <- df %>%
  filter(mgrs_tile_normalized == target_tile) %>%
  mutate(across(all_of(corner_cols), as.numeric)) %>%
  filter(if_all(all_of(corner_cols), ~ is.finite(.x)))

if (nrow(df_32upu) == 0) {
  stop("No rows have mgrs_tile == ", target_tile, ".")
}

sif_polygons_32upu <- df_32upu %>%
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
  st_make_valid() %>%
  st_transform(utm_crs)

sif_centroids_32upu_utm <- suppressWarnings(st_centroid(sif_polygons_32upu))

mgrs <- readRDS(mgrs_rds)
if (!inherits(mgrs, "sf")) {
  stop("mgrs_de.rds must contain an sf object.")
}
if (!"mgrs_tile" %in% names(mgrs)) {
  stop("mgrs_de.rds does not contain a mgrs_tile column.")
}

mgrs_32upu_utm <- mgrs %>%
  mutate(mgrs_tile_normalized = normalize_mgrs_tile(mgrs_tile)) %>%
  filter(mgrs_tile_normalized == target_tile) %>%
  select(mgrs_tile_normalized) %>%
  st_make_valid() %>%
  group_by(mgrs_tile_normalized) %>%
  summarise(.groups = "drop") %>%
  st_transform(utm_crs)

if (nrow(mgrs_32upu_utm) == 0) {
  stop("No ", target_tile, " geometry exists in mgrs_de.rds.")
}

centroid_inside_tile <- lengths(
  st_within(sif_centroids_32upu_utm, mgrs_32upu_utm)
) > 0

sif_centroids_32upu_utm <- sif_centroids_32upu_utm %>%
  mutate(
    centroid_inside_mgrs_geometry = centroid_inside_tile,
    point_class = if_else(
      centroid_inside_mgrs_geometry,
      "Inside mgrs_de 32UPU",
      "Outside mgrs_de 32UPU"
    )
  )

membership_summary <- sif_centroids_32upu_utm %>%
  st_drop_geometry() %>%
  count(point_class, name = "n_centroids") %>%
  mutate(
    pct_centroids = n_centroids / sum(n_centroids)
  )

message("Original ", target_tile, " rows plotted: ", nrow(df_32upu))
print(membership_summary, n = Inf)

centroid_xy <- st_coordinates(sif_centroids_32upu_utm)
centroid_membership <- sif_centroids_32upu_utm %>%
  st_drop_geometry() %>%
  transmute(
    source_csv_row,
    mgrs_tile,
    mgrs_tile_normalized,
    Delta_Date,
    product_path,
    centroid_x_utm = centroid_xy[, "X"],
    centroid_y_utm = centroid_xy[, "Y"],
    centroid_inside_mgrs_geometry,
    point_class
  )

write_csv(
  centroid_membership,
  file.path(output_dir, "mgrs_32upu_centroid_membership.csv")
)

germany <- gisco_get_countries(
  country = "Germany",
  resolution = "10",
  epsg = 4326
) %>%
  st_make_valid()

sif_centroids_32upu <- st_transform(sif_centroids_32upu_utm, 4326)
mgrs_32upu <- st_transform(mgrs_32upu_utm, 4326)

point_colors <- c(
  "Inside mgrs_de 32UPU" = "#2166ac",
  "Outside mgrs_de 32UPU" = "#f28e2b"
)

germany_plot <- ggplot() +
  geom_sf(
    data = germany,
    fill = "grey96",
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_sf(
    data = mgrs_32upu,
    fill = scales::alpha("#d73027", 0.18),
    color = "#d73027",
    linewidth = 1.1
  ) +
  geom_sf(
    data = sif_centroids_32upu,
    aes(color = point_class),
    size = 0.6,
    alpha = 0.55
  ) +
  scale_color_manual(values = point_colors, name = NULL) +
  coord_sf(crs = st_crs(4326), datum = NA, expand = FALSE) +
  labs(
    title = "Original 32UPU SIF centroids across Germany",
    subtitle = "Red polygon: 32UPU geometry from mgrs_de.rds",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

point_bbox <- st_bbox(sif_centroids_32upu)
tile_bbox <- st_bbox(mgrs_32upu)
zoom_x <- c(
  min(point_bbox[["xmin"]], tile_bbox[["xmin"]]) - 0.12,
  max(point_bbox[["xmax"]], tile_bbox[["xmax"]]) + 0.12
)
zoom_y <- c(
  min(point_bbox[["ymin"]], tile_bbox[["ymin"]]) - 0.08,
  max(point_bbox[["ymax"]], tile_bbox[["ymax"]]) + 0.08
)

zoom_plot <- ggplot() +
  geom_sf(
    data = germany,
    fill = "grey96",
    color = "grey70",
    linewidth = 0.25
  ) +
  geom_sf(
    data = mgrs_32upu,
    fill = scales::alpha("#d73027", 0.18),
    color = "#d73027",
    linewidth = 1.2
  ) +
  geom_sf(
    data = sif_centroids_32upu,
    aes(color = point_class),
    size = 1,
    alpha = 0.65
  ) +
  scale_color_manual(values = point_colors, name = NULL) +
  coord_sf(
    xlim = zoom_x,
    ylim = zoom_y,
    crs = st_crs(4326),
    datum = NA,
    expand = FALSE
  ) +
  labs(
    title = "32UPU geometry diagnostic",
    subtitle = paste0(
      sum(centroid_inside_tile), " of ", length(centroid_inside_tile),
      " centroids fall inside the mgrs_de.rds geometry"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

combined_plot <- germany_plot + zoom_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(germany_plot)
print(zoom_plot)
print(combined_plot)

ggsave(
  file.path(output_dir, "mgrs_32upu_germany.png"),
  germany_plot,
  width = 8,
  height = 9,
  dpi = 220,
  bg = "white"
)

ggsave(
  file.path(output_dir, "mgrs_32upu_zoom.png"),
  zoom_plot,
  width = 9,
  height = 8,
  dpi = 220,
  bg = "white"
)

ggsave(
  file.path(output_dir, "mgrs_32upu_germany_and_zoom.png"),
  combined_plot,
  width = 15,
  height = 8,
  dpi = 220,
  bg = "white"
)

message("Wrote diagnostic outputs to: ", output_dir)
