library(tidyverse)
library(sf)
library(terra)
library(patchwork)

# Compare the true 100 km MGRS geometries with the overlapping 109.8 km
# Sentinel-2 product extents. All geometries are plotted in the CRS of the
# Sentinel reference tifs.

sf::sf_use_s2(FALSE)

mgrs_rds <- "data/mgrs_de.rds"
mgrs_tif_dir <- "data/temp_data/mgrs_tifs"
output_dir <- "eda_images/mgrs_tile_geometry_comparison"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_mgrs_tile <- function(x) {
  tile <- stringr::str_trim(as.character(x))
  tile <- stringr::str_remove(tile, "^T")
  if_else(is.na(tile) | tile == "", NA_character_, paste0("T", tile))
}

extract_mgrs_tile_from_path <- function(path) {
  stringr::str_extract(basename(path), "T\\d{2}[A-Z]{3}")
}

mgrs_tif_paths <- list.files(
  mgrs_tif_dir,
  pattern = "\\.tif$",
  full.names = TRUE
)

if (length(mgrs_tif_paths) == 0) {
  stop("No Sentinel reference tifs found in: ", mgrs_tif_dir)
}

reference_raster <- terra::rast(mgrs_tif_paths[[1]])
sentinel_crs <- sf::st_crs(terra::crs(reference_raster))

sentinel_tile_extents <- map_dfr(mgrs_tif_paths, function(tif_path) {
  raster <- terra::rast(tif_path)

  if (sf::st_crs(terra::crs(raster)) != sentinel_crs) {
    stop("Sentinel reference tifs do not all use the same CRS: ", tif_path)
  }

  raster_extent <- unname(as.vector(terra::ext(raster)))
  tile_id <- extract_mgrs_tile_from_path(tif_path)

  if (length(raster_extent) != 4 || anyNA(raster_extent)) {
    stop("Invalid raster extent in: ", tif_path)
  }
  if (is.na(tile_id)) {
    stop("Could not extract an MGRS tile ID from: ", basename(tif_path))
  }

  # terra::ext() vector order is xmin, xmax, ymin, ymax. Build the polygon
  # explicitly rather than relying on bbox coercion between terra and sf.
  extent_polygon <- st_polygon(
    list(rbind(
      c(raster_extent[[1]], raster_extent[[3]]),
      c(raster_extent[[2]], raster_extent[[3]]),
      c(raster_extent[[2]], raster_extent[[4]]),
      c(raster_extent[[1]], raster_extent[[4]]),
      c(raster_extent[[1]], raster_extent[[3]])
    ))
  )

  st_sf(
    mgrs_tile_t = tile_id,
    source = "Sentinel-2 raster extent",
    tif_name = basename(tif_path),
    geometry = st_sfc(extent_polygon, crs = sentinel_crs)
  )
}) %>%
  filter(!is.na(mgrs_tile_t)) %>%
  arrange(mgrs_tile_t) %>%
  distinct(mgrs_tile_t, .keep_all = TRUE)

mgrs_raw <- readRDS(mgrs_rds)

if (!inherits(mgrs_raw, "sf")) {
  stop("mgrs_de.rds must contain an sf object.")
}
if (!"mgrs_tile" %in% names(mgrs_raw)) {
  stop("mgrs_de.rds does not contain a mgrs_tile column.")
}

selected_tiles <- sentinel_tile_extents$mgrs_tile_t

mgrs_tile_geometries <- mgrs_raw %>%
  mutate(mgrs_tile_t = normalize_mgrs_tile(mgrs_tile)) %>%
  filter(mgrs_tile_t %in% selected_tiles) %>%
  select(mgrs_tile_t) %>%
  st_make_valid() %>%
  group_by(mgrs_tile_t) %>%
  summarise(.groups = "drop") %>%
  st_transform(sentinel_crs) %>%
  mutate(source = "MGRS 100 km geometry") %>%
  arrange(mgrs_tile_t)

missing_mgrs_geometries <- setdiff(
  selected_tiles,
  mgrs_tile_geometries$mgrs_tile_t
)

if (length(missing_mgrs_geometries) > 0) {
  stop(
    "Missing geometries in mgrs_de.rds for: ",
    paste(missing_mgrs_geometries, collapse = ", ")
  )
}

comparison_sf <- bind_rows(
  mgrs_tile_geometries %>% select(mgrs_tile_t, source),
  sentinel_tile_extents %>% select(mgrs_tile_t, source)
) %>%
  mutate(
    source = factor(
      source,
      levels = c("MGRS 100 km geometry", "Sentinel-2 raster extent")
    )
  )

source_colors <- c(
  "MGRS 100 km geometry" = "#d73027",
  "Sentinel-2 raster extent" = "#2166ac"
)

source_fills <- scales::alpha(source_colors, 0.08)

tile_labels <- mgrs_tile_geometries %>%
  suppressWarnings(st_point_on_surface())

overview_plot <- ggplot() +
  geom_sf(
    data = comparison_sf,
    aes(color = source, fill = source),
    linewidth = 0.9,
    alpha = 0.5
  ) +
  geom_sf_text(
    data = tile_labels,
    aes(label = mgrs_tile_t),
    size = 3.2,
    fontface = "bold",
    color = "#222222"
  ) +
  scale_color_manual(values = source_colors, name = NULL) +
  scale_fill_manual(values = source_fills, name = NULL) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(
    title = "True MGRS geometries and Sentinel-2 raster extents",
    subtitle = "Red: non-overlapping 100 km MGRS geometry | Blue: overlapping Sentinel-2 product extent",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

tile_plot_list <- map(selected_tiles, function(tile_id) {
  tile_comparison <- comparison_sf %>%
    filter(mgrs_tile_t == tile_id)

  ggplot(tile_comparison) +
    geom_sf(
      aes(color = source, fill = source),
      linewidth = 1.1,
      alpha = 0.5
    ) +
    scale_color_manual(values = source_colors, name = NULL) +
    scale_fill_manual(values = source_fills, name = NULL) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = tile_id, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      axis.text = element_text(size = 6),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
})

per_tile_plot <- wrap_plots(tile_plot_list, ncol = 3) +
  plot_annotation(
    title = "Per-tile MGRS and Sentinel-2 extent comparison",
    subtitle = "Red: MGRS 100 km geometry | Blue: Sentinel-2 raster extent"
  )

print(overview_plot)
print(per_tile_plot)

ggsave(
  filename = file.path(output_dir, "mgrs_vs_sentinel_extents_overview.png"),
  plot = overview_plot,
  width = 11,
  height = 10,
  dpi = 220,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "mgrs_vs_sentinel_extents_per_tile.png"),
  plot = per_tile_plot,
  width = 13,
  height = 12,
  dpi = 220,
  bg = "white"
)

message("Plotted tiles: ", paste(selected_tiles, collapse = ", "))
message("Wrote plots to: ", output_dir)
