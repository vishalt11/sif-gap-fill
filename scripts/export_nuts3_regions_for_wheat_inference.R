library(sf)
library(dplyr)

target_mgrs_tiles <- c("32UNA", "32UQV", "32UPU", "32UPA", "32UPV")

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

# mgrs_bavaria_plot <- ggplot() +
#   geom_sf(data = germany_states_mgrs_bavaria, fill = "grey96", color = "grey65", linewidth = 0.25) +
#   geom_sf(data = bavaria_nuts1, fill = "#e8f3ff", color = "#1f78b4", linewidth = 0.7, alpha = 0.55) +
#   geom_sf(data = bavaria_nuts3, fill = NA, color = "#1f78b4", linewidth = 0.25) +
#   geom_sf(data = target_mgrs_bboxes, fill = NA, color = "#e31a1c", linewidth = 0.9) +
#   geom_sf_text(data = target_mgrs_labels, aes(label = mgrs_tile), color = "#e31a1c", size = 3.5, fontface = "bold") +
#   coord_sf(ylim = c(47, 51), xlim = c(8, 14), expand = FALSE) +
#   labs(title = "Selected MGRS Tile Borders With Bavaria NUTS 3", x = NULL, y = NULL) +
#   theme_minimal() +
#   theme(panel.grid.major = element_line(linewidth = 0.15, color = "grey85"), legend.position = "none")
# 
# print(mgrs_bavaria_plot)

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

nuts3_regions_50pct_in_mgrs <- bavaria_nuts3 %>%
  left_join(nuts3_tile_overlap %>% filter(nuts3_area_in_mgrs_pct >= 50), by = c("nuts_id", "nuts3")) %>%
  filter(!is.na(mgrs_tile)) %>%
  st_make_valid()

#--
preferred_existing_tiles <- c("32UNA", "32UQV", "32UPU")
new_tiles <- c("32UPA", "32UPV")

nuts3_regions_50pct_in_mgrs <- nuts3_regions_50pct_in_mgrs %>%
  mutate(
    tile_priority = case_when(
      mgrs_tile %in% preferred_existing_tiles ~ 2L,
      mgrs_tile %in% new_tiles ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  arrange(
    nuts_id,
    desc(tile_priority),
    desc(nuts3_area_in_mgrs_pct),
    mgrs_tile
  ) %>%
  group_by(nuts_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-tile_priority)

stopifnot(!anyDuplicated(nuts3_regions_50pct_in_mgrs$nuts_id))

saveRDS(nuts3_regions_50pct_in_mgrs, 'nuts3_regions_50pct_in_mgrs.rds')

mgrs_nuts3_50pct_plot <- ggplot() +
  geom_sf(data = bavaria_nuts1, fill = "grey96", color = "grey55", linewidth = 0.4) +
  geom_sf(data = bavaria_nuts3, fill = NA, color = "grey78", linewidth = 0.18) +
  geom_sf(data = nuts3_regions_50pct_in_mgrs, aes(fill = mgrs_tile), color = "#1f78b4", linewidth = 0.25, alpha = 0.65) +
  geom_sf(data = target_mgrs_bboxes, fill = NA, color = "#e31a1c", linewidth = 0.5) +
  geom_sf_text(data = target_mgrs_labels, aes(label = mgrs_tile), color = "#e31a1c", size = 3.5, fontface = "bold") +
  coord_sf(ylim = c(47, 51), xlim = c(8, 14), expand = FALSE) +
  labs(title = "Bavaria NUTS 3 Regions With At Least 50 Percent Area Inside Selected MGRS Tiles", fill = "MGRS tile", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_line(linewidth = 0.15, color = "grey85"), legend.position = "right")

print(mgrs_nuts3_50pct_plot)

ggsave("eda_images/pdf/mgrs_nuts3_50pct_plot.pdf", plot = mgrs_nuts3_50pct_plot, width = 10, height = 6, units = "in")


input_path <- "data/nuts3_regions_50pct_in_mgrs.rds"
output_path <- "data/nuts3_regions_50pct_in_mgrs.gpkg"
output_layer <- "nuts3_regions_50pct_in_mgrs"

nuts3_regions <- readRDS(input_path) %>%
  st_make_valid() %>%
  select(
    nuts_id,
    nuts3,
    mgrs_tile,
    nuts3_area_m2,
    overlap_area_m2,
    nuts3_area_in_mgrs_pct,
    geometry
  )

st_write(
  nuts3_regions,
  output_path,
  layer = output_layer,
  delete_dsn = TRUE,
  quiet = FALSE
)

message("Saved ", nrow(nuts3_regions), " NUTS3 regions to ", output_path)
