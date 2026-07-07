library(tidyverse)
library(sf)
library(terra)

crop_type_dir <- "data/crop_type_tif"

#"NIEDERSACHSEN"                          "SACHSEN-ANHALT"                                         "BAYERN" 
#c("32UND", "32UNC", "32UPC")             c("32UPC", "32UPB", "32UQC", "32UQB")       c("32UNA", "32UPA", "32UPV", "32UQV", "32UPU", "32UQU")
#"data/ns_sif_mgrs_crop_composition.rds"  "data/sa_sif_mgrs_crop_composition.rds"     "data/ba_sif_mgrs_crop_composition.rds"




NUTS_NAME <- "SCHLESWIG-HOLSTEIN" 
target_mgrs_tiles <- c("32UNE")
output_file <- "data/sh_sif_mgrs_crop_composition.rds"



ns_sf <- readRDS("data/sif_sf_months2_7_cleaned.rds")
mgrs_de <- readRDS("data/mgrs_de.rds")

crop_classes <- readr::read_delim(
  file.path(crop_type_dir, "LEGEND_CropTypes.txt"),
  delim = "\t",
  show_col_types = FALSE
)
colnames(crop_classes) <- c("code", "label")
crop_classes <- crop_classes %>%
  mutate(
    code = as.integer(code),
    count_col = paste0(
      "crop_count_",
      label %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )

mgrs_target <- mgrs_de %>%
  filter(mgrs_tile %in% target_mgrs_tiles)

if (nrow(mgrs_target) != length(target_mgrs_tiles)) {
  missing_tiles <- setdiff(target_mgrs_tiles, mgrs_target$mgrs_tile)
  stop("Missing MGRS tiles in mgrs_de: ", paste(missing_tiles, collapse = ", "))
}

ns_sf <- ns_sf %>%
  filter(state == NUTS_NAME) %>%
  mutate(
    sif_id = row_number(),
    sif_year = as.integer(format(as.Date(Delta_Date), "%Y"))
  ) %>%
  filter(sif_year <= 2024)

ns_centers <- ns_sf %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

ns_mgrs_points <- st_join(
  ns_centers,
  mgrs_target %>% select(mgrs_tile),
  join = st_within,
  left = FALSE
)

message("Retained ", nrow(ns_mgrs_points), " SA SIF points in target MGRS tiles.")

make_sif_polygon <- function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
  st_polygon(list(rbind(
    c(lon1, lat1),
    c(lon2, lat2),
    c(lon3, lat3),
    c(lon4, lat4),
    c(lon1, lat1)
  )))
}

sif_attrs <- ns_mgrs_points %>%
  st_drop_geometry()

sif_polygon_geometry <- pmap(
  list(
    sif_attrs$Lon_corner1, sif_attrs$Lat_corner1,
    sif_attrs$Lon_corner2, sif_attrs$Lat_corner2,
    sif_attrs$Lon_corner3, sif_attrs$Lat_corner3,
    sif_attrs$Lon_corner4, sif_attrs$Lat_corner4
  ),
  make_sif_polygon
)

sif_polygons <- st_sf(
  sif_attrs,
  geometry = do.call(st_sfc, c(sif_polygon_geometry, list(crs = 4326)))
) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

count_crop_pixels_for_year <- function(year, polygons_year, tiles_sf, crop_classes) {
  message("Processing crop type raster for ", year, " with ", nrow(polygons_year), " SIF polygons.")

  raster_file <- file.path(crop_type_dir, paste0("croptypes_", year, ".tif"))

  if (!file.exists(raster_file)) {
    warning("No crop type raster found for ", year, ": ", raster_file)
    return(tibble(sif_id = polygons_year$sif_id))
  }

  crop_raster <- rast(raster_file)
  levels(crop_raster) <- data.frame(value = crop_classes$code, crop = crop_classes$label)

  polygons_utm <- polygons_year %>%
    st_transform(crs(crop_raster)) %>%
    mutate(extract_id = row_number())

  tiles_union <- st_sf(
    geometry = st_sfc(st_union(st_geometry(tiles_sf)), crs = st_crs(tiles_sf))
  ) %>%
    st_transform(crs(crop_raster))

  crop_raster_tiles <- crop(crop_raster, vect(tiles_union))

  extracted <- terra::extract(crop_raster_tiles, vect(polygons_utm))
  value_col <- setdiff(names(extracted), "ID")[1]

  id_lookup <- polygons_utm %>%
    st_drop_geometry() %>%
    transmute(ID = extract_id, sif_id)

  counts_long <- extracted %>%
    rename(crop_value = all_of(value_col)) %>%
    mutate(
      crop_value = as.character(crop_value),
      crop_code = suppressWarnings(as.integer(crop_value)),
      crop_code = if_else(
        is.na(crop_code),
        crop_classes$code[match(crop_value, crop_classes$label)],
        crop_code
      )
    ) %>%
    filter(!is.na(crop_code)) %>%
    left_join(id_lookup, by = "ID") %>%
    count(sif_id, crop_code, name = "n")

  expand_grid(
    sif_id = polygons_year$sif_id,
    crop_code = crop_classes$code
  ) %>%
    left_join(counts_long, by = c("sif_id", "crop_code")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    left_join(crop_classes %>% select(code, count_col), by = c("crop_code" = "code")) %>%
    select(sif_id, count_col, n) %>%
    pivot_wider(
      names_from = count_col,
      values_from = n,
      values_fill = list(n = 0L),
      values_fn = list(n = sum)
    ) %>%
    mutate(crop_pixel_count = rowSums(across(starts_with("crop_count_")))) %>%
    relocate(crop_pixel_count, .after = sif_id)
}

crop_counts <- sif_polygons %>%
  st_drop_geometry() %>%
  distinct(sif_year) %>%
  arrange(sif_year) %>%
  pull(sif_year) %>%
  map_dfr(function(year) {
    polygons_year <- sif_polygons %>% filter(sif_year == year)
    count_crop_pixels_for_year(year, polygons_year, mgrs_target, crop_classes)
  })

ns_sif_mgrs_crop_composition <- ns_mgrs_points %>%
  left_join(crop_counts, by = "sif_id")

ns_sif_mgrs_crop_composition <- ns_sif_mgrs_crop_composition %>%
  mutate(
    ww_pct = if_else(
      crop_pixel_count > 0,
      crop_count_winter_wheat / crop_pixel_count,
      NA_real_
    )
  )

saveRDS(ns_sif_mgrs_crop_composition, output_file)

message("Saved crop composition output to ", output_file)

#-------------------------------------
library(tidyverse)

df <- readRDS('data/sh_sif_mgrs_crop_composition.rds')
df <- df %>% sf::st_drop_geometry()
write.csv(df, "data/sh_sif_mgrs_crop_composition.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# Visual verification of sif polygon's crop composition

library(leaflet)

z <- ns_sif_mgrs_crop_composition %>%
  filter(sif_id == 87)

crop_classes <- readr::read_delim(
  "data/crop_type_tif/LEGEND_CropTypes.txt",
  delim = "\t",
  show_col_types = FALSE
)
colnames(crop_classes) <- c("code", "label")
crop_classes <- crop_classes %>%
  mutate(code = as.integer(code))

ctr <- rast("data/crop_type_tif/croptypes_2018.tif")
levels(ctr) <- data.frame(value = crop_classes$code, crop = crop_classes$label)

tile_818 <- mgrs_de %>%
  filter(mgrs_tile == z$mgrs_tile)

z_poly <- st_polygon(list(rbind(
  c(z$Lon_corner1, z$Lat_corner1),
  c(z$Lon_corner2, z$Lat_corner2),
  c(z$Lon_corner3, z$Lat_corner3),
  c(z$Lon_corner4, z$Lat_corner4),
  c(z$Lon_corner1, z$Lat_corner1)
))) %>%
  st_sfc(crs = 4326) %>%
  st_sf(sif_id = z$sif_id, geometry = .) %>%
  st_make_valid()

z_center <- z %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# This is the 2 km crop area around the SIF polygon.
# Since ctr is EPSG:32632, buffer units are meters.
z_buffer_utm <- st_buffer(st_transform(z_poly, crs(ctr)), 2000)

ctr_crop <- crop(ctr, vect(z_buffer_utm))
ctr_crop <- mask(ctr_crop, vect(z_buffer_utm))

ctr_crop_ll <- project(ctr_crop, "EPSG:4326", method = "near")

crop_pal <- colorFactor(
  palette = viridisLite::viridis(nrow(crop_classes)),
  domain = crop_classes$code,
  na.color = "transparent"
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "Basemap") %>%
  addRasterImage(
    ctr_crop_ll,
    colors = crop_pal,
    opacity = 0.75,
    project = FALSE,
    group = "Crop type 2018"
  ) %>%
  addPolygons(
    data = tile_818,
    fill = FALSE,
    color = "firebrick",
    weight = 2,
    group = "MGRS tile",
    label = ~mgrs_tile
  ) %>%
  addPolygons(
    data = z_poly,
    fill = FALSE,
    color = "cyan",
    weight = 3,
    group = "SIF polygon",
    label = ~paste("SIF ID:", sif_id)
  ) %>%
  # addCircleMarkers(
  #   data = z_center,
  #   radius = 5,
  #   color = "black",
  #   fillColor = "yellow",
  #   fillOpacity = 1,
  #   weight = 1,
  #   group = "SIF center",
  #   label = ~paste("SIF ID:", sif_id)
  # ) %>%
  addLegend(
    position = "bottomright",
    pal = crop_pal,
    values = crop_classes$code,
    title = "Crop code",
    opacity = 0.8
  ) %>%
  addLayersControl(
    overlayGroups = c("Crop type 2018", "MGRS tile", "SIF polygon", "SIF center"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  fitBounds(
    lng1 = st_bbox(st_transform(z_buffer_utm, 4326))[["xmin"]],
    lat1 = st_bbox(st_transform(z_buffer_utm, 4326))[["ymin"]],
    lng2 = st_bbox(st_transform(z_buffer_utm, 4326))[["xmax"]],
    lat2 = st_bbox(st_transform(z_buffer_utm, 4326))[["ymax"]]
  )


# Plot of mgrs tiles + sif data points
# if (interactive()) {
#   # Quick plot check: target MGRS tiles and retained Niedersachsen SIF centers.
#   germany_states <- giscoR::gisco_get_nuts(
#     country = "DE",
#     nuts_level = 1,
#     resolution = "03",
#     epsg = 4326
#   ) %>%
#     select(state = NUTS_NAME, geometry)
# 
#   ggplot() +
#     geom_sf(data = germany_states, fill = "grey95", color = "grey40", linewidth = 0.25) +
#     geom_sf(
#       data = mgrs_target,
#       fill = NA,
#       color = "firebrick",
#       linewidth = 0.45
#     ) +
#     geom_sf(
#       data = ns_mgrs_points,
#       size = 0.1,
#       color = "darkgreen"
#     ) +
#     coord_sf(
#       xlim = c(6.5, 11.8),
#       ylim = c(51.2, 54),
#       expand = FALSE
#     )
# }
