library(tidyverse)
library(sf)
library(terra)
library(leaflet)
library(giscoR)

df <- readRDS("data/ns_sif_mgrs_crop_composition.rds") %>%
  st_drop_geometry() %>%
  mutate(Delta_Date = as.Date(Delta_Date)) %>%
  filter(format(Delta_Date, "%Y-%m") == "2020-06")

df <- df %>% filter(ww_pct > 0.5)

sif_polygons_2020 <- df %>%
  filter(
    if_all(
      c(starts_with("Lat_corner"), starts_with("Lon_corner"), Daily_SIF_740nm),
      ~ !is.na(.x)
    )
  ) %>%
  mutate(
    geometry = pmap(
      list(
        Lon_corner1, Lat_corner1,
        Lon_corner2, Lat_corner2,
        Lon_corner3, Lat_corner3,
        Lon_corner4, Lat_corner4
      ),
      function(lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4) {
        st_polygon(list(matrix(
          c(
            lon1, lat1,
            lon2, lat2,
            lon3, lat3,
            lon4, lat4,
            lon1, lat1
          ),
          ncol = 2,
          byrow = TRUE
        )))
      }
    )
  ) %>%
  st_as_sf(crs = 4326) %>%
  st_make_valid()

target_ns_32upc <- st_read("data/target_ns_32upc.kml", quiet = TRUE) %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_transform(4326) %>%
  st_make_valid()

sif_polygons_2020 <- sif_polygons_2020[
  lengths(st_intersects(sif_polygons_2020, target_ns_32upc)) > 0,
]

niedersachsen <- gisco_get_nuts(
  country = "DE",
  nuts_level = 1,
  resolution = "01",
  epsg = 4326
) %>%
  filter(NUTS_ID == "DE9") %>%
  select(NUTS_NAME, NUTS_ID, geometry)

crop_classes <- readr::read_delim(
  "data/crop_type_tif/LEGEND_CropTypes.txt",
  delim = "\t",
  show_col_types = FALSE
)
colnames(crop_classes) <- c("code", "label")

crop_raster_2020 <- rast("data/crop_type_tif/croptypes_2020.tif")
levels(crop_raster_2020) <- data.frame(
  value = crop_classes$code,
  crop = crop_classes$label
)

target_crop_crs <- target_ns_32upc %>%
  st_transform(crs(crop_raster_2020)) %>%
  vect()

crop_raster_target <- crop(crop_raster_2020, target_crop_crs) %>%
  mask(target_crop_crs)

crop_raster_leaflet <- raster::raster(crop_raster_target)

wasp_product_id <- "SENTINEL2X_20200615-000000-000_L3A_T32UPC_C_V4-0"
wasp_product_dir <- file.path(
  "data/geodes_wasp_zips/32UPC/2020",
  wasp_product_id
)

reflectance_quantification_value <- 10000
reflectance_nodata <- -10000
land_flag_value <- 4

band_ids_10m <- c("B2", "B3", "B4", "B8")
band_ids_20m <- c("B5", "B6", "B7", "B8A", "B11", "B12")
band_ids <- c(band_ids_10m, band_ids_20m)

wasp_file <- function(type, band_or_res = NULL) {
  filename <- if (is.null(band_or_res)) {
    paste0(wasp_product_id, "_", type, ".tif")
  } else {
    paste0(wasp_product_id, "_", type, "_", band_or_res, ".tif")
  }

  subdir <- if (type %in% c("FLG", "WGT", "DTS")) "MASKS" else ""
  file.path(wasp_product_dir, subdir, filename)
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path, call. = FALSE)
  }

  path
}

normalized_difference <- function(a, b) {
  denominator <- a + b
  ifel(denominator != 0, (a - b) / denominator, NA)
}

load_scaled_land_band <- function(band_id, flg, target_vect) {
  band <- rast(check_file_exists(wasp_file("FRC", band_id)))
  band <- crop(band, target_vect)

  flg_crop <- crop(flg, band) %>%
    resample(band, method = "near")

  band <- ifel(
    flg_crop == land_flag_value & band != reflectance_nodata,
    band / reflectance_quantification_value,
    NA
  )

  mask(band, target_vect)
}

upsample_to_10m <- function(band_20m, template_10m, target_vect) {
  band_10m <- resample(band_20m, template_10m, method = "bilinear")
  mask(band_10m, target_vect)
}

flg_r1 <- rast(check_file_exists(wasp_file("FLG", "R1")))
flg_r2 <- rast(check_file_exists(wasp_file("FLG", "R2")))

target_wasp_crs <- target_ns_32upc %>%
  st_transform(crs(flg_r1)) %>%
  vect()

bands_10m <- band_ids_10m %>%
  set_names() %>%
  map(load_scaled_land_band, flg = flg_r1, target_vect = target_wasp_crs)

template_10m <- bands_10m$B8

bands_20m_upsampled <- band_ids_20m %>%
  set_names() %>%
  map(
    ~ load_scaled_land_band(.x, flg = flg_r2, target_vect = target_wasp_crs) %>%
      upsample_to_10m(template_10m, target_wasp_crs)
  )

band_layers_10m <- c(bands_10m, bands_20m_upsampled)

index_layers_10m <- list(
  ndvi = normalized_difference(band_layers_10m$B8, band_layers_10m$B4),
  ndre = normalized_difference(band_layers_10m$B8, band_layers_10m$B5),
  ndre8a = normalized_difference(band_layers_10m$B8A, band_layers_10m$B6),
  psri = ifel(
    band_layers_10m$B6 != 0,
    (band_layers_10m$B4 - band_layers_10m$B2) / band_layers_10m$B6,
    NA
  ),
  osavi = ifel(
    band_layers_10m$B8 + band_layers_10m$B4 + 0.16 != 0,
    1.16 * (band_layers_10m$B8 - band_layers_10m$B4) /
      (band_layers_10m$B8 + band_layers_10m$B4 + 0.16),
    NA
  ),
  ndwi = normalized_difference(band_layers_10m$B3, band_layers_10m$B8),
  nirv = normalized_difference(band_layers_10m$B8, band_layers_10m$B4) *
    band_layers_10m$B8,
  tcari = ifel(
    band_layers_10m$B4 != 0,
    3 * (
      (band_layers_10m$B5 - band_layers_10m$B4) -
        0.2 * (band_layers_10m$B5 - band_layers_10m$B3) *
          (band_layers_10m$B5 / band_layers_10m$B4)
    ),
    NA
  )
)

sif_palette <- colorNumeric(
  palette = grDevices::hcl.colors(256, "viridis"),
  domain = sif_polygons_2020$Daily_SIF_740nm,
  na.color = "#8f8f8f"
)

crop_palette <- colorFactor(
  palette = grDevices::hcl.colors(nrow(crop_classes), "Spectral"),
  domain = crop_classes$code,
  na.color = "transparent"
)

sif_labels <- sprintf(
  "<strong>SIF:</strong> %.4f<br/><strong>Date:</strong> %s<br/><strong>MGRS:</strong> %s<br/><strong>WW share:</strong> %.1f%%",
  sif_polygons_2020$Daily_SIF_740nm,
  format(sif_polygons_2020$Delta_Date, "%Y-%m-%d"),
  sif_polygons_2020$mgrs_tile,
  100 * sif_polygons_2020$ww_pct
) %>%
  lapply(htmltools::HTML)

target_bbox <- st_bbox(target_ns_32upc)

add_continuous_raster <- function(map, raster_layer, group, palette, opacity) {
  raster_values <- terra::values(raster_layer, mat = FALSE)
  raster_palette <- colorNumeric(
    palette = palette,
    domain = raster_values,
    na.color = "transparent"
  )

  map %>%
    addRasterImage(
      raster::raster(raster_layer),
      colors = raster_palette,
      opacity = opacity,
      project = TRUE,
      maxBytes = 50 * 1024 * 1024,
      group = group
    )
}

band_layer_groups <- paste("Band", names(band_layers_10m))
index_layer_groups <- toupper(names(index_layers_10m))

ns_sif_leaflet <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addRasterImage(
    crop_raster_leaflet,
    colors = crop_palette,
    opacity = 0.45,
    project = TRUE,
    maxBytes = 50 * 1024 * 1024,
    group = "Crop types 2020"
  )

for (band_id in names(band_layers_10m)) {
  ns_sif_leaflet <- add_continuous_raster(
    ns_sif_leaflet,
    band_layers_10m[[band_id]],
    group = paste("Band", band_id),
    palette = grDevices::gray.colors(256, start = 0, end = 1),
    opacity = 0.65
  )
}

for (index_id in names(index_layers_10m)) {
  ns_sif_leaflet <- add_continuous_raster(
    ns_sif_leaflet,
    index_layers_10m[[index_id]],
    group = toupper(index_id),
    palette = grDevices::hcl.colors(256, "RdYlGn"),
    opacity = 0.65
  )
}

ns_sif_leaflet <- ns_sif_leaflet %>%
  addPolygons(
    data = sif_polygons_2020,
    fillColor = ~ sif_palette(Daily_SIF_740nm),
    fillOpacity = 0.75,
    color = "#ffffff",
    weight = 0.6,
    opacity = 0.9,
    label = sif_labels,
    group = "2020 SIF polygons"
  ) %>%
  addPolygons(
    data = target_ns_32upc,
    fill = FALSE,
    color = "#2dc0fb",
    weight = 3,
    opacity = 1,
    group = "Target KML"
  ) %>%
  addPolygons(
    data = niedersachsen,
    fill = FALSE,
    color = "#202020",
    weight = 2,
    opacity = 1,
    group = "Niedersachsen border"
  ) %>%
  addLegend(
    pal = sif_palette,
    values = sif_polygons_2020$Daily_SIF_740nm,
    title = "Daily SIF 740 nm",
    position = "bottomright"
  ) %>%
  addLayersControl(
    overlayGroups = c(
      "Crop types 2020",
      band_layer_groups,
      index_layer_groups,
      "2020 SIF polygons",
      "Target KML",
      "Niedersachsen border"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup(c(band_layer_groups, index_layer_groups)) %>%
  addScaleBar(position = "bottomleft") %>%
  fitBounds(
    lng1 = target_bbox[["xmin"]],
    lat1 = target_bbox[["ymin"]],
    lng2 = target_bbox[["xmax"]],
    lat2 = target_bbox[["ymax"]]
  )

ns_sif_leaflet

