library(sf)
library(terra)
library(raster)
library(leaflet)
library(htmlwidgets)

#-------------------------------------------------------------------------------
# Compare annual crop-type changes inside the target 32UPC area
#-------------------------------------------------------------------------------

crop_dir <- "data/crop_type_tif"
kml_path <- "data/target_ns_32upc.kml"
legend_path <- file.path(crop_dir, "LEGEND_CropTypes.txt")
output_html <- "data/croptype_change_32upc_leaflet.html"

# Use the KML as a bounding box, as requested. Set TRUE if you later want the
# crop rasters clipped to the exact KML polygon instead of its rectangular extent.
MASK_TO_KML_POLYGON <- FALSE

stopifnot(file.exists(kml_path))
stopifnot(file.exists(legend_path))

crop_files <- list.files(
  crop_dir,
  pattern = "^croptypes_[0-9]{4}\\.tif$",
  full.names = TRUE
)

if (length(crop_files) == 0) {
  stop("No croptypes_YYYY.tif files found in: ", crop_dir)
}

years <- sub("^croptypes_([0-9]{4})\\.tif$", "\\1", basename(crop_files))
file_order <- order(years)
crop_files <- crop_files[file_order]
years <- years[file_order]
year_groups <- paste("Crop type", years)

crop_classes <- read.delim(legend_path, sep = "\t", stringsAsFactors = FALSE)
names(crop_classes) <- c("code", "label")
crop_classes$code <- as.integer(crop_classes$code)

target_raw <- sf::st_read(kml_path, quiet = TRUE)
target_raw <- sf::st_make_valid(target_raw)
target_wgs84 <- sf::st_transform(target_raw, 4326)
target_outline <- sf::st_as_sf(sf::st_union(target_wgs84))
target_outline$name <- "target_ns_32upc"

legend_classes <- crop_classes[crop_classes$code != 0, ]

fallback_palette <- grDevices::hcl.colors(nrow(legend_classes), "Dark 3")
names(fallback_palette) <- as.character(legend_classes$code)

semantic_palette <- c(
  #"11" = "#f2c94c",  # winter_wheat
  "11" = "blue",  # winter_wheat
  "12" = "#c9a227",  # winter_barley
  "13" = "#b08935",  # winter_rye
  "14" = "#d7b56d",  # other_winter_cereals
  "21" = "#ffe08a",  # spring_wheat
  "22" = "#bf8f2f",  # spring_barley
  "23" = "#e6cc98",  # spring_oat
  "30" = "#2f9e44",  # maize
  "40" = "#20c997",  # legumes
  "50" = "#8d6e63",  # potato
  "60" = "#e64980",  # sugar_beet
  "71" = "#ffd43b",  # rapeseed
  "81" = "#51cf66",  # clover/alfalfa
  "82" = "#69db7c",  # arable_grass
  "83" = "#2b8a3e",  # permanent_grassland
  "90" = "#9c36b5",  # vineyard
  "100" = "#e03131", # fruit_trees_and_other_woody_vegetation
  "110" = "#087f5b", # hops
  "111" = "#adb5bd"  # other_agricultural_use
)

crop_palette <- fallback_palette
matched_palette_codes <- intersect(names(crop_palette), names(semantic_palette))
crop_palette[matched_palette_codes] <- semantic_palette[matched_palette_codes]

pal <- leaflet::colorFactor(
  palette = unname(crop_palette),
  domain = as.integer(names(crop_palette)),
  na.color = "transparent"
)

crop_to_leaflet_raster <- function(path) {
  year <- sub("^croptypes_([0-9]{4})\\.tif$", "\\1", basename(path))

  crop_raster <- terra::rast(path)
  levels(crop_raster) <- data.frame(
    value = crop_classes$code,
    crop = crop_classes$label
  )

  target_in_raster_crs <- sf::st_transform(target_raw, terra::crs(crop_raster))
  target_bbox <- sf::st_as_sfc(
    sf::st_bbox(target_in_raster_crs),
    crs = sf::st_crs(target_in_raster_crs)
  )

  crop_target <- terra::crop(
    crop_raster,
    terra::vect(sf::st_as_sf(target_bbox)),
    snap = "out"
  )

  if (MASK_TO_KML_POLYGON) {
    crop_target <- terra::mask(
      crop_target,
      terra::vect(target_in_raster_crs)
    )
  }

  crop_target[crop_target == 0] <- NA

  crop_target_wgs84 <- terra::project(
    crop_target,
    "EPSG:4326",
    method = "near"
  )

  names(crop_target_wgs84) <- paste0("crop_type_", year)
  as(crop_target_wgs84, "Raster")
}

leaflet_rasters <- stats::setNames(
  lapply(crop_files, crop_to_leaflet_raster),
  year_groups
)

croptype_change_map <- leaflet::leaflet(
  options = leaflet::leafletOptions(preferCanvas = TRUE)
)

croptype_change_map <- leaflet::addProviderTiles(
  croptype_change_map,
  leaflet::providers$Esri.WorldImagery,
  group = "Esri World Imagery"
)

for (group_name in names(leaflet_rasters)) {
  croptype_change_map <- leaflet::addRasterImage(
    croptype_change_map,
    leaflet_rasters[[group_name]],
    colors = pal,
    opacity = 0.72,
    group = group_name,
    project = FALSE,
    maxBytes = 64 * 1024 * 1024
  )
}

croptype_change_map <- leaflet::addPolygons(
  croptype_change_map,
  data = target_outline,
  color = "#00d5ff",
  weight = 2,
  opacity = 1,
  fill = FALSE,
  group = "KML bbox"
)

croptype_change_map <- leaflet::addLegend(
  croptype_change_map,
  position = "bottomright",
  colors = unname(crop_palette[as.character(legend_classes$code)]),
  labels = gsub("_", " ", legend_classes$label),
  opacity = 0.85,
  title = "Crop type"
)

croptype_change_map <- leaflet::addLayersControl(
  croptype_change_map,
  baseGroups = "Esri World Imagery",
  overlayGroups = c(year_groups, "KML bbox"),
  options = leaflet::layersControlOptions(collapsed = FALSE)
)

if (length(year_groups) > 1) {
  croptype_change_map <- leaflet::hideGroup(
    croptype_change_map,
    year_groups[-length(year_groups)]
  )
}

target_bbox_wgs84 <- sf::st_bbox(target_outline)
croptype_change_map <- leaflet::fitBounds(
  croptype_change_map,
  lng1 = target_bbox_wgs84[["xmin"]],
  lat1 = target_bbox_wgs84[["ymin"]],
  lng2 = target_bbox_wgs84[["xmax"]],
  lat2 = target_bbox_wgs84[["ymax"]]
)

htmlwidgets::saveWidget(
  croptype_change_map,
  file = output_html,
  selfcontained = FALSE
)

croptype_change_map
