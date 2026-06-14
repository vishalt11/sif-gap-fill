library(sf)
library(dplyr)
library(purrr)
library(ggplot2)
library(giscoR)
library(geographiclib)
#library(mgrs)

# German federal states from GISCO
germany_states <- gisco_get_nuts(
  country = "DE",
  nuts_level = 1,
  resolution = "01",
  epsg = 4326
) %>%
  select(state = NUTS_NAME, geometry)

germany_union <- germany_states %>%
  summarise(geometry = st_union(geometry))

# Make candidate lon/lat points over a slightly expanded Germany bbox
bb <- st_bbox(germany_union)

candidate_points <- expand.grid(
  lon = seq(bb["xmin"] - 1, bb["xmax"] + 1, by = 0.05),
  lat = seq(bb["ymin"] - 1, bb["ymax"] + 1, by = 0.05)
)

# Convert candidate points to 100 km MGRS tile IDs
candidate_mgrs <- geographiclib::mgrs_fwd(
  as.matrix(candidate_points[, c("lon", "lat")]),
  precision = 0L
)

candidate_mgrs <- sort(unique(candidate_mgrs))


make_mgrs_100km_tile <- function(code) {
  info <- geographiclib::mgrs_rev(code)
  
  epsg <- as.integer(sub("EPSG:", "", info$crs[1]))
  
  # Snap to the lower-left 100 km UTM grid corner.
  # This avoids depending on whether mgrs_rev(code) returns a corner or center.
  x0 <- floor(info$x[1] / 100000) * 100000
  y0 <- floor(info$y[1] / 100000) * 100000
  
  coords <- matrix(
    c(
      x0,          y0,
      x0 + 100000, y0,
      x0 + 100000, y0 + 100000,
      x0,          y0 + 100000,
      x0,          y0
    ),
    ncol = 2,
    byrow = TRUE
  )
  
  st_sf(
    mgrs_tile = code,
    grid_zone = info$grid_zone[1],
    square_100km = info$square_100km[1],
    zone = info$zone[1],
    crs_utm = info$crs[1],
    geometry = st_sfc(st_polygon(list(coords)), crs = epsg)
  ) %>%
    st_transform(4326)
}

mgrs_tiles_all <- map_dfr(candidate_mgrs, make_mgrs_100km_tile)

mgrs_tiles_germany <- mgrs_tiles_all %>%
  st_filter(germany_union, .predicate = st_intersects)

#only 32U tiles
mgrs_tiles_germany <- mgrs_tiles_germany %>% filter(grid_zone %in% c('32U'))

saveRDS(mgrs_tiles_germany, '../data/mgrs_de.rds')

ggplot() +
  geom_sf(
    data = mgrs_tiles_germany,
    fill = NA,
    color = "firebrick",
    linewidth = 0.45
  ) +
  geom_sf(
    data = germany_states,
    fill = NA,
    color = "grey25",
    linewidth = 0.25
  ) +
  geom_sf_text(
    data = st_point_on_surface(mgrs_tiles_germany),
    aes(label = mgrs_tile),
    size = 3,
    color = "firebrick"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal()



# mgrs_tiles_germany_clipped <- st_intersection(
#   mgrs_tiles_germany,
#   germany_union
# )
# 
# ggplot() +
#   geom_sf(
#     data = mgrs_tiles_germany_clipped,
#     fill = NA,
#     color = "firebrick",
#     linewidth = 0.45
#   ) +
#   geom_sf(
#     data = germany_states,
#     fill = NA,
#     color = "grey25",
#     linewidth = 0.25
#   ) +
#   coord_sf(expand = FALSE) +
#   theme_minimal()
