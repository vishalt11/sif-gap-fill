library(tidyverse)
library(terra)
library(sf)


#-------------------------------------------------------------------------------
# Find high density wheat farm clusters and focus on these regions
#-------------------------------------------------------------------------------

ctr <- rast('data/crop_type_tif/croptypes_2024.tif')

crop_classes <- readr::read_delim("data/crop_type_tif/LEGEND_CropTypes.txt", delim = "\t")
colnames(crop_classes) <- c("code", "label")
levels(ctr) <- data.frame(value = crop_classes$code, crop = crop_classes$label)

nuts1_de <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 1, resolution = "01",epsg = 4326) %>% select(NUTS_NAME, NUTS_ID, geometry)
nuts1_de  <- subset(nuts1_de, grepl("^DEA", NUTS_ID))

nuts3_de <- giscoR::gisco_get_nuts(country = "DE", nuts_level = 3, resolution = "01",epsg = 4326) %>% select(NUTS_NAME, NUTS_ID, geometry)

bav_utm <- st_transform(nuts1_de, crs(ctr))
bav_utm <- vect(bav_utm)

ctr_bav <- crop(ctr, bav_utm)
ctr_bav <- mask(ctr_bav, bav_utm)

writeRaster(ctr_bav,"data/crop_type_tif/croptypes_2024_nw.tif", overwrite = TRUE)

#-----
ctr_bav <- rast('../data/crop_type_tif/croptypes_2024_nw.tif')
levels(ctr_bav)
plot(ctr_bav)

ctr_bav_1km <- terra::aggregate(ctr_bav, fact = 100, fun  = 'modal', na.rm = TRUE)
res(ctr_bav_1km)

nuts3_utm <- st_transform(nuts3_de, crs(ctr))
nuts3_utm <- vect(nuts3_utm)

plot(ctr_bav_1km)
lines(nuts3_utm, col = "black", lwd = 2)
#-----

wheat <- ctr_bav == 11
wheat_1km  <- terra::aggregate(wheat, fact = 100, fun  = 'mean', na.rm = TRUE)
high_wheat <- wheat_1km >= 0.4   #  >40% wheat

plot(wheat_1km)
plot(high_wheat)

writeRaster(high_wheat, 'data/nw_ww_gt40.tif')

r <- rast('data/tifs/bavaria_ww_gt40.tif')

plot(r)
#-------------------------------------------------------------------------------
# Where are spatially coherent, high-density wheat systems that form real farming regions?
# High Density Wheat farm cluster
#-------------------------------------------------------------------------------

library(spdep)

# points at cell centers = one feature per cell (no dissolving)
pts <- as.points(wheat_1km, na.rm=TRUE) |> st_as_sf()
names(pts)[names(pts) == names(pts)[1]] <- "wheat_share"
pts$wheat_share <- as.numeric(pts$wheat_share)

# k-nearest neighbors (good for regular grids)
coords <- st_coordinates(pts)
nb <- knn2nb(knearneigh(coords, k = 8))
lw <- nb2listw(nb, style="W")

mi <- localmoran(pts$wheat_share, lw)

pts$Ii   <- mi[, "Ii"]
pts$Z    <- mi[, "Z.Ii"]
pts$pval <- mi[, "Pr(z != E(Ii))"]

pts$lisa <- "Not significant"
pts$lisa[pts$wheat_share >= 0.6 & pts$Ii > 0 & pts$pval < 0.05] <- "High–High"
pts$lisa[pts$wheat_share <  0.4 & pts$Ii > 0 & pts$pval < 0.05] <- "Low–Low"

mean(pts$Ii > 0 & pts$pval < 0.05, na.rm = TRUE)

hotspots <- pts[pts$wheat_share >= 0.6 & pts$Ii > 0 & pts$pval < 0.05,]
plot(hotspots, max.plot = 1)

#Red → high–high clusters (wheat hotspots)
#Blue → low–low clusters
#Near zero → spatially random


ggplot(pts) +
  geom_sf(aes(color = Z), size = 0.6) +
  geom_sf(data = nuts3_de, fill = NA, color = "black", linewidth = 0.3) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  labs(color = "Local Moran Z") +
  theme_bw()

# With LISA High–High, you are explicitly keeping the cells that satisfy all three:

# High wheat density (your threshold, e.g. ≥ 60%)

# Surrounded by high wheat neighbors (positive local Moran’s I)

# Statistically significant (p-value threshold, e.g. 0.05)
# Conceptually: “This cell has a lot of wheat, and so do the cells around it — more than expected by chance.”


saveRDS(pts, '../data/local_Moran_wheatclusters.rds')



hh_pts <- pts %>%
  mutate(hh = wheat_share >= 0.6 & Ii > 0 & pval < 0.05)

nuts3_de <- st_transform(nuts3_de, st_crs(hh_pts))

pts_n3 <- st_join(
  hh_pts,
  nuts3_de[, c("NUTS_ID", "NUTS_NAME")],
  join = st_within,
  left = FALSE
)

# 3) summarize per district
nuts3_stats <- pts_n3 %>%
  st_drop_geometry() %>%
  group_by(NUTS_ID, NUTS_NAME) %>%
  summarise(
    n_cells_total = n(),
    n_cells_hh    = sum(hh, na.rm = TRUE),
    hh_share      = n_cells_hh / n_cells_total,
    hh_area_km2   = n_cells_hh * 1  # 1 km² per 1km cell
  ) %>%
  arrange(desc(hh_share))

nuts3_map <- nuts3_de %>% left_join(nuts3_stats, by = c("NUTS_ID","NUTS_NAME"))

saveRDS(nuts3_map, '../data/nuts3_local_moran_HH_summary.rds')