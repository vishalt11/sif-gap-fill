# https://excalidraw.com/#json=7ZHLSFM96vFcDuIt015t8,Vab-cWU0LlP5LirkOegdHA
library(ncdf4) 
library(readr)

library(terra)
library(sf) 
library(rnaturalearth)
library(rnaturalearthdata)
library(leaflet)
library(ggrepel)
library(viridis)
library(giscoR)

library(tidyverse)


#-------------------------------------------------------------------------------
# 1 file test
nc_data <- nc_open('../data/oco2_dump//oco2_LtSIF_200404_B11012Ar_220825231256s.nc4')

attributes(nc_data$var)

lon_corners <- ncvar_get(nc_data,"Longitude_Corners")

#-------------------------------------------------------------------------------
# Measurement mode shape analysis

df <- base::readRDS('../data/SIF_v2_corrected_subset.rds')
df <- drop_na(df)

table(df$Metadata.MeasurementMode)
table(df$Quality_Flag)

plot_data <- df %>% filter(Delta_Time >= as.Date("2020-02-01"), Delta_Time <= as.Date("2020-05-31"))

germany <- rnaturalearth::ne_states(country = "Germany", returnclass = "sf")
bavaria <- germany[germany$name_en == "Bavaria",]

plot_sf <- st_as_sf(plot_data, coords = c("Longitude", "Latitude"), crs = 4326)
# crop the sf to include bavaria only 
plot_bavaria <- plot_sf[bavaria,]

table(plot_bavaria$Metadata.MeasurementMode)

plot_bavaria_poly <- plot_bavaria %>%
  st_drop_geometry() %>%
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
            lon1, lat1   # close polygon
          ),
          ncol = 2,
          byrow = TRUE
        )))
      }
    )
  ) %>%
  st_as_sf(crs = 4326)

st_geometry_type(plot_bavaria_poly)
st_crs(plot_bavaria_poly)
head(plot_bavaria_poly, 2)

sort(unique(as.Date(plot_bavaria_poly$Delta_Time)))

one_mode1 <- plot_bavaria_poly %>%
  filter(Metadata.MeasurementMode == 4) %>%
  slice(1)

d <- st_drop_geometry(one_mode1)

corner_df <- tibble(
  corner = paste0("corner", 1:4),
  lon = unname(unlist(d[paste0("Lon_corner", 1:4)])),
  lat = unname(unlist(d[paste0("Lat_corner", 1:4)]))
)

corner_df
corner_df[chull(corner_df$lon, corner_df$lat), ]

plot(st_geometry(one_mode1), col = "blue", main = "Mode 1 footprint")
points(corner_df$lon, corner_df$lat, pch = 19, col = "red")
text(corner_df$lon, corner_df$lat, labels = corner_df$corner)



one_mode1 <- plot_bavaria_poly %>%
  filter(Metadata.MeasurementMode == 0)


ggplot() +
  geom_sf(data = bavaria, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(
    data = one_mode1,
    aes(fill = Daily_SIF_740nm),
    color = "white",
    linewidth = 0.05,
    alpha = 0.8
  ) +
  scale_fill_viridis_c(option = "viridis", na.value = "transparent") +
  coord_sf(
    xlim = c(10.55, 11.05),
    ylim = c(47.45, 47.85),
    expand = FALSE
  ) +
  theme_minimal()

one_mode1 <- one_mode1 %>%
  st_transform(25832) %>%
  mutate(
    area_m2 = as.numeric(st_area(geometry)),
    area_km2 = area_m2 / 1e6
  ) %>%
  st_transform(4326)

summary(one_mode1$area_km2)



