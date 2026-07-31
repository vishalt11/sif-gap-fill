library(sf)
library(dplyr)

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
