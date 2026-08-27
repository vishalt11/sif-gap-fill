library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(furrr)
library(future)

input_file <- 'data/base_sif/sif_sf_1_12_cleaned.rds'
output_file <- "data/sif_sf_1_12_crop_zonal_18_24.rds"
output_csv <- "data/sif_sf_1_12_crop_zonal_18_24.csv"

crop_type_dir <- "data/crop_type_tif"
parallel_workers <- 2
chunk_grid_size_m <- 25000
buffer_cells <- 2

corner_cols <- c(
  "Lat_corner1", "Lat_corner2", "Lat_corner3", "Lat_corner4",
  "Lon_corner1", "Lon_corner2", "Lon_corner3", "Lon_corner4"
)

df <- readRDS(input_file)
df <- df %>% st_drop_geometry()

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

build_sif_sf <- function(data) {
  if (inherits(data, "sf")) {
    sif_sf <- st_as_sf(data)

    if (is.na(st_crs(sif_sf))) {
      st_crs(sif_sf) <- 4326
    }

    return(st_make_valid(sif_sf))
  }

  missing_corner_cols <- setdiff(corner_cols, names(data))
  if (length(missing_corner_cols) > 0) {
    stop("Missing corner columns: ", paste(missing_corner_cols, collapse = ", "))
  }

  data %>%
    mutate(across(all_of(corner_cols), as.numeric)) %>%
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
}

empty_crop_count_table <- function(row_ids, crop_classes, fill_value = NA_integer_) {
  out <- tibble(.extract_row_id = row_ids)

  for (count_col in crop_classes$count_col) {
    out[[count_col]] <- rep(fill_value, length(row_ids))
  }

  out %>%
    mutate(
      crop_pixel_count = rowSums(
        across(all_of(crop_classes$count_col)),
        na.rm = FALSE
      )
    ) %>%
    relocate(crop_pixel_count, .after = .extract_row_id)
}

buffered_polygon_extent <- function(raster, polygons_same_crs, buffer_cells = 2) {
  polygons_vect <- vect(st_as_sf(polygons_same_crs))
  polygon_extent <- ext(polygons_vect)
  raster_extent <- ext(raster)

  x_buffer <- abs(xres(raster)) * buffer_cells
  y_buffer <- abs(yres(raster)) * buffer_cells

  crop_xmin <- max(xmin(polygon_extent) - x_buffer, xmin(raster_extent))
  crop_xmax <- min(xmax(polygon_extent) + x_buffer, xmax(raster_extent))
  crop_ymin <- max(ymin(polygon_extent) - y_buffer, ymin(raster_extent))
  crop_ymax <- min(ymax(polygon_extent) + y_buffer, ymax(raster_extent))

  if (crop_xmin >= crop_xmax || crop_ymin >= crop_ymax) {
    return(NULL)
  }

  ext(crop_xmin, crop_xmax, crop_ymin, crop_ymax)
}

crop_raster_to_polygon_extent <- function(raster, polygons_same_crs, buffer_cells = 2) {
  crop_extent <- buffered_polygon_extent(
    raster,
    polygons_same_crs,
    buffer_cells = buffer_cells
  )

  if (is.null(crop_extent)) {
    return(NULL)
  }

  tryCatch(
    crop(raster, crop_extent),
    error = function(e) {
      warning("Raster crop failed; using full raster. Error: ", conditionMessage(e))
      raster
    }
  )
}

add_spatial_chunk_id <- function(polygons_same_crs, chunk_grid_size_m) {
  centroid_xy <- st_coordinates(st_centroid(st_geometry(polygons_same_crs)))

  polygons_same_crs %>%
    mutate(
      .chunk_x = floor(centroid_xy[, "X"] / chunk_grid_size_m),
      .chunk_y = floor(centroid_xy[, "Y"] / chunk_grid_size_m),
      .chunk_id = paste(.chunk_x, .chunk_y, sep = "_")
    )
}

count_crop_pixels_in_chunk <- function(
    raster_file,
    polygons_chunk,
    crop_classes,
    buffer_cells = 2) {
  row_ids <- polygons_chunk$.extract_row_id

  crop_raster <- rast(raster_file)
  levels(crop_raster) <- data.frame(
    value = crop_classes$code,
    crop = crop_classes$label
  )

  crop_raster_chunk <- crop_raster_to_polygon_extent(
    crop_raster,
    polygons_chunk,
    buffer_cells = buffer_cells
  )

  if (is.null(crop_raster_chunk)) {
    return(empty_crop_count_table(row_ids, crop_classes, fill_value = 0L))
  }

  extracted <- terra::extract(crop_raster_chunk, vect(polygons_chunk))

  if (nrow(extracted) == 0) {
    return(empty_crop_count_table(row_ids, crop_classes, fill_value = 0L))
  }

  value_col <- setdiff(names(extracted), "ID")[1]

  if (is.na(value_col)) {
    return(empty_crop_count_table(row_ids, crop_classes, fill_value = 0L))
  }

  id_lookup <- polygons_chunk %>%
    st_drop_geometry() %>%
    mutate(ID = row_number()) %>%
    select(ID, .extract_row_id)

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
    filter(!is.na(crop_code), crop_code %in% crop_classes$code) %>%
    left_join(id_lookup, by = "ID") %>%
    filter(!is.na(.extract_row_id)) %>%
    count(.extract_row_id, crop_code, name = "n")

  expand_grid(
    .extract_row_id = row_ids,
    crop_code = crop_classes$code
  ) %>%
    left_join(counts_long, by = c(".extract_row_id", "crop_code")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    left_join(
      crop_classes %>% select(code, count_col),
      by = c("crop_code" = "code")
    ) %>%
    select(.extract_row_id, count_col, n) %>%
    pivot_wider(
      names_from = count_col,
      values_from = n,
      values_fill = list(n = 0L),
      values_fn = list(n = sum)
    ) %>%
    mutate(
      crop_pixel_count = rowSums(across(all_of(crop_classes$count_col)))
    ) %>%
    relocate(crop_pixel_count, .after = .extract_row_id)
}

count_crop_pixels_for_year <- function(year, polygons_year, crop_classes) {
  raster_file <- file.path(crop_type_dir, paste0("croptypes_", year, ".tif"))

  if (!file.exists(raster_file)) {
    warning("No crop type raster found for ", year, ": ", raster_file)
    return(empty_crop_count_table(
      polygons_year$.extract_row_id,
      crop_classes,
      fill_value = NA_integer_
    ))
  }

  message("Preparing crop chunks for ", year, " with ", nrow(polygons_year), " SIF polygons.")

  crop_template <- rast(raster_file)
  polygons_crop_crs <- polygons_year %>%
    st_transform(crs(crop_template)) %>%
    add_spatial_chunk_id(chunk_grid_size_m = chunk_grid_size_m)

  polygon_chunks <- split(polygons_crop_crs, polygons_crop_crs$.chunk_id)

  message(
    "Processing ",
    length(polygon_chunks),
    " spatial chunks for crop raster ",
    basename(raster_file)
  )

  future_imap_dfr(
    polygon_chunks,
    function(polygons_chunk, chunk_id) {
      message(
        "Extracting crop pixels for ",
        year,
        " chunk ",
        chunk_id,
        " (",
        nrow(polygons_chunk),
        " polygons)."
      )

      count_crop_pixels_in_chunk(
        raster_file = raster_file,
        polygons_chunk = polygons_chunk,
        crop_classes = crop_classes,
        buffer_cells = buffer_cells
      )
    },
    .progress = TRUE,
    .options = furrr_options(
      seed = FALSE,
      scheduling = 1,
      packages = c(
        "dplyr",
        "tidyr",
        "purrr",
        "tibble",
        "sf",
        "terra",
        "stringr",
        "magrittr"
      )
    )
  )
}

message("Preparing SIF polygons.")

sif_sf <- build_sif_sf(df) %>%
  mutate(.extract_row_id = row_number())

if (!"sif_row_id" %in% names(sif_sf)) {
  sif_sf <- sif_sf %>%
    mutate(sif_row_id = .extract_row_id)
}

if (!"sif_year" %in% names(sif_sf)) {
  if (!"Delta_Date" %in% names(sif_sf)) {
    stop("Missing both sif_year and Delta_Date. One is needed to choose crop rasters.")
  }

  sif_sf <- sif_sf %>%
    mutate(sif_year = year(as.Date(Delta_Date)))
} else {
  sif_sf <- sif_sf %>%
    mutate(sif_year = as.integer(sif_year))
}

extractable_sif_sf <- sif_sf %>%
  filter(!is.na(sif_year))

extractable_sif_sf <- extractable_sif_sf[
  !st_is_empty(st_geometry(extractable_sif_sf)),
]

extract_years <- extractable_sif_sf %>%
  st_drop_geometry() %>%
  distinct(sif_year) %>%
  arrange(sif_year) %>%
  pull(sif_year)

if (length(extract_years) == 0) {
  crop_counts <- empty_crop_count_table(integer(), crop_classes, fill_value = NA_integer_)
} else {
  plan(multisession, workers = parallel_workers)

  crop_counts <- map_dfr(extract_years, function(target_year) {
    polygons_year <- extractable_sif_sf %>%
      filter(sif_year == target_year)

    count_crop_pixels_for_year(
      year = target_year,
      polygons_year = polygons_year,
      crop_classes = crop_classes
    )
  })

  plan(sequential)
}

crop_output_cols <- c(
  "crop_pixel_count",
  crop_classes$count_col,
  "ww_pct"
)

sif_sf <- sif_sf %>%
  select(-any_of(crop_output_cols)) %>%
  left_join(crop_counts, by = ".extract_row_id") %>%
  mutate(
    ww_pct = if_else(
      !is.na(crop_pixel_count) & crop_pixel_count > 0,
      crop_count_winter_wheat / crop_pixel_count,
      NA_real_
    )
  ) %>%
  select(-.extract_row_id)

message("Input rows: ", nrow(df))
message("Rows eligible for crop extraction: ", nrow(extractable_sif_sf))
message("Rows with crop pixels: ", sum(!is.na(sif_sf$crop_pixel_count) & sif_sf$crop_pixel_count > 0))

active_growth_months <- list(
  crop_count_winter_wheat = c(2:7, 10:11),
  crop_count_winter_barley = c(2:6, 10:11),
  crop_count_winter_rye = c(2:7, 10:11),
  crop_count_other_winter_cereals = c(2:7, 10:11),
  crop_count_spring_wheat = 3:8,
  crop_count_spring_barley = 3:8,
  crop_count_spring_oat = 3:8,
  crop_count_maize = 5:10,
  crop_count_legumes = 4:9,
  crop_count_potato = 4:9,
  crop_count_sugar_beet = 4:10,
  crop_count_rapeseed = c(2:7, 9:11),
  crop_count_clover_alfalfa = 3:10,
  crop_count_arable_grass = 3:11,
  crop_count_permanent_grassland = 3:11,
  crop_count_vineyard = 4:10,
  crop_count_fruit_trees_and_other_woody_vegetation = 3:10,
  crop_count_hops = 4:9,
  crop_count_other_agricultural_use = 3:10
)

add_active_growth_pct <- function(data, active_growth_months) {
  data_attrs <- if (inherits(data, "sf")) {
    st_drop_geometry(data)
  } else {
    data
  }

  if ("Delta_Date" %in% names(data_attrs)) {
    sif_month <- month(as.Date(substr(as.character(data_attrs$Delta_Date), 1, 10)))
  } else if ("Delta_Time" %in% names(data_attrs)) {
    sif_month <- month(as.Date(substr(as.character(data_attrs$Delta_Time), 1, 10)))
  } else {
    stop("Missing both Delta_Date and Delta_Time. One is needed for active growth months.")
  }

  active_crop_cols <- intersect(names(active_growth_months), names(data_attrs))
  active_growth_pixel_count <- rep(NA_real_, nrow(data_attrs))

  for (target_month in sort(unique(sif_month[!is.na(sif_month)]))) {
    month_rows <- which(sif_month == target_month)
    month_crop_cols <- active_crop_cols[
      map_lgl(active_growth_months[active_crop_cols], ~ target_month %in% .x)
    ]

    if (length(month_crop_cols) == 0) {
      active_growth_pixel_count[month_rows] <- 0
    } else {
      active_growth_pixel_count[month_rows] <- rowSums(
        as.matrix(data_attrs[month_rows, month_crop_cols, drop = FALSE]),
        na.rm = FALSE
      )
    }
  }

  data %>%
    mutate(
      active_growth_pixel_count = active_growth_pixel_count,
      active_growth_pct = if_else(
        !is.na(crop_pixel_count) & crop_pixel_count > 0,
        active_growth_pixel_count / crop_pixel_count,
        NA_real_
      )
    )
}

sif_sf <- add_active_growth_pct(sif_sf, active_growth_months)

sif_sf[is.na(sif_sf$ww_pct),]$ww_pct <- 0
sif_sf[is.na(sif_sf$active_growth_pct),]$active_growth_pct <- 0

saveRDS(sif_sf, output_file)

sif_sf %>%
  st_drop_geometry() %>%
  write_csv(output_csv)

message("Saved crop composition RDS to ", output_file)
message("Saved crop composition CSV to ", output_csv)


unique(lubridate::year(sif_sf$Delta_Date))
sif_sf <- sif_sf[lubridate::year(sif_sf$Delta_Date) %in% c('2019', '2020', '2021', '2022', '2023', '2024'),]

saveRDS(sif_sf, 'data/sif_sf_1_12_crop_zonal_19_24.rds')
sif_sf %>%
  st_drop_geometry() %>%
  write_csv('data/sif_sf_1_12_crop_zonal_19_24.csv')
