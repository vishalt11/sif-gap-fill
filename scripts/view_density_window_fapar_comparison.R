# Compare SNAP 20 m, SNAP averaged onto the GLASS grid, and GLASS nominal 250 m.
# Run with source(), or call plot_window_fapar("YOUR_WINDOW_ID") after sourcing.
# No downloaded rasters are modified; no new raster files are written.

library(terra)
library(sf)
library(ggplot2)
library(patchwork)

ROOT <- "D:/UE/4_Semester/code/scripts"
WINDOW_ID <- "s2_density_4000m_T32UPU_20190725_m1_s00441_db002_w016_lc60"

plot_window_fapar <- function(window_id, root = ROOT) {
  manifest_path <- file.path(
    root, "data", "density_aggregation", "sen2_spataggr_fapar",
    "density_cluster_4000m_aggregate_manifest.csv"
  )
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("aggregation_id", "Delta_Date", "window_crs",
                "cell_xmin", "cell_ymin", "cell_xmax", "cell_ymax")
  if (!all(required %in% names(manifest))) {
    stop("Manifest is missing required columns: ",
         paste(setdiff(required, names(manifest)), collapse = ", "))
  }
  row <- manifest[which(manifest$aggregation_id == window_id), required, drop = FALSE]
  if (nrow(row) != 1L) stop("Expected one manifest row for ", window_id, "; found ", nrow(row))

  target_date <- as.Date(row$Delta_Date)
  if (is.na(target_date)) stop("Invalid Delta_Date for ", window_id)
  year <- as.integer(format(target_date, "%Y"))
  doy <- 1L + 8L * ((as.integer(format(target_date, "%j")) - 1L) %/% 8L)
  composite_start <- as.Date(sprintf("%04d-01-01", year)) + doy - 1L
  composite_end <- min(composite_start + 7L, as.Date(sprintf("%04d-12-31", year)))

  snap_path <- file.path(root, "data", "sentinel2_fapar_composite_snap",
                         paste0(window_id, "_fapar.tif"))
  if (!file.exists(snap_path)) stop("Missing SNAP FAPAR: ", snap_path)
  snap20 <- terra::rast(snap_path)
  if (terra::nlyr(snap20) != 1L) stop("Expected a single-band SNAP FAPAR raster")
  if (!terra::same.crs(snap20, row$window_crs)) stop("SNAP CRS differs from manifest window_crs")
  expected_extent <- as.numeric(unlist(row[c("cell_xmin", "cell_xmax", "cell_ymin", "cell_ymax")]))
  if (any(!is.finite(expected_extent)) ||
      !isTRUE(all.equal(as.vector(terra::ext(snap20)), expected_extent,
                       tolerance = 1e-6, check.attributes = FALSE)) ||
      any(abs(terra::res(snap20) - 20) > 1e-6)) {
    stop("SNAP raster does not match the manifest extent / expected 20 m grid")
  }
  # Both products already represent physical FAPAR. terra applies any file
  # scale/offset metadata on reading: do NOT divide these values again.
  snap20 <- terra::clamp(snap20, lower = 0, upper = 1, values = FALSE)
  names(snap20) <- "FAPAR"

  tiles <- c("h18v03", "h18v04")
  glass_paths <- vapply(tiles, function(tile) {
    folder <- file.path(root, "data", "glass_geotiff", "fapar", tile, year)
    pattern <- sprintf("A%04d%03d\\.%s\\..*\\.tif$", year, doy, tile)
    matches <- list.files(folder, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(matches) != 1L) {
      stop("Expected exactly one GLASS file for ", tile, ", A", year,
           sprintf("%03d", doy), "; found ", length(matches), " in ", folder)
    }
    matches[[1]]
  }, character(1))
  g1 <- terra::rast(glass_paths[[1]])
  g2 <- terra::rast(glass_paths[[2]])
  if (terra::nlyr(g1) != 1L || terra::nlyr(g2) != 1L ||
      !terra::same.crs(g1, g2) || any(abs(terra::res(g1) - terra::res(g2)) > 1e-6) ||
      any(abs(terra::origin(g1) - terra::origin(g2)) > 1e-5)) {
    stop("The two GLASS tiles must be single-band and on the same native grid")
  }
  glass_merged <- terra::merge(g1, g2)
  # A densified raster-extent outline is preferable when transforming the
  # window boundary between UTM and the GLASS projection.
  window <- terra::as.polygons(snap20, extent = TRUE)
  window_glass <- terra::project(window, terra::crs(glass_merged))
  glass250 <- terra::crop(glass_merged, terra::ext(window_glass), snap = "out")
  glass250 <- terra::clamp(glass250, lower = 0, upper = 1, values = FALSE)

  # 250/20 is not an integer, and GLASS nominal 250 m can have a different
  # actual native cell size. Use its EXACT CRS, origin, resolution and extent.
  # GDAL mean uses overlap-weighted contributing non-NA source pixels. This
  # combines the required grid alignment / CRS change and averaging in one step.
  # It is not bilinear interpolation or an integer-factor aggregate().
  if (terra::same.crs(snap20, glass250)) {
    snap250 <- terra::resample(snap20, glass250, method = "mean")
  } else {
    snap250 <- terra::project(snap20, glass250, method = "mean")
  }
  # Keep cells touching the window. Plots below clip their geometry to the
  # original 4 km UTM extent. Edge SNAP means use only available in-window data.
  snap250 <- terra::mask(snap250, window_glass, touches = TRUE)
  glass250 <- terra::mask(glass250, window_glass, touches = TRUE)
  names(snap250) <- names(glass250) <- "FAPAR"

  colours <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(256)
  fill_scale <- function() {
    scale_fill_gradientn(colours = colours, limits = c(0, 1),
                         breaks = seq(0, 1, 0.2), na.value = "white", name = "FAPAR")
  }
  map_theme <- theme_void() + theme(
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
  )
  x_limits <- expected_extent[1:2]
  y_limits <- expected_extent[3:4]
  fine_df <- as.data.frame(snap20, xy = TRUE, na.rm = FALSE)
  p20 <- ggplot(fine_df, aes(x = x, y = y, fill = FAPAR)) +
    geom_raster(interpolate = FALSE) + fill_scale() +
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(title = "SNAP | native 20 m") + map_theme

  coarse_panel <- function(r, title) {
    # Transform cell GEOMETRIES for a common display orientation, not their
    # values. This avoids a second raster resampling that would blur GLASS.
    cells <- terra::as.polygons(r, aggregate = FALSE, round = FALSE,
                               values = TRUE, na.rm = FALSE)
    cells <- sf::st_transform(sf::st_as_sf(cells), crs = sf::st_crs(terra::crs(snap20)))
    ggplot(cells) + geom_sf(aes(fill = FAPAR), colour = NA) + fill_scale() +
      coord_sf(crs = sf::st_crs(cells), datum = NA, xlim = x_limits,
               ylim = y_limits, expand = FALSE) +
      labs(title = title) + map_theme
  }
  p250 <- coarse_panel(snap250, "SNAP | mean on GLASS grid")
  pglass <- coarse_panel(glass250, "GLASS | nominal 250 m")
  grid_size <- paste(sprintf("%.3f", terra::res(glass250)), collapse = " x ")
  units <- if (terra::is.lonlat(glass250)) "degrees" else "projected CRS units"
  caption <- paste0(
    "GLASS composite: ", composite_start, " to ", composite_end,
    " (DOY ", sprintf("%03d", doy), "). Native cell size: ", grid_size, " ", units, ".\n",
    "White = NoData. SNAP means ignore NoData; boundary cells have partial support. ",
    "SNAP source dates were selected within target +/-8 days."
  )
  combined <- (p20 | p250 | pglass) +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste("FAPAR comparison | target", target_date),
                    subtitle = window_id, caption = caption) &
    theme(legend.position = "right")
  message("GLASS files: ", paste(basename(glass_paths), collapse = " + "))
  message("GLASS native resolution: ", grid_size, " ", units)
  message("Middle/right panels share the same GLASS grid; no new files written.")
  print(combined)
  invisible(list(plot = combined, snap20 = snap20, snap_on_glass_grid = snap250,
                 glass250 = glass250, window_id = window_id))
}

comparison <- plot_window_fapar(WINDOW_ID)

# For another window after sourcing:
# comparison <- plot_window_fapar("YOUR_WINDOW_ID")
# Optional save:
# ggplot2::ggsave(file.path(ROOT, "fapar_comparison.png"), comparison$plot,
#                 width = 15, height = 6, dpi = 200)
