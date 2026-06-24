library(tidyverse)
library(terra)

# Inspect a GLASS HDF file and export small data-frame views of its structure.
# Keep this FALSE at first: a full 250 m tile can contain tens of millions of cells.
make_full_df <- FALSE
sample_size <- 10000

hdf_path <- "GLASS09D01.V60.A2024193.h18v04.2026046.hdf"
output_dir <- file.path("data", "glass_hdf_inspection")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

t <- rast(hdf_path)
plot(t)
if (!file.exists(hdf_path)) {
  stop("File not found: ", hdf_path)
}

message("Reading HDF subdatasets from: ", hdf_path)
glass_sds <- sds(hdf_path)

message("\nSubdatasets/layers found:")
print(glass_sds)

subdataset_summary <- map_dfr(seq_along(glass_sds), function(i) {
  r <- rast(glass_sds[i])

  tibble(
    subdataset_id = i,
    name = names(r)[1],
    n_layers = nlyr(r),
    n_rows = nrow(r),
    n_cols = ncol(r),
    n_cells = ncell(r),
    crs = crs(r),
    xmin = xmin(r),
    xmax = xmax(r),
    ymin = ymin(r),
    ymax = ymax(r)
  )
})

message("\nSubdataset summary:")
print(subdataset_summary)

write_csv(
  subdataset_summary,
  file.path(output_dir, "glass_hdf_subdataset_summary.csv")
)

sample_one_subdataset <- function(subdataset_id) {
  r <- rast(glass_sds[subdataset_id])
  names(r) <- make.names(names(r), unique = TRUE)

  sample_df <- spatSample(
    r,
    size = sample_size,
    method = "regular",
    cells = TRUE,
    xy = TRUE,
    na.rm = FALSE,
    as.df = TRUE
  ) %>%
    as_tibble() %>%
    mutate(subdataset_id = subdataset_id, .before = 1)

  sample_df
}

sample_df <- map_dfr(seq_along(glass_sds), sample_one_subdataset)

message("\nSample data-frame structure:")
glimpse(sample_df)

message("\nSample rows:")
print(head(sample_df, 20))

write_csv(
  sample_df,
  file.path(output_dir, "glass_hdf_sample_pixels.csv")
)

if (make_full_df) {
  message("\nBuilding full pixel data frame. This may be large and slow.")

  full_df <- map_dfr(seq_along(glass_sds), function(subdataset_id) {
    r <- rast(glass_sds[subdataset_id])
    names(r) <- make.names(names(r), unique = TRUE)

    as.data.frame(
      r,
      cells = TRUE,
      xy = TRUE,
      na.rm = FALSE
    ) %>%
      as_tibble() %>%
      mutate(subdataset_id = subdataset_id, .before = 1)
  })

  message("\nFull data-frame structure:")
  glimpse(full_df)

  saveRDS(
    full_df,
    file.path(output_dir, "glass_hdf_full_pixels.rds")
  )
}

message("\nWrote inspection outputs to: ", output_dir)
