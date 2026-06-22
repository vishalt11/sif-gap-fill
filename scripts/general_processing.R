library(tidyverse)
library(sf)
library(terra)


df <- readRDS('data/ba_sif_mgrs_crop_composition.rds')

head(df,1)

d1 <- 'data/spectral_indices_means_nonveg_masked'
d2 <- 'data/band_means_nonveg_masked'

output_dir <- 'data'
output_rds <- file.path(output_dir, 'sif_wasp_features_nonveg_masked.rds')
output_csv <- file.path(output_dir, 'sif_wasp_features_nonveg_masked.csv')

read_csv_dir <- function(dir_path, source_col) {
  csv_files <- list.files(
    dir_path,
    pattern = '\\.csv$',
    full.names = TRUE
  )

  if (length(csv_files) == 0) {
    stop('No CSV files found in: ', dir_path)
  }

  csv_files %>%
    set_names(basename(.)) %>%
    map_dfr(
      ~ read_csv(
        .x,
        col_types = cols(sif_id = col_character(), .default = col_guess()),
        show_col_types = FALSE
      ),
      .id = source_col
    )
}

spectral_indices <- read_csv_dir(d1, 'spectral_indices_csv')
band_means <- read_csv_dir(d2, 'band_means_csv')

join_cols <- c('Daily_SIF_740nm', 'Delta_Time', 'Latitude', 'Longitude')

if (anyDuplicated(spectral_indices[join_cols]) > 0) {
  message('Duplicate join key values found after row-binding CSV files in: ', d1)
}

if (anyDuplicated(band_means[join_cols]) > 0) {
  message('Duplicate join key values found after row-binding CSV files in: ', d2)
}

band_only_cols <- setdiff(names(band_means), names(spectral_indices))

sif_wasp_features_nonveg_masked <- spectral_indices %>%
  full_join(
    band_means %>% select(all_of(join_cols), all_of(band_only_cols)),
    by = join_cols
  )

colnames(sif_wasp_features_nonveg_masked)

sif_wasp_features_nonveg_masked <- sif_wasp_features_nonveg_masked %>% 
  select(-c('contaminated_20m_fraction_threshold', 'spectral_indices_csv', 'source_file', 'file_id',
            'sif_id', 'band_means_csv', 'wasp_product_id', 'non_crop_low_ndvi_threshold', 'sif_extract_id'))

saveRDS(sif_wasp_features_nonveg_masked, output_rds)
write_csv(sif_wasp_features_nonveg_masked, output_csv)


