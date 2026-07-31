library(tidyverse)

# Join Bavaria winter-wheat yield to the prepared NUTS3 predictor table.
# The names in predictors$nuts3 are expected to match the yield-file names.

predictor_csv <- paste0(
  "data/winter_wheat_yield_model/raw_mixed_sif/",
  "nuts3_crop_pure_and_raw_mixed_sif_predictors.csv"
)

yield_dir <- "data/bayern_yield"

output_csv <- paste0(
  "data/winter_wheat_yield_model/raw_mixed_sif/",
  "nuts3_crop_pure_and_raw_mixed_sif_predictors_with_yield.csv"
)

yield_column <- "Winterweizen (einschl. Dinkel und Einkorn)"
source_encoding <- "Windows-1252"

yield_files <- list.files(
  yield_dir,
  pattern = "^[0-9]{4}\\.csv$",
  full.names = TRUE
)

if (length(yield_files) == 0) {
  stop("No yearly yield CSV files found in: ", yield_dir)
}

read_yield_file <- function(path) {
  year <- as.integer(tools::file_path_sans_ext(basename(path)))

  yield_df <- readr::read_delim(
    path,
    delim = ";",
    col_types = cols(.default = col_character()),
    locale = locale(encoding = source_encoding),
    trim_ws = TRUE,
    show_col_types = FALSE
  )

  required_columns <- c("name", yield_column)
  missing_columns <- setdiff(required_columns, names(yield_df))

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in ",
      path,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  yield_df %>%
    transmute(
      year = year,
      name = stringr::str_squish(name),
      ww_yield = readr::parse_double(
        .data[[yield_column]],
        na = c("", "/", "-", "NA"),
        locale = locale(decimal_mark = ",")
      ) / 10
    )
}

# One table containing every available Bavaria yield year.
bayern_yield <- purrr::map_dfr(yield_files, read_yield_file)

duplicate_yield_keys <- bayern_yield %>%
  count(name, year) %>%
  filter(n > 1)

if (nrow(duplicate_yield_keys) > 0) {
  print(duplicate_yield_keys, n = Inf)
  stop("The combined yield data contain duplicate name/year rows.")
}

predictors <- readr::read_csv(
  predictor_csv,
  locale = locale(encoding = source_encoding),
  show_col_types = FALSE
) %>%
  mutate(
    nuts3 = stringr::str_squish(nuts3),
    year = as.integer(year)
  )

required_predictor_columns <- c("nuts3", "year")
missing_predictor_columns <- setdiff(
  required_predictor_columns,
  names(predictors)
)

if (length(missing_predictor_columns) > 0) {
  stop(
    "Predictor table is missing columns: ",
    paste(missing_predictor_columns, collapse = ", ")
  )
}

if ("ww_yield" %in% names(predictors)) {
  predictors <- predictors %>% select(-ww_yield)
}

model_data <- predictors %>%
  left_join(
    bayern_yield,
    by = c("nuts3" = "name", "year")
  )

if (nrow(model_data) != nrow(predictors)) {
  stop("The yield join unexpectedly changed the predictor row count.")
}

unmatched_rows <- model_data %>%
  filter(is.na(ww_yield)) %>%
  distinct(nuts3, year) %>%
  arrange(year, nuts3)

if (nrow(unmatched_rows) > 0) {
  message(
    "Unmatched or suppressed NUTS3/year combinations: ",
    nrow(unmatched_rows)
  )
  print(unmatched_rows, n = Inf)
}

# Save the complete joined table before making the temporary modelling subset.
# The source values are dt/ha; ww_yield is stored here in metric tons/ha.
readr::write_csv(model_data, output_csv, na = "")

message("Combined yearly yield rows: ", nrow(bayern_yield))
message("Input predictor rows: ", nrow(predictors))
message("Rows with ww_yield: ", sum(!is.na(model_data$ww_yield)))
message("Rows without ww_yield: ", sum(is.na(model_data$ww_yield)))
message("Saved complete joined table: ", output_csv)

#-------------------------------------------------------------------------------

months <- c("March", "April", "May", "June", "July")

identifier_columns <- c(
  "nuts_id",
  "nuts3",
  "mgrs_tile",
  "year"
)

calibrated_crop_pure_sif_columns <- paste0(
  "SIF_", months, "_calibrated"
)
nirv_columns <- paste0("mean_nirv_", months)
raw_all_sif_columns <- paste0("raw_all_SIF_", months)
raw_ww05_sif_columns <- paste0("raw_ww05_SIF_", months)
raw_ww10_sif_columns <- paste0("raw_ww10_SIF_", months)
raw_ww30_sif_columns <- paste0("raw_ww30_SIF_", months)

model_columns <- c(
  identifier_columns,
  calibrated_crop_pure_sif_columns,
  nirv_columns,
  raw_all_sif_columns,
  raw_ww05_sif_columns,
  raw_ww10_sif_columns,
  raw_ww30_sif_columns,
  "ww_pct",
  "ww_yield"
)

missing_model_columns <- setdiff(model_columns, names(model_data))

if (length(missing_model_columns) > 0) {
  stop(
    "The joined table is missing requested modelling columns: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

# Keep all rows and all missing values for inspection before choosing a
# missing-data strategy.
model_data <- model_data %>%
  select(all_of(model_columns)) %>%
  arrange(year, nuts_id)

message("Reduced modelling rows retained for NA inspection: ", nrow(model_data))

print(summary(model_data$ww_yield))

message("Missing values by retained column:")

model_data <- model_data %>% filter(!is.na(ww_yield))

print(sort(colSums(is.na(model_data))))

#-------------------------------------------------------------------------------
# Seasonal raw-SIF summaries and training-only mean imputation

row_mean_observed_months <- function(data, columns) {
  values <- data %>%
    select(all_of(columns)) %>%
    as.matrix()

  available <- is.finite(values)
  n_available <- rowSums(available)
  value_sum <- rowSums(ifelse(available, values, 0))

  ifelse(n_available > 0, value_sum / n_available, NA_real_)
}

# Calculate these before imputation so each seasonal value represents only the
# raw OCO-2 months genuinely observed for that region-year.
model_data$raw_all_seasonal_sif <- row_mean_observed_months(
  model_data,
  raw_all_sif_columns
)
model_data$raw_ww05_seasonal_sif <- row_mean_observed_months(
  model_data,
  raw_ww05_sif_columns
)
model_data$raw_ww10_seasonal_sif <- row_mean_observed_months(
  model_data,
  raw_ww10_sif_columns
)
model_data$raw_ww30_seasonal_sif <- row_mean_observed_months(
  model_data,
  raw_ww30_sif_columns
)

raw_seasonal_columns <- c(
  "raw_all_seasonal_sif",
  "raw_ww05_seasonal_sif",
  "raw_ww10_seasonal_sif",
  "raw_ww30_seasonal_sif"
)

imputation_columns <- c(
  calibrated_crop_pure_sif_columns,
  nirv_columns,
  raw_all_sif_columns,
  raw_ww05_sif_columns,
  raw_ww10_sif_columns,
  raw_ww30_sif_columns,
  raw_seasonal_columns,
  "ww_pct"
)

# Estimate imputation values from training years only. These same fixed means
# are applied to 2024 without using any 2024 values to estimate them.
training_imputation_means <- model_data %>%
  filter(year %in% 2019:2022) %>%
  summarise(
    across(
      all_of(imputation_columns),
      ~ mean(.x, na.rm = TRUE)
    )
  )

invalid_imputation_columns <- imputation_columns[
  !is.finite(
    unlist(
      training_imputation_means[1, imputation_columns],
      use.names = FALSE
    )
  )
]

if (length(invalid_imputation_columns) > 0) {
  stop(
    "No finite 2019-2022 mean is available for: ",
    paste(invalid_imputation_columns, collapse = ", ")
  )
}

missing_before_imputation <- sort(colSums(is.na(model_data)))

for (column in imputation_columns) {
  missing <- is.na(model_data[[column]])
  model_data[[column]][missing] <- training_imputation_means[[column]]
}

message("Training-only means used for predictor imputation:")
print(training_imputation_means)

message("Missing values before imputation:")
print(missing_before_imputation)

message("Missing values after predictor imputation:")
print(sort(colSums(is.na(model_data))))


write_csv(model_data, 'data/winter_wheat_yield_model/ww_model_data.csv')





