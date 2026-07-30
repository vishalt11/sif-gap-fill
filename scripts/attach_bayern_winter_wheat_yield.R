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
      )
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
# ww_yield remains in dt/ha; write_csv writes numeric decimals with a dot.
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
  "SIF_",
  months,
  "_calibrated"
)

nirv_columns <- paste0("mean_nirv_", months)
raw_ww10_sif_columns <- paste0("raw_ww10_SIF_", months)

model_columns <- c(
  identifier_columns,
  calibrated_crop_pure_sif_columns,
  "ww_pct",
  nirv_columns,
  raw_ww10_sif_columns,
  "ww_yield"
)

missing_model_columns <- setdiff(model_columns, names(model_data))

if (length(missing_model_columns) > 0) {
  stop(
    "The joined table is missing requested modelling columns: ",
    paste(missing_model_columns, collapse = ", ")
  )
}

rows_without_yield <- sum(is.na(model_data$ww_yield))

model_data <- model_data %>%
  filter(!is.na(ww_yield)) %>%
  select(all_of(model_columns))

message("Rows removed because ww_yield is NA: ", rows_without_yield)
message("Temporary modelling rows: ", nrow(model_data))
message("The temporary reduced model_data object has not been saved.")

print(summary(model_data$ww_yield))

message("Missing values by retained column:")
sort(colSums(is.na(model_data)))

model_data <- model_data %>% select(nuts_id, nuts3, mgrs_tile, year,
                                    SIF_March_calibrated, SIF_April_calibrated, SIF_May_calibrated, SIF_June_calibrated, SIF_July_calibrated,
                                    mean_nirv_March, mean_nirv_April, mean_nirv_May, mean_nirv_June, mean_nirv_July,
                                    raw_ww10_SIF_March, raw_ww10_SIF_May, raw_ww10_SIF_June, raw_ww10_SIF_July,
                                    ww_pct, ww_yield)
colnames(model_data)

#-------------------------------------------------------------------------------
library(tidyverse)

months <- c("March", "April", "May", "June", "July")

identifier_columns <- c(
  "nuts_id",
  "nuts3",
  "mgrs_tile",
  "year"
)

crop_pure_columns <- paste0("SIF_", months, "_calibrated")
nirv_columns <- paste0("mean_nirv_", months)
raw_sif_columns <- paste0("raw_ww10_SIF_", months)

selected_columns <- c(
  identifier_columns,
  crop_pure_columns,
  nirv_columns,
  raw_sif_columns,
  "ww_pct",
  "ww_yield"
)

model_data <- model_data %>%
  select(all_of(selected_columns))

crop_pure_matrix <- model_data %>%
  select(all_of(crop_pure_columns)) %>%
  as.matrix()

raw_sif_matrix <- model_data %>%
  select(all_of(raw_sif_columns)) %>%
  as.matrix()

# Use only months for which both raw and crop-pure SIF are available.
# This gives the two seasonal predictors identical temporal support.
matched_month_mask <- is.finite(raw_sif_matrix) &
  is.finite(crop_pure_matrix)

n_raw_months_available <- rowSums(is.finite(raw_sif_matrix))
n_matched_months <- rowSums(matched_month_mask)

raw_matched_sum <- rowSums(
  ifelse(matched_month_mask, raw_sif_matrix, 0)
)

crop_pure_matched_sum <- rowSums(
  ifelse(matched_month_mask, crop_pure_matrix, 0)
)

raw_seasonal_sif_matched <- ifelse(
  n_matched_months > 0,
  raw_matched_sum / n_matched_months,
  NA_real_
)

crop_pure_seasonal_sif_matched <- ifelse(
  n_matched_months > 0,
  crop_pure_matched_sum / n_matched_months,
  NA_real_
)

matched_months_used <- apply(
  matched_month_mask,
  1,
  function(available) {
    used <- months[available]
    
    if (length(used) == 0) {
      return(NA_character_)
    }
    
    paste(used, collapse = ";")
  }
)

seasonal_model_data <- model_data %>%
  mutate(
    n_raw_months_available = n_raw_months_available,
    n_matched_months = n_matched_months,
    matched_months_used = matched_months_used,
    raw_seasonal_sif_matched = raw_seasonal_sif_matched,
    crop_pure_seasonal_sif_matched =
      crop_pure_seasonal_sif_matched
  ) %>%
  select(
    all_of(identifier_columns),
    raw_seasonal_sif_matched,
    crop_pure_seasonal_sif_matched,
    n_raw_months_available,
    n_matched_months,
    matched_months_used,
    all_of(nirv_columns),
    ww_pct,
    ww_yield
  )

# Primary dataset: at least one matched month.
seasonal_model_data_min1 <- seasonal_model_data %>%
  filter(n_matched_months >= 1)

# Sensitivity dataset: at least two matched months.
seasonal_model_data_min2 <- seasonal_model_data %>%
  filter(n_matched_months >= 2)

message("All region-year rows: ", nrow(seasonal_model_data))
message(
  "Rows with >= 1 matched month: ",
  nrow(seasonal_model_data_min1)
)
message(
  "Rows with >= 2 matched months: ",
  nrow(seasonal_model_data_min2)
)

print(table(seasonal_model_data$n_matched_months, useNA = "ifany"))
print(colSums(is.na(seasonal_model_data_min2)))

#-------------------------------------------------------------------------------

nirv_columns <- paste0(
  "mean_nirv_",
  c("March", "April", "May", "June", "July")
)

train_data <- seasonal_model_data_min1 %>%
  filter(year %in% 2019:2022)

test_data <- seasonal_model_data_min1 %>%
  filter(year == 2024)

# Record which NIRv values were originally missing.
for (column in nirv_columns) {
  missing_column <- paste0(column, "_missing")
  
  train_data[[missing_column]] <- as.integer(
    is.na(train_data[[column]])
  )
  
  test_data[[missing_column]] <- as.integer(
    is.na(test_data[[column]])
  )
}

# Calculate month-specific NIRv means from training data only.
nirv_training_means <- train_data %>%
  summarise(
    across(
      all_of(nirv_columns),
      ~ mean(.x, na.rm = TRUE)
    )
  )

# Apply the training means to both training and test data.
for (column in nirv_columns) {
  imputation_value <- nirv_training_means[[column]]
  
  train_data[[column]][is.na(train_data[[column]])] <-
    imputation_value
  
  test_data[[column]][is.na(test_data[[column]])] <-
    imputation_value
}

print(nirv_training_means)
sort(colSums(is.na(train_data)))
sort(colSums(is.na(test_data)))


final_df <- rbind(train_data, test_data)

write_csv(final_df, 'data/winter_wheat_yield_model/ww_model_data.csv')
