library(tidyverse)
library(mgcv)
library(ranger)
library(xgboost)
library(earth)

set.seed(42)

# model_data_path <- "data/final_model.csv"
# target_col <- "sif"
train_fraction <- 0.8
# 
# model_df <- read_csv(model_data_path, show_col_types = FALSE) %>%
#   mutate(across(everything(), as.numeric)) %>%
#   drop_na()

#model_df <- model_df %>% select(-c(mean_ndvi, mean_ndre, mean_ndre8a, mean_psri, mean_osavi, mean_ndwi, mean_nirv, mean_tcari))

model_df <- read.csv('data/extracted_modis_data/df_3state_fp_par_big.csv')
colnames(model_df)
model_df1 <- read.csv('data/extracted_modis_data/df_3state_evi_big.csv')

common_cols <- intersect(names(model_df), names(model_df1))

model_df <- model_df %>%
  left_join(
    model_df1 %>%
      select(all_of(common_cols), mean_evi),
    by = common_cols
  )

model_df$month <- lubridate::month(as.Date(model_df$Delta_Date))
model_df <- model_df %>% filter(Metadata.MeasurementMode == 0)

model_df <- model_df %>% 
  select(c(Daily_SIF_740nm, Latitude, Longitude, mean_fapar, mean_par, apar, month,mean_evi))



target_col <- "Daily_SIF_740nm"

# if (!target_col %in# if (!target_col %in# if (!target_col %in% names(model_df)) {
#   stop("Missing target column: ", target_col)
# }
# 
# if (nrow(model_df) < 2) {
#   stop("Need at least two complete modeling rows.")
# }
# 
feature_cols <- setdiff(names(model_df), target_col)
# if (length(feature_cols) == 0) {
#   stop("No predictor columns found.")
# }

calc_metrics <- function(test_data, predictions, model_name) {
  residual <- test_data[[target_col]] - predictions
  tss <- sum((test_data[[target_col]] - mean(test_data[[target_col]], na.rm = TRUE))^2, na.rm = TRUE)

  tibble(
    model = model_name,
    n_train = nrow(train_df),
    n_test = nrow(test_data),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    r_squared = 1 - sum(residual^2, na.rm = TRUE) / tss
  )
}

model_df <- model_df %>% drop_na()

train_n <- floor(train_fraction * nrow(model_df))
train_n <- min(max(train_n, 1), nrow(model_df) - 1)
train_rows <- sample(seq_len(nrow(model_df)), size = train_n)

train_df <- model_df[train_rows, ]
test_df <- model_df[-train_rows, ]

model_metrics <- list()

#-------------------------------------------------------------------------------
# Linear regression baseline

sif_lm <- lm(Daily_SIF_740nm ~ ., data = train_df)

lm_predictions <- as.numeric(predict(sif_lm, newdata = test_df))
model_metrics$linear_regression <- calc_metrics(
  test_df,
  lm_predictions,
  "linear_regression"
)

print(model_metrics$linear_regression)
summary(sif_lm)


#-------------------------------------------------------------------------------
# GAM

make_gam_smooth <- function(col_name, data, default_k = 10) {
  n_unique <- n_distinct(data[[col_name]], na.rm = TRUE)

  if (n_unique < 3) {
    return(NULL)
  }

  smooth_k <- min(default_k, n_unique)
  sprintf("s(%s, k = %d)", col_name, smooth_k)
}

spatial_cols <- c("Latitude", "Longitude")
non_spatial_cols <- setdiff(feature_cols, spatial_cols)

gam_terms <- map(
  non_spatial_cols,
  ~ make_gam_smooth(
    .x,
    train_df,
    default_k = if_else(.x == "month", 6L, 10L)
  )
) %>%
  compact() %>%
  unlist(use.names = FALSE)

if (all(spatial_cols %in% feature_cols) &&
    n_distinct(train_df$Latitude, na.rm = TRUE) >= 3 &&
    n_distinct(train_df$Longitude, na.rm = TRUE) >= 3) {
  gam_terms <- c(gam_terms, "te(Latitude, Longitude, k = c(10, 10))")
} else {
  gam_terms <- c(
    gam_terms,
    map(intersect(spatial_cols, feature_cols), ~ make_gam_smooth(.x, train_df)) %>%
      compact() %>%
      unlist(use.names = FALSE)
  )
}

if (length(gam_terms) == 0) {
  stop("No GAM terms could be created from the predictor columns.")
}

gam_formula <- as.formula(paste(target_col, "~", paste(gam_terms, collapse = " + ")))

sif_gam <- bam(
  gam_formula,
  data = train_df,
  method = "fREML",
  select = TRUE,
  discrete = TRUE
)

gam_predictions <- as.numeric(predict(sif_gam, newdata = test_df))
model_metrics$gam <- calc_metrics(
  test_df,
  gam_predictions,
  "gam"
)

print(model_metrics$gam)
summary(sif_gam)


#-------------------------------------------------------------------------------
# Random forest

rf_mtry <- max(1, floor(sqrt(length(feature_cols))))

sif_rf <- ranger(
  Daily_SIF_740nm ~ .,
  data = train_df,
  num.trees = 500,
  mtry = rf_mtry,
  min.node.size = 5,
  importance = "permutation",
  seed = 42
)

rf_predictions <- as.numeric(predict(sif_rf, data = test_df)$predictions)
model_metrics$random_forest <- calc_metrics(
  test_df,
  rf_predictions,
  "random_forest"
)

print(model_metrics$random_forest)
print(sort(sif_rf$variable.importance, decreasing = TRUE))


#-------------------------------------------------------------------------------
# XGBoost


xgb_formula <- as.formula(paste("~", paste(feature_cols, collapse = " + "), "- 1"))
xgb_train_matrix <- model.matrix(xgb_formula, data = train_df)
xgb_test_matrix <- model.matrix(xgb_formula, data = test_df)

xgb_train <- xgb.DMatrix(
  data = xgb_train_matrix,
  label = train_df[[target_col]]
)

xgb_params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.05,
  max_depth = 4,
  min_child_weight = 5,
  subsample = 0.8,
  colsample_bytree = 0.8
)

sif_xgb <- xgb.train(
  params = xgb_params,
  data = xgb_train,
  nrounds = 300,
  verbose = 0
)

xgb_predictions <- as.numeric(predict(sif_xgb, newdata = xgb_test_matrix))
model_metrics$xgboost <- calc_metrics(
  test_df,
  xgb_predictions,
  "xgboost"
)

print(model_metrics$xgboost)
print(xgb.importance(model = sif_xgb))


#-------------------------------------------------------------------------------
# MARS

sif_mars <- earth(
  sif ~ .,
  data = train_df,
  degree = 2
)

mars_predictions <- as.numeric(predict(sif_mars, newdata = test_df))
model_metrics$mars <- calc_metrics(
  test_df,
  mars_predictions,
  "mars"
)

print(model_metrics$mars)
summary(sif_mars)


#-------------------------------------------------------------------------------
# Model comparison

model_comparison <- bind_rows(model_metrics) %>%
  arrange(rmse)

print(model_comparison)
