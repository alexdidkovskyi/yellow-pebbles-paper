#' @title XGBoost Regression Modeling for Tracer Movement
#' @description Functions for analyzing and predicting tracer displacement and virtual velocity
#' @author Yellow Pebbles Research Team

#' @title Analyze Dataset with XGBoost Regression
#' @description A single function to analyze either displacement or velocity with XGBoost
#' @param reg_data Data frame containing features and target
#' @param target_col Name of the target column (e.g., "graph_dist" or "graph_vel")
#' @param feature_cols Vector of feature column names
#' @param n_iterations Number of iterations for cross-validation
#' @param nrounds Maximum number of boosting rounds
#' @param nfolds Number of cross-validation folds
#' @param seed Random seed for reproducibility
#' @return List with model, metrics, and predictions
#' @export
analyze_regression <- function(
  reg_data,
  target_col,
  feature_cols,
  n_iterations = 1,
  nrounds = 500,
  nfolds = 5,
  seed = 42
) {
  # Set seed for reproducibility
  set.seed(seed)

  # Prepare data
  data_features <- reg_data %>% dplyr::select(all_of(feature_cols))
  target_values <- reg_data[[target_col]]

  # Check if dataset is empty or has too few rows
  if (nrow(reg_data) < 10) {
    warning(paste0(
      "Dataset for ",
      target_col,
      " analysis has only ",
      nrow(reg_data),
      " rows."
    ))
    return(list(
      model = NULL,
      metrics = list(rmse = NA, mae = NA, wasserstein = NA),
      predictions = NULL,
      actual_values = target_values,
      n_observations = length(target_values)
    ))
  }

  # Baseline performance
  baseline_rmse <- MLmetrics::RMSE(target_values, mean(target_values))
  baseline_mae <- MLmetrics::MAE(target_values, mean(target_values))

  cat("Baseline RMSE (mean prediction):", baseline_rmse, "\n")
  cat("Baseline MAE (mean prediction):", baseline_mae, "\n")

  # Storage for metrics
  vec_RMSE <- numeric(n_iterations)
  vec_MAE <- numeric(n_iterations)
  vec_Wass <- numeric(n_iterations)

  # Store last predictions
  last_predictions <- NULL
  best_iteration <- NULL

  # XGBoost parameters
  params <- list(
    max_depth = 6,
    eta = 0.03,
    eval_metric = "rmse",
    objective = "reg:linear",
    nthread = 2
  )

  # Iterate cross-validation
  for (i in 1:n_iterations) {
    # Create DMatrix
    dtrain <- xgboost::xgb.DMatrix(
      data = data_features %>% data.matrix(),
      label = target_values
    )

    # Run cross-validation
    xgb_cv_result <- xgboost::xgb.cv(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      nfold = nfolds,
      early_stopping_rounds = 60,
      print_every_n = 50,
      prediction = TRUE
    )

    # Calculate performance metrics
    pred <- xgb_cv_result$pred
    pred[pred < 0] <- 0
    vec_RMSE[i] <- MLmetrics::RMSE(target_values, pred)
    vec_MAE[i] <- MLmetrics::MAE(target_values, pred)

    tryCatch(
      {
        vec_Wass[i] <- transport::wasserstein1d(target_values, pred, p = 2)
      },
      error = function(e) {
        vec_Wass[i] <- NA
        warning("Wasserstein distance calculation failed: ", e$message)
      }
    )

    # Save predictions from last iteration
    if (i == n_iterations) {
      last_predictions <- pred
      best_iteration <- xgb_cv_result$best_iteration
    }
  }

  # Train a final model for feature importance
  final_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = ifelse(is.null(best_iteration), 110, best_iteration),
    verbose = 0
  )

  # Return results
  results <- list(
    model = final_model,
    metrics = list(
      rmse = mean(vec_RMSE, na.rm = TRUE),
      mae = mean(vec_MAE, na.rm = TRUE),
      wasserstein = mean(vec_Wass, na.rm = TRUE)
    ),
    baseline = list(rmse = baseline_rmse, mae = baseline_mae),
    predictions = last_predictions,
    actual_values = target_values,
    n_observations = length(target_values)
  )

  return(results)
}


# Function for training and evaluating classification models
analyze_classification <- function(
  data,
  target_col,
  feature_cols,
  nfolds = 5,
  nrounds = 500,
  seed = 42
) {
  # Set seed for reproducibility
  set.seed(seed)

  # Prepare data
  data_features <- data %>% dplyr::select(all_of(feature_cols))
  target_values <- data[[target_col]]

  # Check if dataset is empty or has too few rows
  if (nrow(data) < 10) {
    warning(paste0(
      "Dataset for classification has only ",
      nrow(data),
      " rows."
    ))
    return(list(
      model = NULL,
      metrics = list(auc = NA, accuracy = NA, kappa = NA),
      probabilities = NULL,
      actual_values = target_values,
      n_observations = length(target_values)
    ))
  }

  # Baseline performance - predict the majority class
  majority_class <- as.numeric(mean(target_values) >= 0.5)
  baseline_accuracy <- mean(target_values == majority_class)

  cat("Baseline accuracy (majority class):", baseline_accuracy, "\n")

  # Storage for metrics
  accuracy_values <- numeric(nfolds)
  auc_values <- numeric(nfolds)
  kappa_values <- numeric(nfolds)

  # Store predictions
  all_probabilities <- numeric(length(target_values))

  # XGBoost parameters for binary classification
  params <- list(
    max_depth = 6,
    eta = 0.03,
    eval_metric = "auc",
    objective = "binary:logistic",
    nthread = 2
  )

  # Run cross-validation
  xgb_cv_result <- xgboost::xgb.cv(
    params = params,
    data = xgboost::xgb.DMatrix(
      data = data_features %>% as.matrix(),
      label = target_values
    ),
    nrounds = nrounds,
    nfold = nfolds,
    early_stopping_rounds = 60,
    print_every_n = 50,
    prediction = TRUE
  )

  # Calculate metrics
  probabilities <- xgb_cv_result$pred
  all_probabilities <- probabilities

  # Calculate performance metrics
  auc <- MLmetrics::AUC(probabilities, target_values)

  # Convert probabilities to class predictions using 0.5 threshold
  predicted_classes <- as.numeric(probabilities >= 0.5)
  accuracy <- mean(predicted_classes == target_values)
  confusion_matrix <- caret::confusionMatrix(
    factor(predicted_classes, levels = c(0, 1)),
    factor(target_values, levels = c(0, 1))
  )
  kappa <- confusion_matrix$overall["Kappa"]

  # Train a final model for feature importance
  final_model <- xgboost::xgb.train(
    params = params,
    data = xgboost::xgb.DMatrix(
      data = data_features %>% as.matrix(),
      label = target_values
    ),
    nrounds = xgb_cv_result$best_iteration,
    verbose = 0
  )

  # Return results
  results <- list(
    model = final_model,
    metrics = list(
      auc = auc,
      accuracy = accuracy,
      kappa = kappa
    ),
    baseline = list(accuracy = baseline_accuracy),
    probabilities = all_probabilities,
    actual_values = target_values,
    n_observations = length(target_values),
    confusion_matrix = confusion_matrix
  )

  return(results)
}
