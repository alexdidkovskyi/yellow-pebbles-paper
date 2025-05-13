# Complete Analysis of River Bed Load Transport
# This script demonstrates the full analysis workflow

system("ulimit -c unlimited")

# Load dependencies
source("src/01_initial_data_processing.R")
source("src/02_markov_chains.R")
source("src/03_xgboost_models.R")
source("src/04_validation.R")
source("src/utils/plotting.R")


# Set seed for reproducibility
set.seed(1234)

# -------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# -------------------------------------------------------------------------

# Load tracer data
cat("Loading data...\n")
list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())


# -------------------------------------------------------------------------
# 2. Markov Chain Analysis
# -------------------------------------------------------------------------

cat("Performing Markov Chain analysis...\n")

# Find optimal thresholds using Jenks Natural Breaks
thresholds <- find_optimal_thresholds(df_models_, 
                                      displacement_col = "graph_dist",
                                      velocity_col = "graph_vel")

displacement_threshold <- thresholds$displacement_threshold
velocity_threshold <- thresholds$velocity_threshold

cat("Displacement threshold:", displacement_threshold, "m \n")
cat("Velocity threshold:", velocity_threshold, "m/s \n")

outlier_counts <- df_models_ %>%
  filter(graph_dist > 0) %>%
  count(
    dist_outlier = graph_dist >= displacement_threshold,
    vel_outlier = graph_vel >= velocity_threshold
  ) %>%
  mutate(
    dist_outlier = ifelse(dist_outlier, "Long Displacement", "Typical displacement"),
    vel_outlier = ifelse(vel_outlier, "High Velocity", "Typical velocity")
  ) %>%
  rename(Count = n)

print(outlier_counts)


displacement_density_plot <- ggplot(df_models_ %>% filter(graph_dist > 0), aes(x = graph_dist)) +
  geom_density() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_vline(xintercept = displacement_threshold, linetype = "dashed") +
  labs(x = "Displacement [m]", y = "Density") +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggsave("data/results/figures/displacement_density_plot.png", displacement_density_plot, width = 8, height = 6)


velocity_density_plot <- ggplot(df_models_ %>% filter(graph_vel > 0), aes(x = graph_vel)) +
  geom_density() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_vline(xintercept = velocity_threshold, linetype = "dashed") +
  labs(x = "Virtual velocity [m/h]", y = "Density") +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggsave("data/results/figures/velocity_density_plot.png", velocity_density_plot, width = 8, height = 6)


# Create Markov Chain for binary movement
mc_two_state <- create_two_state_mc(df_models_, 
                                    displacement_col = "graph_dist",
                                    exclude_pairs = c("9-10", "18-19"))

# Create Markov Chains for displacement and velocity
mc_displacement <- create_displacement_mc(df_models_,
                                          displacement_col = "graph_dist",
                                          long_threshold = displacement_threshold,
                                          exclude_pairs = c("9-10", "18-19"))

mc_velocity <- create_velocity_mc(df_models_,
                                  velocity_col = "graph_vel",
                                  high_threshold = velocity_threshold,
                                  exclude_pairs = c("9-10", "18-19"))

# Print transition matrices
cat("\nTwo-state movement Markov Chain:\n")
print(mc_two_state)

cat("\nDisplacement Markov Chain:\n")
print(mc_displacement)

cat("\nVelocity Markov Chain:\n")
print(mc_velocity)

# -------------------------------------------------------------------------
# 3. Effect of Flood History
# -------------------------------------------------------------------------

cat("Analyzing effect of flood history...\n")

# Run the analysis
results <- compare_consecutive_events(df_models_)

# Display and summarize results
summary <- summarize_comparison_results(results)

# Create a publication-ready table
pub_table <- create_publication_table(results)
print(pub_table)


# Save results to Excel
save_comparison_results(results, "data/results/movement_testing_by_events_pairs.xlsx")

# Save publication table to Excel
openxlsx::write.xlsx(pub_table, "data/results/movement_testing_by_events_pairs_publication_table.xlsx")

# -------------------------------------------------------------------------
# 4. XGBoost Modeling
# -------------------------------------------------------------------------

cat("Training XGBoost models...\n")

# Define feature columns for modeling
feature_cols <- c(
  # Event characteristics
  "duration", "mean_h", "mean_Q",
  
  # Location
  "X_Start", "Y_Start",
  
  # Morphological units
  "TipoStartBanks", "TipoStartBars_sedimentbuildupzones", "TipoStartCascade",
  "TipoStartPlanebed", "TipoStartPools", "TipoStartRun_rapid"
)
feature_cols_with_dimensions <- c(feature_cols, c("pebble_PC1", "pebble_PC2"))

# Create separate datasets for typical and extreme values
typical_displacement_data <- df_models_[df_models_$graph_dist <= displacement_threshold, ] %>% filter(graph_vel > 0)
extreme_displacement_data <- df_models_[df_models_$graph_dist > displacement_threshold, ]
typical_velocity_data <- df_models_[df_models_$graph_vel <= velocity_threshold, ] %>% filter(graph_vel > 0)
extreme_velocity_data <- df_models_[df_models_$graph_vel > velocity_threshold, ] 

# Train typical displacement model
cat("Training typical displacement model...\n")
typical_displacement_results <- analyze_regression(
  reg_data = typical_displacement_data,
  target_col = "graph_dist",
  feature_cols = feature_cols
)

# Train extreme displacement model
cat("Training extreme displacement model...\n")
extreme_displacement_results <- analyze_regression(
  reg_data = extreme_displacement_data,
  target_col = "graph_dist", 
  feature_cols = feature_cols
)

# Train typical velocity model
cat("Training typical velocity model...\n")
typical_velocity_results <- analyze_regression(
  reg_data = typical_velocity_data,
  target_col = "graph_vel",
  feature_cols = feature_cols
)

# Train extreme velocity model
cat("Training extreme velocity model...\n")
extreme_velocity_results <- analyze_regression(
  reg_data = extreme_velocity_data,
  target_col = "graph_vel",
  feature_cols = feature_cols
)

# Train typical velocity & displacement model
cat("Training typical velocity & displacement model...\n")
typical_velocity_displacement_results <- analyze_regression(
  reg_data = typical_displacement_data %>% filter(graph_vel <= velocity_threshold),
  target_col = "graph_vel",
  feature_cols = feature_cols
)

# Train typical displacement model with tracer dimensions
cat("Training typical displacement model with tracer dimensions...\n")
typical_displacement_with_dims_results <- analyze_regression(
  reg_data = typical_displacement_data,
  target_col = "graph_dist",
  feature_cols = feature_cols_with_dimensions
)

# Train extreme displacement model with tracer dimensions
cat("Training extreme displacement model with tracer dimensions...\n")
extreme_displacement_with_dims_results <- analyze_regression(
  reg_data = extreme_displacement_data,
  target_col = "graph_dist", 
  feature_cols = feature_cols_with_dimensions
)

# Train typical velocity model with tracer dimensions
cat("Training typical velocity model with tracer dimensions...\n")
typical_velocity_with_dims_results <- analyze_regression(
  reg_data = typical_velocity_data,
  target_col = "graph_vel",
  feature_cols = feature_cols_with_dimensions
)

# Train extreme velocity model with tracer dimensions
cat("Training extreme velocity model with tracer dimensions...\n")
extreme_velocity_with_dims_results <- analyze_regression(
  reg_data = extreme_velocity_data,
  target_col = "graph_vel",
  feature_cols = feature_cols_with_dimensions
)


# Combine results
model_results <- list(
  typical_displacement = typical_displacement_results,
  extreme_displacement = extreme_displacement_results,
  typical_velocity = typical_velocity_results,
  extreme_velocity = extreme_velocity_results,
  typical_displacement_with_dims = typical_displacement_with_dims_results,
  extreme_displacement_with_dims = extreme_displacement_with_dims_results,
  typical_velocity_with_dims = typical_velocity_with_dims_results,
  extreme_velocity_with_dims = extreme_velocity_with_dims_results
)

# Display model results with baseline comparisons
cat("\nModel results:\n")

# Typical displacement
cat("Typical displacement model:\n")
cat("  RMSE =", model_results$typical_displacement$metrics$rmse, 
    "(Baseline:", model_results$typical_displacement$baseline$rmse, ")\n")
cat("  MAE =", model_results$typical_displacement$metrics$mae, 
    "(Baseline:", model_results$typical_displacement$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$typical_displacement$metrics$wasserstein, "\n\n")

# Typical displacement with dimensions
cat("Typical displacement model with tracer dimensions:\n")
cat("  RMSE =", model_results$typical_displacement_with_dims$metrics$rmse, 
    "(Baseline:", model_results$typical_displacement_with_dims$baseline$rmse, ")\n")
cat("  MAE =", model_results$typical_displacement_with_dims$metrics$mae, 
    "(Baseline:", model_results$typical_displacement_with_dims$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$typical_displacement_with_dims$metrics$wasserstein, "\n\n")


# Extreme displacement
cat("Extreme displacement model:\n")
cat("  RMSE =", model_results$extreme_displacement$metrics$rmse, 
    "(Baseline:", model_results$extreme_displacement$baseline$rmse, ")\n")
cat("  MAE =", model_results$extreme_displacement$metrics$mae, 
    "(Baseline:", model_results$extreme_displacement$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$extreme_displacement$metrics$wasserstein, "\n\n")


# Extreme displacement with dimensions
cat("Extreme displacement model with tracer dimensions:\n")
cat("  RMSE =", model_results$extreme_displacement_with_dims$metrics$rmse, 
    "(Baseline:", model_results$extreme_displacement_with_dims$baseline$rmse, ")\n")
cat("  MAE =", model_results$extreme_displacement_with_dims$metrics$mae, 
    "(Baseline:", model_results$extreme_displacement_with_dims$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$extreme_displacement_with_dims$metrics$wasserstein, "\n\n")


# Typical velocity
cat("Typical velocity model:\n")
cat("  RMSE =", model_results$typical_velocity$metrics$rmse, 
    "(Baseline:", model_results$typical_velocity$baseline$rmse, ")\n")
cat("  MAE =", model_results$typical_velocity$metrics$mae, 
    "(Baseline:", model_results$typical_velocity$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$typical_velocity$metrics$wasserstein, "\n\n")


# Typical velocity with dimensions
cat("Typical velocity model with tracer dimensions:\n")
cat("  RMSE =", model_results$typical_velocity_with_dims$metrics$rmse, 
    "(Baseline:", model_results$typical_velocity_with_dims$baseline$rmse, ")\n")
cat("  MAE =", model_results$typical_velocity_with_dims$metrics$mae, 
    "(Baseline:", model_results$typical_velocity_with_dims$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$typical_velocity_with_dims$metrics$wasserstein, "\n\n")


# Extreme velocity
cat("Extreme velocity model:\n")
cat("  RMSE =", model_results$extreme_velocity$metrics$rmse, 
    "(Baseline:", model_results$extreme_velocity$baseline$rmse, ")\n")
cat("  MAE =", model_results$extreme_velocity$metrics$mae, 
    "(Baseline:", model_results$extreme_velocity$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$extreme_velocity$metrics$wasserstein, "\n")


# Extreme velocity with dimensions
cat("Extreme velocity model with tracer dimensions:\n")
cat("  RMSE =", model_results$extreme_velocity_with_dims$metrics$rmse, 
    "(Baseline:", model_results$extreme_velocity_with_dims$baseline$rmse, ")\n")
cat("  MAE =", model_results$extreme_velocity_with_dims$metrics$mae, 
    "(Baseline:", model_results$extreme_velocity_with_dims$baseline$mae, ")\n")
cat("  Wasserstein =", model_results$extreme_velocity_with_dims$metrics$wasserstein, "\n")

# Get feature importance
typical_disp_importance <- xgb.importance(model = model_results$typical_displacement$model)
typical_vel_importance <- xgb.importance(model = model_results$typical_velocity$model)
extreme_disp_importance <- xgb.importance(model = model_results$extreme_displacement$model)
extreme_vel_importance <- xgb.importance(model = model_results$extreme_velocity$model)

# Print top 5 features for each model
cat("\nTop 5 features for typical displacement model:\n")
print(head(typical_disp_importance[order(-typical_disp_importance$Gain), c("Feature", "Gain")], 5))

cat("\nTop 5 features for typical velocity model:\n")
print(head(typical_vel_importance[order(-typical_vel_importance$Gain), c("Feature", "Gain")], 5))

cat("\nTop 5 features for extreme displacement model:\n")
print(head(extreme_disp_importance[order(-extreme_disp_importance$Gain), c("Feature", "Gain")], 5))

cat("\nTop 5 features for extreme velocity model:\n")
print(head(extreme_vel_importance[order(-extreme_vel_importance$Gain), c("Feature", "Gain")], 5))

# Create visualizations of actual vs predicted distributions
# For typical displacement
typical_disp_data <- data.frame(
  Type = rep(c("Actual", "Predicted"), 
             c(length(typical_displacement_results$actual_values), 
               length(typical_displacement_results$predictions))),
  Value = c(typical_displacement_results$actual_values, 
            typical_displacement_results$predictions)
)

typical_disp_plot <- ggplot(typical_disp_data, aes(x = Value, fill = Type)) + 
  geom_density(size = 1.2) +
  theme_bw() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  labs(title = "Typical Displacement: Actual vs Predicted", 
       x = "Displacement [m]", 
       y = "Density")

# For typical velocity
typical_vel_data <- data.frame(
  Type = rep(c("Actual", "Predicted"), 
             c(length(typical_velocity_results$actual_values), 
               length(typical_velocity_results$predictions))),
  Value = c(typical_velocity_results$actual_values, 
            typical_velocity_results$predictions)
)

typical_vel_plot <- ggplot(typical_vel_data, aes(x = Value, fill = Type)) + 
  geom_density(size = 1.2) +
  theme_bw() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  labs(title = "Typical Velocity: Actual vs Predicted", 
       x = "Velocity [m/h]", 
       y = "Density")

# Save plots
ggsave("data/results/figures/typical_displacement_comparison.png", plot = typical_disp_plot)
ggsave("data/results/figures/typical_velocity_comparison.png", plot = typical_vel_plot)

# Create visualizations for extreme values
# For extreme displacement
extreme_disp_data <- data.frame(
  Type = rep(c("Actual", "Predicted"), 
             c(length(extreme_displacement_results$actual_values), 
               length(extreme_displacement_results$predictions))),
  Value = c(extreme_displacement_results$actual_values, 
            extreme_displacement_results$predictions)
)

extreme_disp_plot <- ggplot(extreme_disp_data, aes(x = Value, fill = Type)) + 
  geom_density(size = 1.2) +
  theme_bw() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  labs(title = "Extreme Displacement: Actual vs Predicted", 
       x = "Displacement [m]", 
       y = "Density")

# For extreme velocity
extreme_vel_data <- data.frame(
  Type = rep(c("Actual", "Predicted"), 
             c(length(extreme_velocity_results$actual_values), 
               length(extreme_velocity_results$predictions))),
  Value = c(extreme_velocity_results$actual_values, 
            extreme_velocity_results$predictions)
)

extreme_vel_plot <- ggplot(extreme_vel_data, aes(x = Value, fill = Type)) + 
  geom_density(size = 1.2) +
  theme_bw() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(text = element_text(size = 14)) +
  labs(title = "Extreme Velocity: Actual vs Predicted", 
       x = "Velocity [m/h]", 
       y = "Density")

# Save extreme value plots
ggsave("data/results/figures/extreme_displacement_comparison.png", plot = extreme_disp_plot)
ggsave("data/results/figures/extreme_velocity_comparison.png", plot = extreme_vel_plot)

# Create a summary table of model metrics
metrics_summary <- data.frame(
  Model = c(
    "Typical Displacement", 
    "Extreme Displacement", 
    "Typical Velocity", 
    "Extreme Velocity", 
    "Typical Displacement with Tracer Dimensions",
    "Extreme Displacement with Tracer Dimensions",
    "Typical Velocity with Tracer Dimensions",
    "Extreme Velocity with Tracer Dimensions"
  ),
  Observations = c(
    length(typical_displacement_results$actual_values),
    length(extreme_displacement_results$actual_values),
    length(typical_velocity_results$actual_values),
    length(extreme_velocity_results$actual_values),
    length(typical_displacement_with_dims_results$actual_values),
    length(extreme_displacement_with_dims_results$actual_values),
    length(typical_velocity_with_dims_results$actual_values),
    length(extreme_velocity_with_dims_results$actual_values)
  ),
  RMSE = c(
    typical_displacement_results$metrics$rmse,
    extreme_displacement_results$metrics$rmse,
    typical_velocity_results$metrics$rmse,
    extreme_velocity_results$metrics$rmse,
    typical_displacement_with_dims_results$metrics$rmse,
    extreme_displacement_with_dims_results$metrics$rmse,
    typical_velocity_with_dims_results$metrics$rmse,
    extreme_velocity_with_dims_results$metrics$rmse
  ),
  MAE = c(
    typical_displacement_results$metrics$mae,
    extreme_displacement_results$metrics$mae,
    typical_velocity_results$metrics$mae,
    extreme_velocity_results$metrics$mae,
    typical_displacement_with_dims_results$metrics$mae,
    extreme_displacement_with_dims_results$metrics$mae,
    typical_velocity_with_dims_results$metrics$mae,
    extreme_velocity_with_dims_results$metrics$mae
  ),
  Wasserstein = c(
    typical_displacement_results$metrics$wasserstein,
    extreme_displacement_results$metrics$wasserstein,
    typical_velocity_results$metrics$wasserstein,
    extreme_velocity_results$metrics$wasserstein,
    typical_displacement_with_dims_results$metrics$wasserstein,
    extreme_displacement_with_dims_results$metrics$wasserstein,
    typical_velocity_with_dims_results$metrics$wasserstein,
    extreme_velocity_with_dims_results$metrics$wasserstein
  )
)

# Add baseline metrics for comparison
baseline_metrics <- data.frame(
  Model = c(
    "Baseline: Typical Displacement", 
    "Baseline: Extreme Displacement", 
    "Baseline: Typical Velocity", 
    "Baseline: Extreme Velocity"
  ),
  Observations = c(
    length(typical_displacement_results$actual_values),
    length(extreme_displacement_results$actual_values),
    length(typical_velocity_results$actual_values),
    length(extreme_velocity_results$actual_values)
  ),
  RMSE = c(
    typical_displacement_results$baseline$rmse,
    extreme_displacement_results$baseline$rmse,
    typical_velocity_results$baseline$rmse,
    extreme_velocity_results$baseline$rmse
  ),
  MAE = c(
    typical_displacement_results$baseline$mae,
    extreme_displacement_results$baseline$mae,
    typical_velocity_results$baseline$mae,
    extreme_velocity_results$baseline$mae
  ),
  Wasserstein = NA
)

# Combine model and baseline metrics
all_metrics <- rbind(metrics_summary, baseline_metrics)
all_metrics[, c("RMSE", "MAE", "Wasserstein")] <- round(all_metrics[, c("RMSE", "MAE", "Wasserstein")], 4)

# Print metrics table
cat("\n====== XGBoost Model Performance Metrics ======\n\n")
print(all_metrics, row.names = FALSE)

# Save results to CSV
write.csv(all_metrics, "data/results/all_error_metrics.csv", row.names = FALSE)

# Calculate improvement over baseline
improvement_metrics <- data.frame(
  Model = metrics_summary$Model,
  RMSE_Improvement = round((baseline_metrics$RMSE - metrics_summary$RMSE) / baseline_metrics$RMSE * 100, 2),
  MAE_Improvement = round((baseline_metrics$MAE - metrics_summary$MAE) / baseline_metrics$MAE * 100, 2)
)

cat("\n====== Improvement Over Baseline (%) ======\n\n")
print(improvement_metrics, row.names = FALSE)

# Add comparison between with/without tracer dimensions
cat("\n====== Comparison: With vs. Without Tracer Dimensions ======\n\n")

dimension_comparison <- data.frame(
  Model_Pair = c(
    "Typical Displacement", 
    "Extreme Displacement", 
    "Typical Velocity", 
    "Extreme Velocity"
  ),
  RMSE_Difference = c(
    typical_displacement_results$metrics$rmse - typical_displacement_with_dims_results$metrics$rmse,
    extreme_displacement_results$metrics$rmse - extreme_displacement_with_dims_results$metrics$rmse,
    typical_velocity_results$metrics$rmse - typical_velocity_with_dims_results$metrics$rmse,
    extreme_velocity_results$metrics$rmse - extreme_velocity_with_dims_results$metrics$rmse
  ),
  RMSE_Percent_Change = c(
    (typical_displacement_results$metrics$rmse - typical_displacement_with_dims_results$metrics$rmse) / typical_displacement_results$metrics$rmse * 100,
    (extreme_displacement_results$metrics$rmse - extreme_displacement_with_dims_results$metrics$rmse) / extreme_displacement_results$metrics$rmse * 100,
    (typical_velocity_results$metrics$rmse - typical_velocity_with_dims_results$metrics$rmse) / typical_velocity_results$metrics$rmse * 100,
    (extreme_velocity_results$metrics$rmse - extreme_velocity_with_dims_results$metrics$rmse) / extreme_velocity_results$metrics$rmse * 100
  ),
  MAE_Difference = c(
    typical_displacement_results$metrics$mae - typical_displacement_with_dims_results$metrics$mae,
    extreme_displacement_results$metrics$mae - extreme_displacement_with_dims_results$metrics$mae,
    typical_velocity_results$metrics$mae - typical_velocity_with_dims_results$metrics$mae,
    extreme_velocity_results$metrics$mae - extreme_velocity_with_dims_results$metrics$mae
  ),
  MAE_Percent_Change = c(
    (typical_displacement_results$metrics$mae - typical_displacement_with_dims_results$metrics$mae) / typical_displacement_results$metrics$mae * 100,
    (extreme_displacement_results$metrics$mae - extreme_displacement_with_dims_results$metrics$mae) / extreme_displacement_results$metrics$mae * 100,
    (typical_velocity_results$metrics$mae - typical_velocity_with_dims_results$metrics$mae) / typical_velocity_results$metrics$mae * 100,
    (extreme_velocity_results$metrics$mae - extreme_velocity_with_dims_results$metrics$mae) / extreme_velocity_results$metrics$mae * 100
  )
)

# Round numeric columns for dimension comparison
dimension_comparison[, c("RMSE_Difference", "RMSE_Percent_Change", "MAE_Difference", "MAE_Percent_Change")] <- 
  round(dimension_comparison[, c("RMSE_Difference", "RMSE_Percent_Change", "MAE_Difference", "MAE_Percent_Change")], 4)

# Add interpretation 
dimension_comparison$Interpretation <- ifelse(
  dimension_comparison$RMSE_Percent_Change > 0,
  "Adding tracer dimensions IMPROVES performance",
  "Adding tracer dimensions WORSENS performance"
)

print(dimension_comparison, row.names = FALSE)



# -------------------------------------------------------------------------
# 5. Model Validation
# -------------------------------------------------------------------------

cat("Validating models...\n")

# Reformat model results to expected structure for validation functions
models <- list(
  typical_displacement = typical_displacement_results,
  extreme_displacement = extreme_displacement_results,
  typical_velocity = typical_velocity_results,
  extreme_velocity = extreme_velocity_results,
  typical_velocity_displacement = typical_velocity_displacement_results
)


# Cross-validation metrics
cv_results <- validate_models(models)

cat("\nCross-validation metrics:\n")
print(cv_results$cv_metrics)

# Train movement classification model
cat("Training movement classification model...\n")

# Create binary target for movement classification (moved/not moved)
df_models_$moved <- as.numeric(df_models_$graph_dist > 0)

# Train model to predict movement probability
movement_classification_results <- analyze_classification(
  data = df_models_,
  target_col = "moved",
  feature_cols = feature_cols
)

# Get feature importance
movement_importance <- xgb.importance(model = movement_classification_results$model)

# Print top 5 features for movement classification
cat("\nTop 5 features for movement classification model:\n")
print(head(movement_importance[order(-movement_importance$Gain), c("Feature", "Gain")], 5))

# Add classification metrics to summary
classification_metrics <- data.frame(
  Model = "Movement Classification",
  Observations = movement_classification_results$n_observations,
  AUC = round(movement_classification_results$metrics$auc, 4),
  Accuracy = round(movement_classification_results$metrics$accuracy, 4),
  Kappa = round(movement_classification_results$metrics$kappa, 4)
)

# Print classification metrics
cat("\n====== Movement Classification Model Performance ======\n\n")
print(classification_metrics, row.names = FALSE)

# Add the movement classification model to the models list
models$movement_classification <- movement_classification_results


# Validate on non-consecutive data
cat("\nValidating on non-consecutive observations...\n")

# Create the proper dataset for non-consecutive validation
non_consecutive_dataset <- df_ch %>% 
  left_join(df_loc_) %>% 
  mutate(Event = EventoEnd) %>% 
  left_join(baseline_weather %>% rename(Event = EventoStart)) %>% 
  drop_na()

# Run the validation
non_consecutive_results <- validate_on_non_consecutive_data(
  tracer_data = non_consecutive_dataset,
  baseline_weather = baseline_weather,
  models = models,
  feature_cols = feature_cols,  # Pass the same feature columns used for training
  displacement_col = "graph_dist",
  velocity_col = "graph_vel"
)

# Display non-consecutive validation results
cat("\nNon-consecutive tracer validation results:\n")
print(non_consecutive_results$metrics)

# Create visualization for non-consecutive validation
cat("Creating non-consecutive validation visualizations...\n")

non_consecutive_plot <- plot_non_consecutive_density(non_consecutive_results)
ggsave("data/results/figures/non_consecutive_density_validation.png", plot = non_consecutive_plot)
# -------------------------------------------------------------------------
# 6. Bed Load Volume Estimation
# -------------------------------------------------------------------------
cat("Estimating bed load volumes...\n")

# Parameters for volume estimation
parameters <- list(
  active_width = 6,                # River width (m)
  porosity = 0.3,                  # Porosity coefficient 
  extrapolation_factor = 2.1,      # Factor for unobserved events
  typical_proportion = 0.88,       # Weight for typical events (88%)
  extreme_proportion = 0.12,       # Weight for extreme events (12%)
  d90 = 0.10,                      # d90 grain size (m)
  mild_event_depth = 0.10,         # Active layer depth for mild events (d90)
  typical_event_depth = 0.15,      # Active layer depth for typical events (1.5*d90)
  intense_event_depth = 0.20,      # Active layer depth for intense events (2*d90)
  basin_volume_per_year = 600      # m³ per year from retention basin
)

# Estimate bed load volumes - passing feature_cols to ensure feature matching
volumes <- estimate_bed_load_volumes(
  data = df_models_,
  pool_grid = pool_grid,
  models = models,
  parameters = parameters,
  monitoring_period_years = 3,
  feature_cols = feature_cols  # Pass the same features used for training
)

cat("\nBed load volume estimates:\n")
cat("Typical velocity contribution:", round(volumes$typical_volume, 2), "m³\n")
cat("Extreme velocity contribution:", round(volumes$extreme_volume, 2), "m³\n")
cat("Total estimated volume:", round(volumes$total_volume, 2), "m³\n")
cat("Expected volume from retention basin:", round(volumes$expected_volume, 2), "m³\n")
cat("Ratio of estimated to expected:", 
    round(volumes$total_volume / volumes$expected_volume, 2), "\n")

# Save results to CSV
write.csv(volumes$volume_by_event, "data/results/volume_by_event.csv", row.names = FALSE)

cat("Bed load volume estimation completed.\n")

# -------------------------------------------------------------------------
# 7. Comprehensive Visualizations
# -------------------------------------------------------------------------

cat("Creating visualizations...\n")

# Visualize Markov Chain transitions
plot_markov_transitions(mc_two_state, mc_displacement, mc_velocity)

# Plot feature importance
plot_feature_importance(
  typical_disp_importance,
  typical_vel_importance,
  extreme_disp_importance,
  extreme_vel_importance
)


# Plot bed load volume contributions by event magnitudes
plot_bed_load_contributions(volumes)

# -------------------------------------------------------------------------
# 8. Save Results
# -------------------------------------------------------------------------

cat("Saving results...\n")

# Create results directory if it doesn't exist
dir.create("data/results/models", showWarnings = FALSE)

# Save models
saveRDS(models, "data/results/models/xgboost_models.rds")

# Save processed data
saveRDS(df_models_, "data/results/data/processed_models_data.rds")

# Save Markov Chains
saveRDS(
  list(
    two_state = mc_two_state,
    displacement = mc_displacement,
    velocity = mc_velocity
  ),
  "data/results/data/markov_chains.rds"
)


# Save all validation results
saveRDS(
  list(
    cv_results = cv_results,
    non_consecutive = non_consecutive_results
  ),
  "data/results/data/validation_results.rds"
)

cat("\nAnalysis complete! Results saved to the 'results' directory.\n")

