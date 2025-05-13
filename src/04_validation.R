# R/04_validation.R

#' @title Validate XGBoost Models
#' @description Calculates validation metrics for all models
#' @param models List containing model results 
#' @return List of validation metrics and wasserstein distances
#' @export
validate_models <- function(models) {
  # Calculate Wasserstein distance
  wasserstein_distances <- list(
    typical_displacement = calculate_wasserstein_distance(
      models$typical_displacement$predictions,
      models$typical_displacement$actual_values
    ),
    extreme_displacement = calculate_wasserstein_distance(
      models$extreme_displacement$predictions,
      models$extreme_displacement$actual_values
    ),
    typical_velocity = calculate_wasserstein_distance(
      models$typical_velocity$predictions,
      models$typical_velocity$actual_values
    ),
    extreme_velocity = calculate_wasserstein_distance(
      models$extreme_velocity$predictions,
      models$extreme_velocity$actual_values
    )
  )
  
  # Compile CV metrics
  cv_metrics <- data.frame(
    Model = c("Typical Displacement", "Extreme Displacement", 
              "Typical Velocity", "Extreme Velocity"),
    RMSE = c(
      models$typical_displacement$metrics$rmse,
      models$extreme_displacement$metrics$rmse,
      models$typical_velocity$metrics$rmse,
      models$extreme_velocity$metrics$rmse
    ),
    MAE = c(
      models$typical_displacement$metrics$mae,
      models$extreme_displacement$metrics$mae,
      models$typical_velocity$metrics$mae,
      models$extreme_velocity$metrics$mae
    ),
    Wasserstein = c(
      wasserstein_distances$typical_displacement,
      wasserstein_distances$extreme_displacement,
      wasserstein_distances$typical_velocity,
      wasserstein_distances$extreme_velocity
    )
  )
  
  return(list(
    cv_metrics = cv_metrics,
    wasserstein_distances = wasserstein_distances
  ))
}

#' @title Calculate Wasserstein Distance
#' @description Calculates the Wasserstein distance between two distributions
#' @param predicted Vector of predicted values
#' @param actual Vector of actual values
#' @param p Order of the Wasserstein distance (default: 2)
#' @return Wasserstein distance value
#' @export
calculate_wasserstein_distance <- function(predicted, actual, p = 2) {
  # Use transport package for reliable calculation
  tryCatch({
    transport::wasserstein1d(predicted, actual, p = p)
  }, error = function(e) {
    warning("Wasserstein distance calculation failed: ", e$message)
    return(NA)
  })
}


#' @title Validate Models on Non-Consecutive Tracer Observations
#' @description Tests model performance on tracers with non-consecutive observations
#' @param tracer_data Full dataset containing all required columns
#' @param models List containing all model objects
#' @param feature_cols Vector of feature column names used by all models
#' @param displacement_col Name of displacement column
#' @param velocity_col Name of velocity column
#' @return List with validation metrics and predictions
#' @export
validate_on_non_consecutive_data <- function(tracer_data, 
                                             baseline_weather,
                                             models,
                                             feature_cols,
                                             displacement_col = "graph_dist",
                                             velocity_col = "graph_vel") {
  
  # Extract tracers with non-consecutive observations
  non_consecutive_tracers <- tracer_data %>% 
    filter(X_Start >= X_End | Y_Start >= Y_End | graph_dist <= 1) %>%
    filter(EventoStart < EventoEnd - 1) %>%
    filter(graph_dist > 0) %>%
    select(IDREF, EventoStart, EventoEnd, graph_dist) %>%
    mutate(tracer_id = seq_len(nrow(.)))
  
  # If no non-consecutive data found
  if(nrow(non_consecutive_tracers) == 0) {
    warning("No non-consecutive tracer observations found for validation.")
    return(list(
      metrics = data.frame(
        RMSE = NA,
        MAE = NA,
        Wasserstein_Distance = NA
      ),
      predictions = NULL,
      actual = NULL
    ))
  }
  # Extract tracer data for non-consecutive observations
  nonconsecutive_tracers <- tracer_data %>% 
    filter(X_Start >= X_End | Y_Start >= Y_End | graph_dist <= 1) %>%
    filter(EventoStart < EventoEnd - 1) %>%
    filter(graph_dist > 0) %>%
    select(pebble_PC1, pebble_PC2, X_Start, Y_Start, duration, mean_Q, mean_h, 
           starts_with("TipoStart"), IDREF, EventoStart, EventoEnd) %>%
    mutate(tracer_id = seq_len(n()))
  
  # Create expanded dataset with one row per intervening event
  expanded_events <- lapply(
    split(nonconsecutive_tracers, nonconsecutive_tracers$tracer_id),
    function(tracer) {
      start_event <- tracer$EventoStart + 1
      end_event <- tracer$EventoEnd
      
      # Replicate tracer info for each event between start and end
      tracer[rep(1, end_event - start_event + 1), ] %>% 
        mutate(Event = start_event:end_event)
    }
  ) %>% 
    bind_rows() %>% 
    select(-duration, -mean_Q, -mean_h)
  
  # Get unique hydrological data
  hydrology_data <- tracer_data %>% 
    select(Event = EventoStart, duration, mean_Q, mean_h, weather_PC1, weather_PC2) %>%
    distinct()
  
  # Join with hydrological data
  events_with_hydrology <- expanded_events %>% 
    select(-EventoStart, -EventoEnd) %>% 
    left_join(baseline_weather, by = c("Event"= "EventoStart"))
  
  # Add morphological unit indicator variables
  events_with_hydrology <- events_with_hydrology %>% 
    bind_cols(
      model.matrix(~ TipoStart - 1, data = events_with_hydrology) %>%
        as.data.frame()
    )
  
  # Prepare prediction data using feature_cols
  prediction_matrix <- events_with_hydrology %>%
    select(all_of(intersect(feature_cols, colnames(events_with_hydrology)))) %>%
    as.matrix()
  
  # Make predictions for each step in the tracer's path
  # 1. Movement probability prediction
  movement_probabilities <- predict(
    models$movement_classification$model, 
    prediction_matrix
  )
  
  # 2. Typical velocity prediction
  typical_velocities <- predict(
    models$typical_velocity_displacement$model, 
    prediction_matrix
  )
  typical_velocities[typical_velocities < 0] <- 0
  
  # Calculate predicted displacement for each event
  event_displacements <- events_with_hydrology %>% 
    select(IDREF, tracer_id) %>% 
    mutate(
      typical_event_displacement = events_with_hydrology$duration * typical_velocities * movement_probabilities,
    )
  
  # Sum displacements by tracer to get total movement
  tracer_total_displacements <- event_displacements %>%
    group_by(IDREF, tracer_id) %>% 
    summarize(
      typical_total_displacement = sum(typical_event_displacement),
      .groups = "drop"
    ) %>%
    arrange(tracer_id)
  
  combined_displacements <- tracer_total_displacements$typical_total_displacement 
  # Calculate validation metrics
  rmse <- MLmetrics::RMSE(combined_displacements, non_consecutive_tracers$graph_dist)
  mae <- MLmetrics::MAE(combined_displacements, non_consecutive_tracers$graph_dist)
  
  # Calculate Wasserstein distance for distribution similarity
  wasserstein_dist <- transport::wasserstein1d(
    combined_displacements, 
    non_consecutive_tracers$graph_dist,
    p = 2
  )
  
  # Create validation results dataframe
  validation_results <- data.frame(
    IDREF = non_consecutive_tracers$IDREF,
    ActualDisplacement = non_consecutive_tracers$graph_dist,
    PredictedDisplacement = combined_displacements,
    PredictionError = non_consecutive_tracers$graph_dist - combined_displacements,
    EventStart = non_consecutive_tracers$EventoStart,
    EventEnd = non_consecutive_tracers$EventoEnd,
    EventGap = non_consecutive_tracers$EventoEnd - non_consecutive_tracers$EventoStart - 1
  )
  
  # Return comprehensive results
  return(list(
    metrics = data.frame(
      RMSE = rmse,
      MAE = mae,
      Wasserstein_Distance = wasserstein_dist
    ),
    predictions = combined_displacements,
    actual = non_consecutive_tracers$graph_dist,
    validation_data = validation_results
  ))
}


plot_non_consecutive_density <- function(non_consecutive_results) {
  # Prepare data for plotting
  plot_data <- data.frame(
    Value = c(
      non_consecutive_results$validation_data$ActualDisplacement,
      non_consecutive_results$validation_data$PredictedDisplacement
    ),
    Type = factor(c(
      rep("Actual Displacement", nrow(non_consecutive_results$validation_data)),
      rep("Predicted Displacement", nrow(non_consecutive_results$validation_data))
    ))
  )
  
  # Create density plot
  non_consecutive_plot <- ggplot(plot_data, aes(x = Value, fill = Type)) + 
    geom_density(alpha = 0.5, size = 12) +
    theme_bw() +
    scale_fill_viridis_d(alpha = 0.5) +
    theme(text = element_text(size = 14)) +
    labs(
      title = "Non-consecutive Tracer: Actual vs Predicted Displacement", 
      x = "Displacement [m]", 
      y = "Density"
    ) +
    annotate(
      "text", 
      x = max(plot_data$Value) * 0.7, 
      y = max(density(plot_data$Value)$y) * 0.9,
      label = paste(
        "RMSE:", 
        round(non_consecutive_results$metrics$RMSE, 2),
        "\nMAE:", 
        round(non_consecutive_results$metrics$MAE, 2)
      )
    )
}


#' @title Estimate Bed Load Volumes
#' @description Estimates sediment transport volumes using virtual velocity approach
#' following Haschenburger & Church (1998), Liebault & Laronne (2008), and Cazzador et al. (2020)
#' @param data Full dataset with tracer information
#' @param pool_grid Grid dataset covering the entire reach
#' @param models List containing typical_velocity, extreme_velocity, and movement_classification models
#' @param parameters List of parameters for volume calculation (optional)
#' @param feature_cols Vector of feature column names to use for prediction
#' @param monitoring_period_years Total monitoring period in years (optional)
#' @return List with volume estimates and event-specific details
#' @export
estimate_bed_load_volumes <- function(data, pool_grid, models, parameters, feature_cols, monitoring_period_years = NULL) {
  # Default parameters if not provided, based on paper methodology
  if (missing(parameters)) {
    parameters <- list(
      active_width = 6,                # River width (m)
      porosity = 0.3,                  # Porosity coefficient 
      extrapolation_factor = 2.1,      # Factor for unobserved events
      typical_proportion = 0.87,       # Weight for typical events (87%)
      extreme_proportion = 0.13,       # Weight for extreme events (13%)
      d90 = 0.10,                      # d90 grain size (m)
      mild_event_depth = 0.10,         # Active layer depth for mild events (d90)
      typical_event_depth = 0.15,      # Active layer depth for typical events (1.5*d90)
      intense_event_depth = 0.20,      # Active layer depth for intense events (2*d90)
      basin_volume_per_year = 600      # m³ per year from retention basin
    )
  }
  
  # Direct event classification map
  event_classification <- c(
    "1" = "typical",  # Event 1
    "2" = "mild",     # Event 2
    "3" = "typical",  # Event 3
    "4" = "typical",  # Event 4
    "5" = "intense",  # Event 5
    "6" = "typical",     # Event 6
    "7" = "typical",  # Event 7
    "8" = "typical",  # Event 8
    "9" = "typical",  # Event 9
    "10" = "mild",    # Event 10
    "11" = "typical", # Event 11
    "12" = "typical", # Event 12
    "13" = "typical", # Event 13
    "14" = "typical",    # Event 14
    "15" = "typical", # Event 15
    "16" = "mild", # Event 16
    "17" = "intense",    # Event 17
    "18" = "typical", # Event 18
    "19" = "typical", # Event 19
    "20" = "typical", # Event 20
    "21" = "typical",    # Event 21
    "22" = "typical", # Event 22
    "23" = "typical", # Event 23
    "24" = "intense",    # Event 24
    "25" = "intense", # Event 25
    "26" = "intense", # Event 26
    "27" = "typical"  # Event 27
  )
  
  # Get unique events
  events <- data %>% 
    dplyr::select(EventoEnd, duration) %>% 
    dplyr::distinct() %>%
    dplyr::arrange(EventoEnd)
  
  # Add event classification to events data
  events <- events %>%
    dplyr::mutate(event_type = event_classification[as.character(EventoEnd)])
  
  # Assign active layer depth based on event type
  events <- events %>%
    dplyr::mutate(active_depth = dplyr::case_when(
      event_type == "intense" ~ parameters$intense_event_depth,
      event_type == "mild" ~ parameters$mild_event_depth,
      TRUE ~ parameters$typical_event_depth
    ))
  
  # Create standard tracer for predictions
  standard_tracer <- data.frame(
    pebble_PC1 = median(data$pebble_PC1, na.rm = TRUE),
    pebble_PC2 = median(data$pebble_PC2, na.rm = TRUE)
  )
  
  # Prepare results dataframe
  volume_estimates <- data.frame(
    EventoEnd = numeric(),
    Event_Type = character(),
    Duration = numeric(),
    Active_Depth = numeric(),
    Probability_Movement = numeric(),
    Typical_Velocity = numeric(),
    Extreme_Velocity = numeric(),
    Weighted_Velocity = numeric(),
    Typical_Volume = numeric(),
    Extreme_Volume = numeric(),
    Total_Volume = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get morphological unit columns
  morphological_units <- grep("^TipoStart", names(pool_grid), value = TRUE)
  
  # Prepare base prediction grid - rename columns to match model expectations
  pred_base <- pool_grid
  names(pred_base)[names(pred_base) == "x"] <- "X_Start"
  names(pred_base)[names(pred_base) == "y"] <- "Y_Start"
  
  # Add standard tracer characteristics
  pred_base$pebble_PC1 <- standard_tracer$pebble_PC1
  pred_base$pebble_PC2 <- standard_tracer$pebble_PC2
  
  # Loop through each event
  for (i in 1:nrow(events)) {
    event_id <- events$EventoEnd[i]
    event_duration <- events$duration[i]
    event_active_depth <- events$active_depth[i]
    event_type <- events$event_type[i]
    
    # Create prediction grid for this event
    pred_grid <- pred_base
    pred_grid$duration <- event_duration
    
    # Get event data for any needed metrics
    event_data <- data %>% dplyr::filter(EventoEnd == event_id)
    
    # Add any additional features from event data needed for prediction
    for (col in setdiff(feature_cols, names(pred_grid))) {
      if (col %in% names(event_data)) {
        pred_grid[[col]] <- median(event_data[[col]], na.rm = TRUE)
      }
    }
    
    # Select features required by the models
    pred_matrix <- as.matrix(pred_grid[, feature_cols])
    
    # Make predictions for all grid cells
    movement_probabilities <- predict(models$movement_classification$model, pred_matrix)
    typical_velocities <- predict(models$typical_velocity$model, pred_matrix)
    extreme_velocities <- predict(models$extreme_velocity$model, pred_matrix)
    
    # Ensure values are in proper range
    movement_probabilities <- pmax(0, pmin(1, movement_probabilities))
    typical_velocities <- pmax(0, typical_velocities)
    extreme_velocities <- pmax(0, extreme_velocities)
    
    # Calculate morphological sector statistics
    morpho_stats <- data.frame(
      unit = morphological_units,
      cell_count = sapply(morphological_units, function(u) sum(pred_grid[[u]] == 1)),
      stringsAsFactors = FALSE
    )
    
    # Calculate mean values for each morphological unit
    for (j in 1:nrow(morpho_stats)) {
      unit_cells <- which(pred_grid[[morpho_stats$unit[j]]] == 1)
      if (length(unit_cells) > 0) {
        morpho_stats$prob_movement[j] <- mean(movement_probabilities[unit_cells], na.rm = TRUE)
        morpho_stats$typical_velocity[j] <- mean(typical_velocities[unit_cells], na.rm = TRUE)
        morpho_stats$extreme_velocity[j] <- mean(extreme_velocities[unit_cells], na.rm = TRUE)
      } else {
        morpho_stats$prob_movement[j] <- 0
        morpho_stats$typical_velocity[j] <- 0
        morpho_stats$extreme_velocity[j] <- 0
      }
    }
    
    # Calculate area proportions
    morpho_stats$area_prop <- morpho_stats$cell_count / sum(morpho_stats$cell_count)
    
    # Calculate weighted averages
    avg_prob_movement <- sum(morpho_stats$prob_movement * morpho_stats$area_prop, na.rm = TRUE)
    avg_typical_velocity <- sum(morpho_stats$typical_velocity * morpho_stats$area_prop, na.rm = TRUE)
    avg_extreme_velocity <- sum(morpho_stats$extreme_velocity * morpho_stats$area_prop, na.rm = TRUE)
    
    # Weighted average velocity based on data proportions
    weighted_velocity <- (avg_typical_velocity * parameters$typical_proportion) + 
      (avg_extreme_velocity * parameters$extreme_proportion)
    
    # Calculate volumes for this event - separate by type
    typical_volume <- avg_typical_velocity * 
      avg_prob_movement * 
      event_duration * 
      event_active_depth * 
      parameters$active_width * 
      (1 - parameters$porosity)
    
    extreme_volume <- avg_extreme_velocity * 
      avg_prob_movement * 
      event_duration * 
      event_active_depth * 
      parameters$active_width * 
      (1 - parameters$porosity)
    
    # Weighted total volume
    total_volume <- weighted_velocity * 
      avg_prob_movement * 
      event_duration * 
      event_active_depth * 
      parameters$active_width * 
      (1 - parameters$porosity)
    
    # Add to results
    volume_estimates <- rbind(volume_estimates, data.frame(
      EventoEnd = event_id,
      Event_Type = event_type,
      Duration = event_duration,
      Active_Depth = event_active_depth,
      Probability_Movement = avg_prob_movement,
      Typical_Velocity = avg_typical_velocity,
      Extreme_Velocity = avg_extreme_velocity,
      Weighted_Velocity = weighted_velocity,
      Typical_Volume = typical_volume,
      Extreme_Volume = extreme_volume,
      Total_Volume = total_volume
    ))
  }
  
  # Sum volumes for all events
  total_typical_volume <- sum(volume_estimates$Typical_Volume, na.rm = TRUE)
  total_extreme_volume <- sum(volume_estimates$Extreme_Volume, na.rm = TRUE)
  total_volume <- sum(volume_estimates$Total_Volume, na.rm = TRUE)
  
  # Apply extrapolation factor for unobserved events
  extrapolated_typical_volume <- total_typical_volume * parameters$extrapolation_factor
  extrapolated_extreme_volume <- total_extreme_volume * parameters$extrapolation_factor
  extrapolated_volume <- total_volume * parameters$extrapolation_factor
  
  # Apply model proportions to get weighted estimates
  weighted_typical_volume <- extrapolated_typical_volume * parameters$typical_proportion
  weighted_extreme_volume <- extrapolated_extreme_volume * parameters$extreme_proportion
  
  # Expected volume calculation
  expected_volume <- NULL
  if (!is.null(monitoring_period_years)) {
    expected_volume <- parameters$basin_volume_per_year * monitoring_period_years
  }
  
  return(list(
    typical_volume = weighted_typical_volume,
    extreme_volume = weighted_extreme_volume,
    total_volume = extrapolated_volume,
    expected_volume = expected_volume,
    volume_by_event = volume_estimates,
    monitoring_period_years = monitoring_period_years
  ))
}
