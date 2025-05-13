# R/utils/plotting.R

plot_markov_transitions <- function(mc_two_state, mc_displacement, mc_velocity) {
  # Plot Markov Chain transitions as heatmaps
  
  # Convert transition matrices to data frames
  two_state_df <- as.data.frame(mc_two_state@transitionMatrix)
  two_state_df$From <- rownames(two_state_df)
  two_state_df <- two_state_df %>%
    pivot_longer(cols = !From, names_to = "To", values_to = "Probability")
  
  displacement_df <- as.data.frame(mc_displacement@transitionMatrix)
  displacement_df$From <- rownames(displacement_df)
  displacement_df <- displacement_df %>%
    pivot_longer(cols = !From, names_to = "To", values_to = "Probability")
  
  velocity_df <- as.data.frame(mc_velocity@transitionMatrix)
  velocity_df$From <- rownames(velocity_df)
  velocity_df <- velocity_df %>%
    pivot_longer(cols = !From, names_to = "To", values_to = "Probability")
  
  # Create heatmaps
  p1 <- ggplot(two_state_df, aes(x = To, y = From, fill = Probability)) +
    geom_tile() +
    geom_text(aes(label = round(Probability, 2)), color = "white") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = "Two-State Markov Chain Transitions",
         x = "To", y = "From") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(displacement_df, aes(x = To, y = From, fill = Probability)) +
    geom_tile() +
    geom_text(aes(label = round(Probability, 2)), color = "white") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = "Displacement State Transitions",
         x = "To", y = "From") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p3 <- ggplot(velocity_df, aes(x = To, y = From, fill = Probability)) +
    geom_tile() +
    geom_text(aes(label = round(Probability, 2)), color = "white") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(title = "Velocity State Transitions",
         x = "To", y = "From") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plots
  ggsave("data/results/figures/two_state_transitions.png", p1, width = 8, height = 6)
  ggsave("data/results/figures/displacement_transitions.png", p2, width = 10, height = 6)
  ggsave("data/results/figures/velocity_transitions.png", p3, width = 10, height = 6)
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}

plot_feature_importance <- function(typical_disp_importance, typical_vel_importance,
                                    extreme_disp_importance, extreme_vel_importance) {
  # Plot feature importance for all models
  
  # Prepare data frames
  typical_disp_df <- typical_disp_importance %>%
    arrange(desc(Gain)) %>%
    head(10) %>%
    mutate(Model = "Typical Displacement")
  
  typical_vel_df <- typical_vel_importance %>%
    arrange(desc(Gain)) %>%
    head(10) %>%
    mutate(Model = "Typical Velocity")
  
  extreme_disp_df <- extreme_disp_importance %>%
    arrange(desc(Gain)) %>%
    head(10) %>%
    mutate(Model = "Extreme Displacement")
  
  extreme_vel_df <- extreme_vel_importance %>%
    arrange(desc(Gain)) %>%
    head(10) %>%
    mutate(Model = "Extreme Velocity")
  
  # Combine data
  all_importance <- rbind(
    typical_disp_df[, c("Feature", "Gain", "Model")],
    typical_vel_df[, c("Feature", "Gain", "Model")],
    extreme_disp_df[, c("Feature", "Gain", "Model")],
    extreme_vel_df[, c("Feature", "Gain", "Model")]
  )
  
  # Create plot
  p <- ggplot(all_importance, aes(x = reorder(Feature, Gain), y = Gain, fill = Model)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Model, scales = "free_y") +
    coord_flip() +
    labs(title = "Feature Importance Across Models",
         x = "Feature", y = "Gain") +
    theme_minimal() +
    scale_fill_viridis_d()
  
  # Save plot
  ggsave("data/results/figures/feature_importance.png", p, width = 12, height = 8)
  
  return(p)
}

plot_bed_load_contributions <- function(volumes) {
  # Plot bed load volume contributions
  
  # Prepare data
  volume_data <- data.frame(
    Category = c("Typical", "Extreme"),
    Volume = c(volumes$typical_volume, volumes$extreme_volume),
    Percentage = c(
      volumes$typical_volume / volumes$total_volume * 100,
      volumes$extreme_volume / volumes$total_volume * 100
    )
  )
  
  # Create bar plot
  p1 <- ggplot(volume_data, aes(x = Category, y = Volume, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Volume, 1), " m³\n", 
                                 round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    labs(title = "Bed Load Volume Contributions",
         subtitle = paste("Total Volume:", round(volumes$total_volume, 1), "m³"),
         x = "", y = "Volume (m³)") +
    theme_minimal() +
    scale_fill_viridis_d()
  
  # Create pie chart
  p2 <- ggplot(volume_data, aes(x = "", y = Percentage, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = "Proportion of Bed Load Volume",
         subtitle = paste("Total Volume:", round(volumes$total_volume, 1), "m³"),
         x = "", y = "") +
    theme_void() +
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5)) +
    scale_fill_viridis_d()
  
  # Compare with expected volume
  comparison_data <- data.frame(
    Source = c("Model Estimate", "Retention Basin"),
    Volume = c(volumes$total_volume, volumes$expected_volume)
  )
  
  p3 <- ggplot(comparison_data, aes(x = Source, y = Volume, fill = Source)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Volume, 1), " m³")), 
              position = position_stack(vjust = 0.5)) +
    labs(title = "Model Validation: Estimated vs. Measured Volume",
         subtitle = paste("Ratio:", round(volumes$total_volume / volumes$expected_volume, 2)),
         x = "", y = "Volume (m³)") +
    theme_minimal() +
    scale_fill_manual(values = c("Model Estimate" = "#440154", "Retention Basin" = "#21908C"))
  
  # Save plots
  ggsave("data/results/figures/bed_load_volumes.png", p1, width = 8, height = 6)
  ggsave("data/results/figures/bed_load_proportions.png", p2, width = 8, height = 6)
  ggsave("data/results/figures/volume_validation.png", p3, width = 8, height = 6)
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}