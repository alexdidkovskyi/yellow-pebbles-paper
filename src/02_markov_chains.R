# R/02_markov_chains.R

#' Markov Chain Analysis for Pebble Movement
#'
#' This script provides functions for analyzing pebble transport as a Markov process,
#' examining how pebbles transition between different states of movement across events.
#'
#' @author Yellow Pebbles Research Team
#' @date March 2025

#' Define "not in" operator for convenience
'%!in%' <- function(x, y) !('%in%'(x, y))

#' Find optimal thresholds using Jenks Natural Breaks
#'
#' @param data Data frame containing movement metrics
#' @param displacement_col Column name for displacement measurements
#' @param velocity_col Column name for velocity measurements
#' @return List with optimal thresholds for displacement and velocity
find_optimal_thresholds <- function(data, displacement_col = "graph_dist", velocity_col = "graph_vel") {
  # For displacement
  displacement_breaks <- BAMMtools::getJenksBreaks(data[[displacement_col]][data[[displacement_col]] >0], 3)
  displacement_threshold <- displacement_breaks[2]
  
  # For velocity
  velocity_breaks <- BAMMtools::getJenksBreaks(data[[velocity_col]][data[[velocity_col]] > 0], 3)
  velocity_threshold <- velocity_breaks[2]
  
  return(list(
    displacement_threshold = displacement_threshold,
    velocity_threshold = velocity_threshold
  ))
}

#' Create a 2-state Markov Chain (moved/not moved)
#'
#' @param data Data frame with pebble movement data
#' @param displacement_col Column name for displacement
#' @param exclude_pairs Event pairs to exclude (e.g. "9-10")
#' @return Markov chain object
create_two_state_mc <- function(data, displacement_col = "graph_dist", exclude_pairs = NULL) {
  # First, prepare data by filtering and organizing
  df_MC <- data %>%
    filter(EventoStart == EventoEnd - 1) %>%
    filter(X_Start >= X_End | Y_Start >= Y_End)
  
  # Filter out specific event pairs if provided
  if (!is.null(exclude_pairs)) {
    df_MC <- df_MC %>%
      arrange(IDREF, Event) %>%
      group_by(IDREF) %>%
      mutate(
        prev_event = lag(Event),
        event_pair = paste0(prev_event, "-", Event)
      ) %>%
      filter(!(event_pair %in% exclude_pairs)) %>%
      ungroup()
  }
  
  # Group data by IDREF
  list_of_movements_by_IDREF <- df_MC %>%
    arrange(IDREF, Event) %>%
    group_by(IDREF) %>%
    group_split()
  
  # Create sequences for each pebble's movement history
  list_of_movements <- lapply(list_of_movements_by_IDREF, function(x) {
    return(case_when(
      x[[displacement_col]] > 0 ~ 'Moved',
      x[[displacement_col]] == 0 ~ 'Not Moved'
    ))
  })
  
  # Fit Markov chain using MLE
  mc_fit <- markovchain::markovchainFit(data = list_of_movements, method = 'mle')
  
  return(mc_fit$estimate)
}

#' Create a 3-state Markov Chain for displacement (zero/typical/long)
#'
#' @param data Data frame with pebble movement data
#' @param displacement_col Column name for displacement
#' @param long_threshold Threshold for long displacement
#' @param exclude_pairs Event pairs to exclude
#' @return Markov chain object
create_displacement_mc <- function(data, displacement_col = "graph_dist", 
                                   long_threshold, exclude_pairs = NULL) {
  # First, prepare data by filtering and organizing
  df_MC <- data %>%
    filter(EventoStart == EventoEnd - 1) %>%
    filter(X_Start >= X_End | Y_Start >= Y_End)
  
  # Filter out specific event pairs if provided
  if (!is.null(exclude_pairs)) {
    df_MC <- df_MC %>%
      arrange(IDREF, Event) %>%
      group_by(IDREF) %>%
      mutate(
        prev_event = lag(Event),
        event_pair = paste0(prev_event, "-", Event)
      ) %>%
      filter(!(event_pair %in% exclude_pairs)) %>%
      ungroup()
  }
  
  # Group data by IDREF
  list_of_movements_by_IDREF <- df_MC %>%
    arrange(IDREF, Event) %>%
    group_by(IDREF) %>%
    group_split()
  
  # Create sequences for each pebble's movement history
  list_of_movements <- lapply(list_of_movements_by_IDREF, function(x) {
    return(case_when(
      x[[displacement_col]] >= long_threshold ~ 'Long Displacement',
      x[[displacement_col]] < long_threshold & x[[displacement_col]] > 0 ~ 'Typical Displacement',
      x[[displacement_col]] == 0 ~ 'Zero Displacement'
    ))
  })
  
  # Fit Markov chain using MLE
  mc_fit <- markovchain::markovchainFit(data = list_of_movements, method = 'mle')
  
  return(mc_fit$estimate)
}

#' Create a 3-state Markov Chain for velocity (zero/typical/high)
#'
#' @param data Data frame with pebble movement data
#' @param velocity_col Column name for velocity
#' @param high_threshold Threshold for high velocity
#' @param exclude_pairs Event pairs to exclude
#' @return Markov chain object
create_velocity_mc <- function(data, velocity_col = "graph_vel", 
                               high_threshold, exclude_pairs = NULL) {
  # First, prepare data by filtering and organizing
  df_MC <- data %>%
    filter(EventoStart == EventoEnd - 1) %>%
    filter(X_Start >= X_End | Y_Start >= Y_End)
  
  # Filter out specific event pairs if provided
  if (!is.null(exclude_pairs)) {
    df_MC <- df_MC %>%
      arrange(IDREF, Event) %>%
      group_by(IDREF) %>%
      mutate(
        prev_event = lag(Event),
        event_pair = paste0(prev_event, "-", Event)
      ) %>%
      filter(!(event_pair %in% exclude_pairs)) %>%
      ungroup()
  }
  
  # Group data by IDREF
  list_of_movements_by_IDREF <- df_MC %>%
    arrange(IDREF, Event) %>%
    group_by(IDREF) %>%
    group_split()
  
  # Create sequences for each pebble's movement history
  list_of_movements <- lapply(list_of_movements_by_IDREF, function(x) {
    return(case_when(
      x[[velocity_col]] >= high_threshold ~ 'High Velocity',
      x[[velocity_col]] < high_threshold & x[[velocity_col]] > 0 ~ 'Typical Velocity',
      x[[velocity_col]] == 0 ~ 'Zero Velocity'
    ))
  })
  
  # Fit Markov chain using MLE
  mc_fit <- markovchain::markovchainFit(data = list_of_movements, method = 'mle')
  
  return(mc_fit$estimate)
}

#' Save Markov Chain results to Excel file
#'
#' @param mc_two_state Two-state Markov Chain object
#' @param mc_displacement Displacement Markov Chain object
#' @param mc_velocity Velocity Markov Chain object
#' @param file_path Path to save the Excel file
#' @return NULL (invisibly)
#' @export
save_markov_chain_results <- function(mc_two_state, 
                                      mc_displacement, 
                                      mc_velocity, 
                                      file_path) {
  # Extract transition matrices directly from the S4 objects
  # For markovchain objects, we can access the transition matrix with @
  two_state_matrix <- mc_two_state@transitionMatrix
  displacement_matrix <- mc_displacement@transitionMatrix
  velocity_matrix <- mc_velocity@transitionMatrix
  
  # Format matrices as data frames with row names as first column
  two_state_df <- data.frame(
    From_State = rownames(two_state_matrix),
    as.data.frame(round(two_state_matrix, 4))
  )
  
  displacement_df <- data.frame(
    From_State = rownames(displacement_matrix),
    as.data.frame(round(displacement_matrix, 4))
  )
  
  velocity_df <- data.frame(
    From_State = rownames(velocity_matrix),
    as.data.frame(round(velocity_matrix, 4))
  )
  
  # Create workbook and save
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Two State Markov Chain")
  openxlsx::addWorksheet(wb, "Displacement Markov Chain")
  openxlsx::addWorksheet(wb, "Velocity Markov Chain")
  
  openxlsx::writeData(wb, "Two State Markov Chain", two_state_df)
  openxlsx::writeData(wb, "Displacement Markov Chain", displacement_df)
  openxlsx::writeData(wb, "Velocity Markov Chain", velocity_df)
  
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  
  invisible(NULL)
}


#' Compare pebble movement and velocity between consecutive events
#'
#' This function performs statistical tests to compare both movement probability and
#' velocity between consecutive flood events, preserving the specific hypothesis
#' testing approach from the original analysis.
#'
#' @param data A data frame containing pebble movement data
#' @param event_column The column name for event identifiers (default: "Event")
#' @param graph_dist_column The column name for displacement (default: "graph_dist")
#' @param graph_vel_column The column name for velocity (default: "graph_vel")
#' @param velocity_threshold Upper threshold for valid velocity values (default: 1.1)
#' @return A list containing two data frames: movement_comparison and velocity_comparison
#' @export
compare_consecutive_events <- function(data, 
                                       event_column = "Event",
                                       graph_dist_column = "graph_dist",
                                       graph_vel_column = "graph_vel",
                                       velocity_threshold = 1.1) {
  
  # Pre-define event categories based on expected relationships
  # These correspond to the specific hypotheses in the original code
  equal_events <- c(3, 7, 11, 12, 14, 19, 20)
  greater_events <- c(1, 5, 6, 15, 17, 21, 26)
  lesser_events <- c(2, 4, 8, 10, 16, 22, 23, 24, 25)
  
  # Initialize storage for test results
  chi_tests <- list()
  perm_tests <- list()
  
  # Perform tests for each event pair (1 to 26)
  for (i in 1:26) {
    # Determine alternative hypothesis based on event categorization
    if (i %in% equal_events) {
      altern_ <- "two.sided"
    } else if (i %in% greater_events) {
      altern_ <- "greater"
    } else if (i %in% lesser_events) {
      altern_ <- "less"
    } else {
      altern_ <- "two.sided"  # Default
    }
    
    # Extract movement data for current and next event
    vec_1 <- data[[graph_dist_column]][data[[event_column]] == i] > 0
    vec_2 <- data[[graph_dist_column]][data[[event_column]] == (i + 1)] > 0
    
    # Perform proportion test for movement probability
    chi_tests[[i]] <- prop.test(
      c(sum(vec_1), sum(vec_2)), 
      c(length(vec_1), length(vec_2)), 
      alternative = altern_
    )
    
    # Extract velocity data for moved tracers only
    vec_1_vel <- data[[graph_vel_column]][
      (data[[event_column]] == i) & 
        (data[[graph_dist_column]] > 0) & 
        (data[[graph_vel_column]] < velocity_threshold)
    ]
    
    vec_2_vel <- data[[graph_vel_column]][
      (data[[event_column]] == (i + 1)) & 
        (data[[graph_dist_column]] > 0) & 
        (data[[graph_vel_column]] < velocity_threshold)
    ]
    
    # Perform t-test for velocity if enough data points are available
    if (length(vec_2_vel) > 1 & length(vec_1_vel) > 1) {
      perm_tests[[i]] <- t.test(vec_1_vel, vec_2_vel, alternative = altern_)
    } else {
      perm_tests[[i]] <- NULL
    }
  }
  
  # Create movement comparison data frame
  movement_df <- data.frame(
    Event = 1:26,
    Flow_PC1_vs_Flow_PC1_next = NA,
    Expected = NA, 
    Num_obs = data %>% 
      group_by(!!sym(event_column)) %>% 
      summarise(n_ = n()) %>% 
      filter(!!sym(event_column) <= 26) %>% 
      pull(n_),
    Moved = data %>% 
      group_by(!!sym(event_column)) %>% 
      summarise(n_ = sum(!!sym(graph_dist_column) > 0)) %>% 
      filter(!!sym(event_column) <= 26) %>% 
      pull(n_),
    Num_obs_next = data %>% 
      group_by(!!sym(event_column)) %>% 
      summarise(n_ = n()) %>% 
      filter(!!sym(event_column) >= 2, !!sym(event_column) <= 27) %>% 
      pull(n_),
    Moved_next = data %>% 
      group_by(!!sym(event_column)) %>% 
      summarise(n_ = sum(!!sym(graph_dist_column) > 0)) %>% 
      filter(!!sym(event_column) >= 2, !!sym(event_column) <= 27) %>% 
      pull(n_),
    Alternative_Hypothesis = sapply(1:26, function(i) chi_tests[[i]]$alternative),
    p_val = sapply(1:26, function(i) round(chi_tests[[i]]$p.value, 4))
  )
  
  # Create velocity comparison data frame
  velocity_df <- data.frame(
    Event = 1:26,
    Flow_PC1_vs_Flow_PC1_next = NA,
    Expected = NA, 
    Alternative_Hypothesis = sapply(1:26, function(i) {
      if (length(perm_tests[[i]]) > 0) perm_tests[[i]]$alternative else NA
    }),
    p_val = sapply(1:26, function(i) {
      if (length(perm_tests[[i]]) > 0) round(perm_tests[[i]]$p.value, 4) else NA
    })
  )
  
  # Set descriptive labels for flow PC1 relationships
  movement_df$Flow_PC1_vs_Flow_PC1_next[equal_events] <- "Approximately equal"
  movement_df$Flow_PC1_vs_Flow_PC1_next[greater_events] <- "Greater than"
  movement_df$Flow_PC1_vs_Flow_PC1_next[lesser_events] <- "Less than"
  
  velocity_df$Flow_PC1_vs_Flow_PC1_next[equal_events] <- "Approximately equal"
  velocity_df$Flow_PC1_vs_Flow_PC1_next[greater_events] <- "Greater than"
  velocity_df$Flow_PC1_vs_Flow_PC1_next[lesser_events] <- "Less than"
  
  # Set expected relationship descriptions
  movement_df$Expected[equal_events] <- "Same probability"
  movement_df$Expected[greater_events] <- "Lower probability"
  movement_df$Expected[lesser_events] <- "Higher probability"
  
  velocity_df$Expected[equal_events] <- "Same Velocity"
  velocity_df$Expected[greater_events] <- "Lower Velocity"
  velocity_df$Expected[lesser_events] <- "Higher Velocity"
  
  # Add event_next column and join with summary statistics for velocity analysis
  velocity_df$Event_next <- velocity_df$Event + 1
  
  # Create velocity summary statistics
  velocity_summary <- data %>% 
    filter(!!sym(graph_dist_column) > 0, !!sym(graph_vel_column) < velocity_threshold) %>% 
    group_by(!!sym(event_column)) %>% 
    summarise(
      Num_obs = n(), 
      Mean_vel = mean(!!sym(graph_vel_column)),
      SD_vel = sd(!!sym(graph_vel_column))
    )
  
  # Join velocity summary statistics to velocity comparison
  velocity_df <- velocity_df %>%
    left_join(velocity_summary, by = c("Event" = event_column)) %>%
    left_join(velocity_summary, by = c("Event_next" = event_column), 
              suffix = c("", "_next"))
  
  # Return both analyses as a list
  return(list(
    movement_comparison = movement_df,
    velocity_comparison = velocity_df
  ))
}
#' Summarize and display event comparison results
#'
#' This function displays a summary of the movement and velocity comparison results,
#' including counts of significant results for different event magnitude relationships.
#'
#' @param results List containing movement_comparison and velocity_comparison data frames
#' @param significance_level P-value threshold for significance (default: 0.05)
#' @return Invisibly returns a summary data frame
#' @export
summarize_comparison_results <- function(results, significance_level = 0.05) {
  # Extract data frames
  movement_df <- results$movement_comparison
  velocity_df <- results$velocity_comparison
  
  # Round p-values to two decimals
  movement_df$p_val <- round(movement_df$p_val, 2)
  velocity_df$p_val <- round(velocity_df$p_val, 2)
  
  # Create a unified results data frame for display
  combined_results <- data.frame(
    Event_N = movement_df$Event,
    Event_N1 = movement_df$Event + 1,
    Magnitude_Relation = case_when(
      movement_df$Flow_PC1_vs_Flow_PC1_next == "Greater than" ~ ">",
      movement_df$Flow_PC1_vs_Flow_PC1_next == "Approximately equal" ~ "=",
      movement_df$Flow_PC1_vs_Flow_PC1_next == "Less than" ~ "<",
      TRUE ~ NA_character_
    ),
    Expected_Movement = movement_df$Expected,
    Movement_pvalue = movement_df$p_val,
    Expected_Velocity = velocity_df$Expected,
    Velocity_pvalue = velocity_df$p_val
  )
  
  # Remove rows with NA values
  combined_results <- combined_results %>% tidyr::drop_na()
  
  # Display results
  cat("\nPairwise comparison of events:\n")
  print(combined_results)
  
  # Count significant results by magnitude relation
  low_to_high <- combined_results$Magnitude_Relation == ">"
  similar <- combined_results$Magnitude_Relation == "="
  high_to_low <- combined_results$Magnitude_Relation == "<"
  
  # Create summary table
  summary_table <- data.frame(
    Magnitude_Relation = c("Low-to-high (>)", "Similar (=)", "High-to-low (<)"),
    Total_Count = c(sum(low_to_high), sum(similar), sum(high_to_low)),
    Significant_Movement = c(
      sum(combined_results$Movement_pvalue[low_to_high] < significance_level, na.rm = TRUE),
      sum(combined_results$Movement_pvalue[similar] < significance_level, na.rm = TRUE),
      sum(combined_results$Movement_pvalue[high_to_low] < significance_level, na.rm = TRUE)
    ),
    Significant_Velocity = c(
      sum(combined_results$Velocity_pvalue[low_to_high] < significance_level, na.rm = TRUE),
      sum(combined_results$Velocity_pvalue[similar] < significance_level, na.rm = TRUE),
      sum(combined_results$Velocity_pvalue[high_to_low] < significance_level, na.rm = TRUE)
    )
  )
  
  # Calculate percentages
  summary_table$Movement_Percent <- round(summary_table$Significant_Movement / summary_table$Total_Count * 100, 1)
  summary_table$Velocity_Percent <- round(summary_table$Significant_Velocity / summary_table$Total_Count * 100, 1)
  
  # Display detailed results
  cat("\nLow-to-high magnitude sequence results:\n")
  cat("  Significant movement p-values:", summary_table$Significant_Movement[1], 
      "out of", summary_table$Total_Count[1], 
      sprintf("(%.1f%%)\n", summary_table$Movement_Percent[1]))
  cat("  Significant velocity p-values:", summary_table$Significant_Velocity[1], 
      "out of", summary_table$Total_Count[1], 
      sprintf("(%.1f%%)\n", summary_table$Velocity_Percent[1]))
  
  cat("\nSimilar magnitude sequence results:\n")
  cat("  Significant movement p-values:", summary_table$Significant_Movement[2], 
      "out of", summary_table$Total_Count[2], 
      sprintf("(%.1f%%)\n", summary_table$Movement_Percent[2]))
  cat("  Significant velocity p-values:", summary_table$Significant_Velocity[2], 
      "out of", summary_table$Total_Count[2], 
      sprintf("(%.1f%%)\n", summary_table$Velocity_Percent[2]))
  
  cat("\nHigh-to-low magnitude sequence results:\n")
  cat("  Significant movement p-values:", summary_table$Significant_Movement[3], 
      "out of", summary_table$Total_Count[3], 
      sprintf("(%.1f%%)\n", summary_table$Movement_Percent[3]))
  cat("  Significant velocity p-values:", summary_table$Significant_Velocity[3], 
      "out of", summary_table$Total_Count[3], 
      sprintf("(%.1f%%)\n", summary_table$Velocity_Percent[3]))
  
  # Display summary table
  cat("\nSummary table:\n")
  print(summary_table)
  
  # Return summary invisibly
  invisible(summary_table)
}

#' Create a publication-ready table of event comparison results
#'
#' @param results List containing movement_comparison and velocity_comparison data frames
#' @param significance_level P-value threshold for significance (default: 0.05)
#' @return A formatted data frame suitable for publication
#' @export
create_publication_table <- function(results, significance_level = 0.05) {
  # Extract data frames
  movement_df <- results$movement_comparison
  velocity_df <- results$velocity_comparison
  
  # Round p-values to two decimals
  movement_df$p_val <- round(movement_df$p_val, 2)
  velocity_df$p_val <- round(velocity_df$p_val, 2)
  
  # Create base table with movement data (these shouldn't have NAs)
  pub_table <- data.frame(
    Event_Pair = paste0(movement_df$Event, "-", movement_df$Event + 1),
    Flow_Relation = movement_df$Flow_PC1_vs_Flow_PC1_next,
    Movement_Prob_N = paste0(round(movement_df$Moved / movement_df$Num_obs * 100, 1), "%"),
    Movement_Prob_N1 = paste0(round(movement_df$Moved_next / movement_df$Num_obs_next * 100, 1), "%"),
    Movement_pvalue = ifelse(movement_df$p_val < significance_level, 
                             paste0("**", movement_df$p_val, "**"), 
                             as.character(movement_df$p_val))
  )
  
  # Add velocity data, allowing NAs in these columns
  pub_table$Velocity_N <- NA
  pub_table$Velocity_N1 <- NA
  pub_table$Velocity_pvalue <- NA
  
  # Match events between dataframes
  for(i in 1:nrow(pub_table)) {
    event_n <- movement_df$Event[i]
    
    # Find matching row in velocity dataframe
    vel_idx <- which(velocity_df$Event == event_n)
    
    if(length(vel_idx) > 0) {
      # Add velocity data where available
      if(!is.na(velocity_df$Mean_vel[vel_idx])) {
        pub_table$Velocity_N[i] <- paste0(
          round(velocity_df$Mean_vel[vel_idx], 2), 
          " ± ", 
          round(velocity_df$SD_vel[vel_idx], 2)
        )
      }
      
      if(!is.na(velocity_df$Mean_vel_next[vel_idx])) {
        pub_table$Velocity_N1[i] <- paste0(
          round(velocity_df$Mean_vel_next[vel_idx], 2), 
          " ± ", 
          round(velocity_df$SD_vel_next[vel_idx], 2)
        )
      }
      
      if(!is.na(velocity_df$p_val[vel_idx])) {
        pub_table$Velocity_pvalue[i] <- ifelse(
          velocity_df$p_val[vel_idx] < significance_level,
          paste0("**", velocity_df$p_val[vel_idx], "**"),
          as.character(velocity_df$p_val[vel_idx])
        )
      }
    }
  }
  
  # For completely missing velocity data, replace with "-"
  pub_table$Velocity_N[is.na(pub_table$Velocity_N)] <- "-"
  pub_table$Velocity_N1[is.na(pub_table$Velocity_N1)] <- "-"
  pub_table$Velocity_pvalue[is.na(pub_table$Velocity_pvalue)] <- "-"
  
  return(pub_table)
}

#' Save comparison results to Excel file with formatting
#'
#' @param results List containing movement_comparison and velocity_comparison data frames
#' @param file_path Path to save the Excel file
#' @param significance_level P-value threshold for significance (default: 0.05)
#' @return NULL (invisibly)
#' @export
save_comparison_results <- function(results, file_path, significance_level = 0.05) {
  # Create workbook
  wb <- openxlsx::createWorkbook()
  
  # Add movement comparison sheet
  openxlsx::addWorksheet(wb, "Movement_Comparison")
  
  # Process movement data
  movement_df <- results$movement_comparison %>% tidyr::drop_na()
  movement_df$p_val <- round(movement_df$p_val, 2)
  
  # Write data
  openxlsx::writeData(wb, "Movement_Comparison", movement_df)
  
  # Create style for significant p-values
  sig_style <- openxlsx::createStyle(textDecoration = "bold")
  
  # Apply conditional formatting for significant p-values
  sig_rows <- which(movement_df$p_val < significance_level)
  if(length(sig_rows) > 0) {
    openxlsx::addStyle(wb, "Movement_Comparison", style = sig_style, 
                       rows = sig_rows + 1, cols = which(colnames(movement_df) == "p_val"))
  }
  
  # Add velocity comparison sheet
  openxlsx::addWorksheet(wb, "Velocity_Comparison")
  
  # Process velocity data
  velocity_df <- results$velocity_comparison %>% tidyr::drop_na()
  velocity_df$p_val <- round(velocity_df$p_val, 2)
  
  # Write data
  openxlsx::writeData(wb, "Velocity_Comparison", velocity_df)
  
  # Apply conditional formatting for significant p-values
  sig_rows <- which(velocity_df$p_val < significance_level)
  if(length(sig_rows) > 0) {
    openxlsx::addStyle(wb, "Velocity_Comparison", style = sig_style, 
                       rows = sig_rows + 1, cols = which(colnames(velocity_df) == "p_val"))
  }
  
  # Add publication-ready table
  pub_table <- create_publication_table(results, significance_level)
  openxlsx::addWorksheet(wb, "Publication_Table")
  openxlsx::writeData(wb, "Publication_Table", pub_table)
  
  # Save workbook
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  
  invisible(NULL)
}

