#' Yellow Pebbles Analysis - Dependencies
#' 
#' This script loads all required libraries for the Yellow Pebbles analysis pipeline.
#' It handles installation of missing packages, sets global options, and verifies
#' that all dependencies are properly loaded. 
#'
#' @author Yellow Pebbles Research Team
#' @note If you encounter memory issues with large datasets, adjust the java.parameters option


# Set default options
options(max.print = 100)
options(scipen = 999)  # Avoid scientific notation
options(stringsAsFactors = FALSE)
options(java.parameters = "-Xmx4096m")  # For large Excel files with openxlsx

#' Define all required packages
packages <- c(
  # Core data processing and manipulation
  "magrittr",     # For pipe operators
  "dplyr",        # Data manipulation
  "tidyr",        # Tidy data
  "readr",        # Reading CSV files
  "purrr",        # Functional programming
  
  # File I/O
  "openxlsx",     # Reading/writing Excel files
  "readxl",       # Alternative Excel reader
  
  # Statistical modeling
  "glmnet",       # For regularized regression
  "xgboost",      # For gradient boosting
  "markovchain",  # For Markov chain analysis
  "stats",        # Base R statistics
  "pROC",         # For ROC curve analysis
  
  # Geospatial analysis
  "sp",           # Spatial data handling
  "sf",           # Simple features for spatial data
  "proj4",        # For coordinate projections
  "raster",       # For raster data handling
  "geosphere",    # For distance calculations
  "igraph",       # For graph-based operations
  "RTriangle",    # For triangulation
  "mgcv",         # For point-in-polygon operations
  "gstat",        # For geostatistical modeling
  
  # Visualization and output
  "ggplot2",      # For plotting
  "plotly",       # For interactive visualizations
  "scales",       # For scaling data
  "BAMMtools",    # For Jenks Natural Breaks
  "viridis"      # For color scales
  
)

#' Install missing packages
missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

#' Load all packages
invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE))
  message("Loaded package: ", pkg)
}))

# Check if all packages are loaded correctly
loaded_packages <- search()
loaded_packages <- gsub("package:", "", loaded_packages[grepl("package:", loaded_packages)])
not_loaded <- packages[!packages %in% loaded_packages]

if(length(not_loaded) > 0) {
  warning("Some packages could not be loaded: ", paste(not_loaded, collapse = ", "))
}

message("Dependencies loaded successfully!")