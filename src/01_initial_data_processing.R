#' Pebble data processing and preparation
#'
#' This script (src/00_data_processing.R) performs preprocessing of pebble tracking data collected from mountain stream
#' environments, integrating data from multiple years (2016-2018). The processing includes:
#' - Loading and merging raw data from different sources and years
#' - Cleaning and standardizing measurements
#' - Calculating derived metrics (e.g., sphericity, nominal diameter)
#' - Computing graph-based distances for pebble movements
#' - Processing weather/flow data and identifying event characteristics
#' - Generating spatial grids and center line data for spatial analysis
#' - Creating triangulation networks for accurate distance calculations
#' - Preparing datasets for subsequent statistical modeling
#'
#' The output consists of multiple datasets saved as R objects (.RDS files)
#' in the data/preprocessed/ directory.
#'
#' @dependencies Dependencies are loaded from src/00_dependencies.R
#' @output preprocessed_observations.RDS - Base pebble and weather observations
#'         baseline_weather_df.RDS - Processed weather event data
#'         all_datasets.RDS - Comprehensive data list for modeling

# Load dependencies and utility functions
source("src/00_dependencies.R")
source("src/utils/utils.R")

#--------------------------------------------------------
# STEP 1: LOAD RAW DATA
#--------------------------------------------------------

# 1.1 Load pebble location data (2016-2017)
message("Loading 2016-2017 location data...")
path_to_16_17_data <- 'data/initial/Rilievo2016'
movement_16_17_files <- list.files(path_to_16_17_data)

locations_16_17 <- lapply(movement_16_17_files, function(file) {
  read_csv(paste0(path_to_16_17_data, "/", file))
}) %>%
  bind_rows()

# 1.2 Load pebble characteristics data (2016-2017)
message("Loading 2016-2017 characteristics data...")
characteristics_columns <- c('IDREF', 'Weight', 'a_axis', 'b_axis', 'c_axis')
characteristics_16_17 <- read.xlsx('data/initial/Pebbles_16_17.xlsx') %>%
  select(all_of(characteristics_columns))
characteristics_16_17[, -1] %<>% mutate_if(is.character, as.numeric)

# 1.3 Load pebble location data (2018)
message("Loading 2018 location data...")
locations_18 <- read.xlsx('data/initial/MERGED2018_XLS.xlsx')

# Set coordinate reference system
CRS.n <- CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

# Convert coordinates
locations_18[, c("X_Start", "Y_Start")] <-
  data.frame(project(
    locations_18[, c("X_Start", "Y_Start")],
    CRS.n,
    inverse = TRUE
  ))
locations_18[, c("X_End", "Y_End")] <-
  data.frame(project(
    locations_18[, c("X_End", "Y_End")],
    CRS.n,
    inverse = TRUE
  ))

# 1.4 Load pebble characteristics data (2018)
message("Loading 2018 characteristics data...")
characteristics_18 <- read.xlsx(
  'data/initial/DB 2018.xlsx',
  sheet = 'Foglio1'
) %>%
  set_colnames(characteristics_columns) %>%
  mutate(
    a_axis = ifelse(a_axis > 200, 200, a_axis), # Fix outliers
    b_axis = ifelse(b_axis > 135, 135, b_axis)
  )

# 1.5 Load weather data
message("Loading weather data...")
df_weather <- preprocess_weather_data("data/initial/Dati idrometeo.xlsx")

#--------------------------------------------------------
# STEP 2: CLEAN AND MERGE DATA
#--------------------------------------------------------

# 2.1 Combine and filter location data
message("Combining location data...")
df_loc <- bind_rows(locations_16_17, locations_18)
df_loc[, 9:12] <- df_loc[, 9:12] %>% round(., 8)
df_ch <- bind_rows(characteristics_16_17, characteristics_18)

# Find common pebble IDs between location and characteristics data
good_IDREF <- intersect(df_loc$IDREF, df_ch$IDREF)

# Filter for valid records and coordinates
df_loc %<>%
  filter(
    IDREF %in% good_IDREF,
    Y_Start >= 45,
    X_Start >= 9,
    Y_End >= 45,
    X_End >= 9
  )
df_ch %<>% filter(IDREF %in% good_IDREF)

# 2.2 Process pebble characteristics
message("Processing pebble characteristics...")
# Ensure axes are ordered properly (a_axis >= b_axis >= c_axis)
df_ch[, c("a_axis", "b_axis", "c_axis")] <-
  sort_axes_by_row(df_ch[, c("a_axis", "b_axis", "c_axis")])

# Calculate derived metrics
df_ch %<>%
  mutate(
    Nominal_diameter = (a_axis * b_axis * c_axis)^(1 / 3),
    Elongation = b_axis / a_axis,
    Platyness = c_axis / a_axis,
    Sphericity = (b_axis * c_axis / a_axis^2)^(1 / 3)
  )

# Calculate PCA components
pca_result <- princomp(df_ch %>% select(a_axis, b_axis, c_axis))
df_ch %<>%
  bind_cols(
    data.frame(pca_result$scores[, 1:2]) %>%
      set_colnames(c('pebble_PC1', 'pebble_PC2'))
  )

# 2.3 Process location data - add missing events
message("Processing location data...")
df_loc <- process_loc_data(df_loc)

#--------------------------------------------------------
# STEP 3: CALCULATE SPATIAL METRICS
#--------------------------------------------------------

# 3.1 Load pool shape
message("Calculating graph distances...")
pools_shape <- load_pool_shape(
  "data/initial/shp",
  "MORFOLOGIA_OKSHP_WGS84"
)

# 3.1 Calculate graph distances
message("Calculating graph distances...")
distances <- df_loc$Distance_m
df_loc$graph_dist <- calculate_graph_dist(
  df_loc[, c("X_Start", "Y_Start", "X_End", "Y_End")],
  df_loc$Distance_m,
  pools_shape
)

# 3.2 Determine polygon ID for each pebble start location
message("Determining polygon IDs...")
df_loc$polygon_id <- get_polygon_id(
  df_loc[, c("X_Start", "Y_Start")],
  pools_shape
)

# 3.3 Add stuck flags
message("Adding movement status flags...")
df_loc %<>%
  group_by(IDREF) %>%
  arrange(EventoStart) %>%
  mutate(
    is_stuck_1 = (lag(EventoEnd, default = -1) == EventoStart) &
      (lag(Distance_m, default = -1) == 0),
    is_stuck_2 = (lag(EventoEnd, 2, default = -1) + 1 == EventoStart) &
      (lag(Distance_m, default = -1) == 0) &
      (lag(Distance_m, 2, default = -1) == 0),
    is_stuck_3 = (lag(EventoEnd, 3, default = -1) + 2 == EventoStart) &
      (lag(Distance_m, default = -1) == 0) &
      (lag(Distance_m, 2, default = -1) == 0) &
      (lag(Distance_m, 3, default = -1) == 0)
  ) %>%
  ungroup()

# 3.4 Standardize location types
message("Standardizing location types...")
df_loc$TipoStart[df_loc$polygon_id == 3] <- 'Run_rapid'
df_loc$TipoStart[df_loc$polygon_id == 7] <- 'Planebed'
df_loc$TipoStart[df_loc$polygon_id == 23] <- 'Cascade'
df_loc$TipoStart[df_loc$polygon_id == 39] <- 'Steep_poolzone'
df_loc$TipoStart[df_loc$polygon_id == 9] <- 'Banks'


df_loc_ <- df_loc
# 3.5 Filter backward movements
message("Filtering backward movements...")
df_loc %<>%
  filter(
    X_Start >= X_End |
      Y_Start >= Y_End |
      graph_dist <= 1
  ) %>%
  filter(EventoStart == EventoEnd - 1)


#--------------------------------------------------------
# STEP 4: PROCESS WEATHER DATA
#--------------------------------------------------------

message("Processing weather data...")
water_depth_threshold <- 35

df_weather_processed <- df_weather %>%
  group_by(EventoStart) %>%
  mutate(
    id = row_number(),
    is_thr = h_CP_cm_ >= water_depth_threshold
  ) %>%
  filter(is_thr == TRUE) %>%
  mutate(
    id_delta = id - lag(id, default = NA),
    id_delta_diff = id_delta > 1,
    number_of_subevents = sum(id_delta_diff, na.rm = TRUE) + 1
  ) %>%
  ungroup()

# Calculate event summaries
baseline_weather <- df_weather_processed %>%
  group_by(EventoStart) %>%
  summarise(
    duration = n(),
    mean_h = mean(h_CP_cm_),
    n_o_s = as.numeric(mean(number_of_subevents) > 1) + 1,
    mean_Q = mean(Q_CP_m3_s_)
  ) %>%
  mutate(id = row_number())

# Cluster events
scaled_data <- baseline_weather %>%
  select(duration, mean_h, mean_Q) %>%
  mutate(across(everything(), rescale))

dist_matrix <- dist(scaled_data)
hcl_result <- hclust(dist_matrix, method = 'ward.D2')
clusters <- cutree(hcl_result, 3)

# Calculate PCA components
pca_result <- princomp(scaled_data)

# Add PCA and clusters to the dataset
baseline_weather %<>%
  bind_cols(
    data.frame(pca_result$scores[, 1:2]) %>%
      set_colnames(c('weather_PC1', 'weather_PC2'))
  ) %>%
  mutate(cl_ = clusters) %>%
  select(-id) # Remove temporary ID column


#--------------------------------------------------------
# STEP 5: GENERATE POOL GRID AND TEST DATA
#--------------------------------------------------------

message("Generating pool grid and test data...")

# 5.1 Load pool shape data using sf
pools_shape <- load_pool_shape(
  dsn = "data/initial/shp",
  layer = "MORFOLOGIA_OKSHP_WGS84"
)

# 5.2 Generate pool grid
grid_data <- generate_pool_grid(pools_shape, total_points_num = 3000, CRS.n)
pool_grid <- grid_data$pool_grid
pools_shape_diff <- grid_data$pools_shape_diff

# 5.3 Create derived datasets
message("Creating derived datasets...")

# Create the combined model dataframe
df_models_ <- df_ch %>%
  left_join(df_loc) %>%
  mutate(Event = EventoEnd) %>%
  left_join(baseline_weather %>% rename(Event = EventoStart)) %>%
  drop_na()

# Add dummy variables for location type and calculate velocity
df_models_ <- df_models_ %>%
  bind_cols(
    model.matrix(~ TipoStart - 1, data = .) %>%
      as.data.frame()
  ) %>%
  mutate(
    graph_vel = graph_dist / duration,
    id_o = 1:nrow(.)
  )

1 +
  1
# Create features dataset
df_ <- df_models_ %>%
  dplyr::select(
    -TipoEnd,
    -TipoStart,
    -X_Start,
    -Y_Start,
    -X_End,
    -Y_End,
    -ID,
    -IDREF,
    -IDNUM,
    -EventoStart,
    -Event,
    -EventoEnd,
    -duration,
    -graph_vel,
    -Distance_m,
    -starts_with('Tipo'),
    -polygon_id,
    -graph_dist
  )

# Create target data dataset
df_pred_ <- df_models_ %>%
  dplyr::select(duration, graph_dist, graph_vel, id_o)

# Create locations dataset with coordinates
df_coordinates_and_type_ <- df_models_ %>%
  dplyr::select(
    TipoEnd,
    TipoStart,
    X_Start,
    Y_Start,
    X_End,
    Y_End,
    EventoStart,
    Event,
    polygon_id,
    EventoEnd,
    id_o,
    starts_with('Tipo')
  )

# Add projected coordinates
df_coordinates_and_type_[, c('X_Start_WGS84', 'Y_Start_WGS84')] <-
  proj4::project(df_coordinates_and_type_[, c("X_Start", "Y_Start")], CRS.n)

# Calculate min/max for scaling
x_min <- min(df_coordinates_and_type_$X_Start_WGS84)
x_max <- max(df_coordinates_and_type_$X_Start_WGS84)
y_min <- min(df_coordinates_and_type_$Y_Start_WGS84)
y_max <- max(df_coordinates_and_type_$Y_Start_WGS84)

# Scale coordinates
df_coordinates_and_type_$X_Start_WGS84_scaled <-
  scale_to_unit_range(df_coordinates_and_type_$X_Start_WGS84, x_min, x_max)
df_coordinates_and_type_$Y_Start_WGS84_scaled <-
  scale_to_unit_range(df_coordinates_and_type_$Y_Start_WGS84, y_min, y_max)

# Scale pool grid coordinates
pool_grid$x_WGS84_scaled <- scale_to_unit_range(pool_grid$x_WGS84, x_min, x_max)
pool_grid$y_WGS84_scaled <- scale_to_unit_range(pool_grid$y_WGS84, y_min, y_max)

# 5.4 Create test pebble data
# Define test cases with varied characteristics
id_0_case <- c(1089, 1025, 1159, 1679, 1640, 1238)

# Create test pebble dataset
df_test_peb_ <- df_[id_0_case, ] %>%
  select(-colnames(baseline_weather)[3:7], -id_o) %>%
  bind_cols(
    baseline_weather %>%
      filter(EventoStart == 15) %>%
      .[rep(1, 6), 2:7]
  )


#--------------------------------------------------------
# STEP 6: GENERATE CENTER LINE DATA
#--------------------------------------------------------

message("Generating center line data...")

# Generate center line data
step_df <- generate_center_line_data(
  "data/initial/Line.kml",
  pools_shape,
  CRS.n
)


#--------------------------------------------------------
# STEP 7: CREATE COMPREHENSIVE DATA LIST
#--------------------------------------------------------

message("Creating comprehensive data list...")

# Create triangulation distance matrix if not already calculated
if (!exists("triangl_dist_btw_points")) {
  message("Calculating triangulation distances...")

  # Extract pool border using utility function
  pool_border <- extract_pool_border(pools_shape)

  # Extract pebble start locations
  pebbles_loc_st <- df_coordinates_and_type_ %>%
    select('X_Start', 'Y_Start') %>%
    dplyr::distinct() %>%
    set_colnames(c('X', 'Y')) %>%
    anti_join(pool_border)

  # Create triangulation coordinates using existing utility functions approach
  coords_for_triang <- unique(rbind(
    as.matrix(pool_border),
    as.matrix(pebbles_loc_st)
  ))

  # Calculate triangulation distances using utility function
  result <- triangulate_and_calculate_distance(coords_for_triang, pool_border)
  triangl_dist <- result$distances

  # Extract distances between pebble locations
  st_ <- nrow(pool_border) + 1
  en_ <- nrow(coords_for_triang)
  triangl_dist_btw_points <- triangl_dist[st_:en_, st_:en_]

  # Cleanup temporary objects
  rm(triangl_dist, st_, en_)
}

# Create data list
data_list <- list(
  dfs = list(
    df_ = df_,
    df_pred_ = df_pred_,
    df_coordinates_and_type_ = df_coordinates_and_type_
  ),
  dfs_test = list(
    df_test_peb_ = df_test_peb_,
    pool_grid = pool_grid,
    step_df = step_df
  ),
  shape = list(
    pools_shape = pools_shape,
    pools_shape_WGS84 = pools_shape_diff
  ),
  triangl_dist_btw_points = triangl_dist_btw_points
)
#--------------------------------------------------------
# STEP 8: SAVE PROCESSED DATA
#--------------------------------------------------------

message("Saving processed data...")

# Save base observations
saveRDS(
  list(
    df_ch = df_ch,
    df_loc = df_loc,
    df_weather = df_weather
  ),
  file = 'data/preprocessed/preprocessed_observations.RDS',
  compress = FALSE
)

# Save weather data
saveRDS(
  baseline_weather,
  file = 'data/preprocessed/baseline_weather_df.RDS',
  compress = FALSE
)

# Save comprehensive data list
saveRDS(
  data_list,
  file = 'data/preprocessed/all_datasets.RDS',
  compress = FALSE
)

message("Data processing complete!")
