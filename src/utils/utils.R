#' Sort pebble axes within each row in descending order
#'
#' This function sorts all values within each row of a data frame in descending order.
#' Used primarily to ensure pebble axis measurements follow a_axis >= b_axis >= c_axis.
#'
#' @param df A data frame containing axes values to be sorted within each row
#' @return A data frame with values sorted within each row
sort_axes_by_row <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      sorted = list(sort(c_across(everything()), decreasing = TRUE))
    ) %>%
    mutate(
      a_axis = sorted[[1]],
      b_axis = sorted[[2]],
      c_axis = sorted[[3]]
    ) %>%
    select(-sorted) %>%
    ungroup()
}

#' Scale numeric values to unit range [0,1]
#'
#' @param x Numeric vector to scale
#' @param x_min Minimum value in the original scale
#' @param x_max Maximum value in the original scale
#' @return Scaled numeric vector with values between 0 and 1
scale_to_unit_range <- function(x, x_min, x_max) {
  (x - x_min) / (x_max - x_min)
}


#' Generate a sequence of event pairs
#'
#' @param start The first event start ID
#' @param end The last event end ID
#' @return A data frame with EventoStart and EventoEnd columns
generate_event_sequence <- function(start, end) {
  data.frame(
    EventoStart = start:(end - 1),
    EventoEnd = (start + 1):end
  )
}

#' Expand event data to include all intermediate events
#'
#' Takes a data frame with EventoStart and EventoEnd columns and expands it to include
#' all intermediate events in separate rows.
#'
#' @param df A data frame containing at least EventoStart and EventoEnd columns
#' @return A data frame with expanded rows for all intermediate events
expand_events <- function(df) {
  df %>%
    rowwise() %>%
    # Generate the event sequence for each row
    mutate(event_sequence = list(generate_event_sequence(EventoStart, EventoEnd))) %>%
    # Unnest the event_sequence column
    unnest(event_sequence, names_sep = '_') %>%
    # Update the original EventoStart and EventoEnd columns
    mutate(
      EventoStart = event_sequence_EventoStart,
      EventoEnd = event_sequence_EventoEnd
    ) %>%
    # Select only the required columns
    select(names(df), -starts_with("event_sequence"))
}

#' Process location data to fill missing intermediate events
#'
#' Takes a data frame of location records and fills in missing events for cases where
#' a pebble didn't move (distance = 0) across multiple event periods.
#'
#' @param df_loc Data frame of pebble location records
#' @return Processed data frame with expanded event sequences
process_loc_data <- function(df_loc) {
  # Filter for cases with non-consecutive events and zero movement
  df_loc_temp <- df_loc %>%
    filter(EventoEnd - EventoStart >= 2, Distance_m == 0)
  
  # Expand events for these cases
  df_loc_expanded <- expand_events(df_loc_temp)
  
  # Combine with the rest of the dataset
  result <- df_loc %>%
    filter(!(EventoEnd - EventoStart >= 2 & Distance_m == 0)) %>%
    bind_rows(df_loc_expanded)
  
  return(result)
}

#' Read pool shape data
#'
#' @param dsn Data source name (directory path)
#' @param layer Layer name
#' @return An sf object containing the pool shape data
load_pool_shape <- function(dsn = "data/initial/shp", 
                             layer = "MORFOLOGIA_OKSHP_WGS84") {
  # Read shape data directly as sf object
  pools_shape <- sf::st_read(dsn = dsn, layer = layer)
  
  # Fix invalid geometries
  pools_shape <- sf::st_make_valid(pools_shape)
  
  return(pools_shape)
}

#' Extract coordinates of the exterior border of unioned polygons
#'
#' This function takes an sf object containing multiple polygons, performs a robust union
#' to find the main outline, extracts the exterior boundary of that outline,
#' and returns its unique coordinates.
#'
#' @param pools_shape An sf object (e.g., from `sf::st_read`) containing
#'   `sfc_POLYGON` or `sfc_MULTIPOLYGON` geometries to be unioned.
#' @param precision The number of decimal places to round the final
#'   coordinates to (default is 8).
#'
#' @return A data frame with columns `X` and `Y` containing the unique,
#'   rounded coordinates of the exterior border of the unioned shape.
extract_pool_border <- function(pools_shape, precision = 8) {
  
  if (!inherits(pools_shape, "sf")) {
    stop("Input 'pools_shape' must be an sf object.")
  }
  unioned_geom <- sf::st_union(pools_shape)
  
  cleaned_geom <- sf::st_buffer(unioned_geom, dist = 0)
  
  polygons_cast <- sf::st_cast(cleaned_geom, "POLYGON")
  
  if (length(polygons_cast) == 0) {
    stop("No valid polygon geometry resulted from the union and cleaning process.")
  } else if (length(polygons_cast) > 1) {
    # If multiple disjoint polygons result, pick the one with the largest area
    message(paste("Union resulted in", length(polygons_cast),
                  "distinct polygons. Selecting the largest by area for the border."))
    main_polygon <- polygons_cast[which.max(sf::st_area(polygons_cast))]
  } else {
    # Only one polygon resulted
    main_polygon <- polygons_cast
  }
  
  pool_exterior_ring_geometry <- sf::st_exterior_ring(main_polygon)
  
  pool_border_coords <- sf::st_coordinates(pool_exterior_ring_geometry) %>%
    as.data.frame() %>%
    dplyr::select(X, Y) %>%
    dplyr::mutate(across(c(X, Y), ~round(., digits = precision))) %>%
    dplyr::distinct()
  
  return(pool_border_coords)
}

#' Extract pool intra coordinates
#'
#' @param pools_shape The pool shape sf object
#' @param pool_border Data frame of pool border coordinates
#' @return A data frame of interior coordinates
extract_pool_intra <- function(pools_shape, pool_border) {
  # Extract all coordinates from all geometries
  pool_intra <- sf::st_coordinates(pools_shape) %>%
    as.data.frame() %>%
    # Select only X and Y columns
    dplyr::select(X, Y) %>%
    # Round coordinates to 8 decimal places and remove duplicates
    dplyr::mutate(across(c(X, Y), ~round(., 8))) %>%
    dplyr::distinct() %>%
    # Rename columns to match pool_border
    dplyr::rename(x = X, y = Y) %>%
    # Remove coordinates that are on the border
    dplyr::anti_join(pool_border, by = c("x" = "X", "y" = "Y"))
  
  return(pool_intra)
}

#' Prepare coordinates for triangulation
#'
#' This function merges pool border, intra-pool, and pebble coordinates,
#' removing duplicates (excluding the pool_border) and ensuring that
#' pebble start/end points not already part of the triangulation set are added.
#'
#' @param pool_border Matrix or data frame of pool boundary coordinates (x, y)
#' @param pool_intra Matrix or data frame of intra-pool structure coordinates (x, y)
#' @param pebble_coords Data frame containing columns X_Start, Y_Start, X_End, Y_End for pebble movements
#' @return Matrix of unique coordinates to use in triangulation
#' @export
build_coords_for_triang <- function(pool_border, pool_intra, pebble_coords) {
  # Ensure correct names for consistency
  if (!all(c("x", "y") %in% colnames(pool_border))) {
    colnames(pool_border) <- c("x", "y")
  }
  if (!all(c("x", "y") %in% colnames(pool_intra))) {
    colnames(pool_intra) <- c("x", "y")
  }
  
  # Start with pool_border as-is (preserve order)
  coords <- rbind(pool_border, pool_intra) %>%
    dplyr::mutate(x = round(x, 8), y = round(y, 8)) %>%
    dplyr::distinct()
  
  pebbles_start <- pebble_coords %>%
    dplyr::select(X_Start, Y_Start) %>%
    dplyr::rename(x = X_Start, y = Y_Start) %>%
    dplyr::mutate(x = round(x, 8), y = round(y, 8)) %>%
    dplyr::distinct() %>%
    dplyr::anti_join(coords, by = c("x", "y"))
  
  pebbles_end <- pebble_coords %>%
    dplyr::select(X_End, Y_End) %>%
    dplyr::rename(x = X_End, y = Y_End) %>%
    dplyr::mutate(x = round(x, 8), y = round(y, 8)) %>%
    dplyr::distinct() %>%
    dplyr::anti_join(coords, by = c("x", "y"))
  
  # Combine with original (non-rounded) pool_border to preserve geometry
  final_coords <- rbind(
    coords,  # untouched
    rbind(pebbles_start, pebbles_end) %>%
      dplyr::mutate(x = round(x, 8), y = round(y, 8)) %>%
      dplyr::distinct()
  ) %>% 
    as.matrix()
  
  return(final_coords)
}

#' Triangulate domain and calculate distance on the graph
#'
#' @param coords_for_triang Coordinates for triangulation
#' @param pool_border Pool border coordinates
#' @return A matrix of distances between points
triangulate_and_calculate_distance<- function(coords_for_triang, pool_border) {
  # Ensure column names are consistent
  if (!all(c("X", "Y") %in% colnames(pool_border))) {
    colnames(pool_border) <- c("X", "Y")
  }
  
  # Create boundary segments for triangulation
  m_border <- cbind(1:nrow(pool_border), c(2:nrow(pool_border), 1))
  
  # Create PSLG (Planar Straight Line Graph) for triangulation
  data_plsg <- RTriangle::pslg(P = coords_for_triang, S = m_border)
  
  # Triangulate the domain
  data_traingl <- RTriangle::triangulate(data_plsg, a = 1e-11)
  
  # Extract edge list
  edge_list <- data_traingl$E
  
  # Create graph from edge list
  graph <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
  
  # Calculate edge weights using Haversine distance (for geographic coordinates)
  edge_weights <- apply(edge_list, 1, function(i) {
    geosphere::distHaversine(data_traingl$P[i[1], 1:2],
                             data_traingl$P[i[2], 1:2])
    
  })
  
  # Calculate shortest paths between all pairs of points
  triangl_dist <- igraph::distances(
    graph, 
    v = igraph::V(graph),
    to = igraph::V(graph),
    mode = c("all"),
    weights = edge_weights, 
    algorithm = "automatic"
  )
  
  return(list(
    distances = triangl_dist,
    triangulation = data_traingl
  ))
}

#' Calculate graph distances for pebble movements
#'
#' @param pebble_coords Data frame with pebble coordinates
#' @param distances Vector of straight-line distances
#' @param pools_shape Pools shape
#' @return Vector of graph distances
calculate_graph_dist <- function(pebble_coords, distances, pools_shape) {

  pool_border <- extract_pool_border(pools_shape)
  pool_intra <- extract_pool_intra(pools_shape, pool_border)
  
  coords_for_triang <- build_coords_for_triang(pool_border, pool_intra, pebble_coords)
  
  result <- triangulate_and_calculate_distance(coords_for_triang, pool_border)
  triangl_dist <- result$distances
  data_traingl <- result$triangulation
  
  # Match pebble coordinates to vertices in the triangulation more efficiently
  pebble_ids <- purrr::map_dfr(1:nrow(pebble_coords), function(row_idx) {
    row <- pebble_coords[row_idx, ]
    
    start_id <- which(data_traingl$P[, 1] == row$X_Start & 
                        data_traingl$P[, 2] == row$Y_Start)[1]
    
    end_id <- which(data_traingl$P[, 1] == row$X_End & 
                      data_traingl$P[, 2] == row$Y_End)[1]
    
    # Return as a one-row tibble
    tibble::tibble(
      row = row_idx,
      start_id = start_id,
      end_id = end_id
    )
  })
  
  # Extract distances using the indices
  graph_dist <- sapply(1:nrow(pebble_ids), function(i) {
    triangl_dist[pebble_ids$start_id[i], pebble_ids$end_id[i]]
  })
  
  # Replace infinite distances with original straight-line distances
  graph_dist[is.infinite(graph_dist)] <- distances[is.infinite(graph_dist)]
  
  return(graph_dist)
}

#' Determine polygon ID for each pebble location
#'
#' @param pebble_points Data frame with pebble point coordinates
#' @param pools_shape Pools shape
#' @return Vector of polygon IDs
get_polygon_id <- function(pebble_points, pools_shape) {
  # Convert pebble points to sf object
  if (!("X" %in% colnames(pebble_points) && "Y" %in% colnames(pebble_points))) {
    # If columns are not named X and Y, rename first two columns
    names(pebble_points)[1:2] <- c("X", "Y")
  }
  
  # Create sf points from coordinates
  pebble_sf <- sf::st_as_sf(
    pebble_points, 
    coords = c("X", "Y"),
    crs = sf::st_crs(pools_shape)
  )
  
  # For each point, determine which polygon contains it
  polygon_id <- sapply(1:nrow(pebble_sf), function(i) {
    point <- pebble_sf[i, ]
    # st_intersects returns a list of indices of geometries that intersect with the point
    intersections <- sf::st_intersects(point, pools_shape)
    
    # Return the first polygon ID that contains the point, or NA if none
    if (length(intersections[[1]]) > 0) {
      return(intersections[[1]][1])
    } else {
      return(NA)
    }
  })
  
  return(polygon_id)
}

#' Generate pool grid based on shape data using sf
#'
#' Creates a grid of points within each pool polygon, with density proportional to area
#'
#' @param pools_shape sf object of the pool shapes
#' @param total_points_num Total number of points to generate across all pools
#' @param CRS.n Coordinate reference system for projection
#' @return A data frame with grid points and their properties
#' @export
generate_pool_grid <- function(pools_shape, total_points_num = 3000, CRS.n) {
  # Calculate area of each polygon
  area_pol <- sf::st_area(pools_shape)
  pool_l <- nrow(pools_shape)
  
  # Distribute points proportionally to area
  points_num <- round(as.numeric(area_pol) / sum(as.numeric(area_pol)) * total_points_num, 0)
  
  # Generate points for each polygon
  pool_grid <- lapply(1:pool_l, function(i) {
    # Extract polygon boundary
    pool_geom <- pools_shape[i, ]
    bbox <- sf::st_bbox(pool_geom)
    
    # Generate more points than needed, as some might fall outside the polygon
    pool_df <- data.frame(
      x = runif(3000, bbox["xmin"], bbox["xmax"]),
      y = runif(3000, bbox["ymin"], bbox["ymax"])
    )
    
    # Convert to sf points
    pts <- sf::st_as_sf(pool_df, coords = c("x", "y"), crs = sf::st_crs(pools_shape))
    
    # Filter points inside the polygon
    pts_in <- pts[sf::st_intersects(pts, pool_geom, sparse = FALSE), ]
    
    # Extract coordinates and add metadata
    coords <- sf::st_coordinates(pts_in)
    result <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      in_out = TRUE,
      poly_id = i,
      TipoStart = as.character(pools_shape$TIPO[i])
    ) %>%
      head(points_num[i])
    
    return(result)
  })
  
  # Combine all points
  pool_grid <- do.call(rbind, pool_grid)
  
  # Create dummy variables for location type
  dummy_vars <- model.matrix(~ TipoStart - 1, data = pool_grid)
  pool_grid <- cbind(pool_grid, dummy_vars)
  
  # Clean up the data
  pool_grid <- pool_grid %>% distinct() %>% na.omit()
  
  # Project coordinates
  pools_shape_diff <- sf::st_transform(pools_shape, CRS.n)
  
  # Create sf points object
  pts_sf <- sf::st_as_sf(pool_grid, coords = c("x", "y"), crs = sf::st_crs(pools_shape))
  pts_transformed <- sf::st_transform(pts_sf, CRS.n)
  coords_transformed <- sf::st_coordinates(pts_transformed)
  
  pool_grid[, c("x_WGS84", "y_WGS84")] <- coords_transformed
  
  return(list(
    pool_grid = pool_grid,
    pools_shape_diff = pools_shape_diff
  ))
}

#' Generate center line data using sf
#'
#' Creates a dense set of points along the center line for functional prediction
#'
#' @param line_path Path to the center line KML file
#' @param pools_shape sf object of the pool shapes
#' @param CRS.n Coordinate reference system for projection
#' @return A data frame with points along the center line
#' @export
generate_center_line_data <- function(line_path = "data/initial/Line.kml", 
                                      pools_shape, 
                                      CRS.n) {
  # Load center line as sf object
  center_line_sf <- sf::st_read(line_path, quiet = TRUE)
  
  # Extract linestring coordinates
  line_coords <- sf::st_coordinates(center_line_sf)
  center_line <- data.frame(x = line_coords[, 1], y = line_coords[, 2])
  
  nr_ <- nrow(center_line) - 1
  
  # Calculate distances between line points
  dist_between_line_points <- sapply(1:nr_, function(i) {
    sf::st_length(sf::st_linestring(matrix(
      c(center_line$x[i], center_line$y[i],
        center_line$x[i+1], center_line$y[i+1]), 
      ncol = 2, byrow = TRUE
    )))
  })
  
  # Generate step points along center line
  step_list <- list()
  for (i in 1:(nrow(center_line) - 1)) {
    dist_r <- 10 * (round(as.numeric(dist_between_line_points[i]), 1))
    x_s <- center_line$x[i]
    x_e <- center_line$x[i + 1]
    y_s <- center_line$y[i]
    y_e <- center_line$y[i + 1]
    
    step_list[[i]] <- sapply(
      seq(0, 1, length.out = dist_r + 1), 
      function(step) {
        return(c(
          x_s + (x_e - x_s) * step,
          y_s + (y_e - y_s) * step
        ))
      }
    ) %>% t()
  }
  
  # Combine step points
  step_df <- do.call('rbind', step_list) %>% 
    data.frame() %>% 
    distinct() %>%
    set_colnames(c('x', 'y'))
  
  # Assign location types to step points (using sf)
  step_points <- sf::st_as_sf(step_df, coords = c("x", "y"), crs = sf::st_crs(pools_shape))
  
  # Initialize TipoStart column
  step_df$TipoStart <- "0"
  
  # For each polygon, identify which points fall inside it
  for (i in 1:nrow(pools_shape)) {
    polygon <- pools_shape[i, ]
    intersects <- sf::st_intersects(step_points, polygon, sparse = FALSE)
    step_df$TipoStart[intersects] <- as.character(pools_shape$TIPO[i])
  }
  
  # Project step points
  step_points_transformed <- sf::st_transform(step_points, CRS.n)
  coords_transformed <- sf::st_coordinates(step_points_transformed)
  step_df[, c('x_WGS84', 'y_WGS84')] <- coords_transformed
  
  return(step_df)
}

#' Process weather data
#'
#' @param file_path Path to the Excel file containing weather data
#' @param start_sheet First sheet index to read
#' @param end_sheet Last sheet index to read
#' @return A data frame of processed weather data
preprocess_weather_data <- function(file_path, start_sheet = 3, end_sheet = 30) {
  # Generate the sheet indices
  sheet_indices <- start_sheet:end_sheet
  
  # Read and preprocess the data from multiple sheets
  df_weather <- purrr::map_dfr(sheet_indices, function(i) {
    readxl::read_excel(
      path = file_path,
      sheet = i,
      skip = 2,
      col_types = "guess"
    ) %>%
      select(2:5) %>%
      set_names(gsub("\\.+", "_", names(.))) %>%
      mutate(EventoStart = i - start_sheet + 1)
  }) %>%
    set_colnames(c("Data", "h_CP_cm_", "h_CP_m_", "Q_CP_m3_s_", "EventoStart")) %>%
    select(-h_CP_m_)
  
  return(df_weather)
}
