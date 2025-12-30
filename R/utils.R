# Internal utility functions for spopt
# These are not exported and not documented

# Validate and clean sf data for regionalization
# Removes empty geometries (with warning) and checks for NA values
validate_regionalization_data <- function(data, check_cols, call_name = "regionalization") {

  n_original <- nrow(data)
  removed_idx <- integer(0)

  # Check for empty geometries
  empty_geom <- sf::st_is_empty(data)
  if (any(empty_geom)) {
    n_empty <- sum(empty_geom)
    removed_idx <- which(empty_geom)
    warning(
      sprintf("%s: Removed %d observation(s) with empty geometries", call_name, n_empty),
      call. = FALSE
    )
    data <- data[!empty_geom, ]
  }

  # Check for NA values in required columns
  if (length(check_cols) > 0) {
    existing_cols <- intersect(check_cols, names(data))

    if (length(existing_cols) > 0) {
      na_info <- lapply(existing_cols, function(col) {
        which(is.na(data[[col]]))
      })
      names(na_info) <- existing_cols

      cols_with_na <- existing_cols[sapply(na_info, length) > 0]

      if (length(cols_with_na) > 0) {
        na_counts <- sapply(na_info[cols_with_na], length)
        col_summary <- paste(
          sprintf("%s (%d)", cols_with_na, na_counts),
          collapse = ", "
        )
        stop(
          sprintf(
            "%s: Found NA values in columns: %s\nPlease filter these rows or impute values before running %s()",
            call_name, col_summary, call_name
          ),
          call. = FALSE
        )
      }
    }
  }

  list(data = data, removed_idx = removed_idx)
}

# Attach spopt metadata to result
attach_spopt_metadata <- function(result, metadata) {
  attr(result, "spopt") <- metadata
  result
}

# Extract attribute columns from sf object
extract_attrs <- function(data, attrs) {
  if (is.null(attrs)) {
    geom_col <- attr(data, "sf_column")
    numeric_cols <- sapply(sf::st_drop_geometry(data), is.numeric)
    attr_names <- names(numeric_cols)[numeric_cols]
  } else if (is.character(attrs)) {
    attr_names <- attrs
  } else {
    stop("`attrs` must be a character vector of column names", call. = FALSE)
  }

  if (length(attr_names) == 0) {
    stop("No attributes found for clustering", call. = FALSE)
  }

  missing <- setdiff(attr_names, names(data))
  if (length(missing) > 0) {
    stop("Columns not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  attr_df <- sf::st_drop_geometry(data)[, attr_names, drop = FALSE]

  non_numeric <- !sapply(attr_df, is.numeric)
  if (any(non_numeric)) {
    stop("Non-numeric columns: ", paste(names(non_numeric)[non_numeric], collapse = ", "),
         call. = FALSE)
  }

  as.matrix(attr_df)
}

# Validate and prepare spatial weights
# bridge_islands: if TRUE, automatically connect disconnected components
# If FALSE (default), errors on disconnected graphs
prepare_weights <- function(data, weights, bridge_islands = FALSE, call_name = NULL) {
  # Suppress spdep warnings - we'll issue our own
  nb <- suppressWarnings({
    if (is.null(weights) || identical(weights, "queen")) {
      sp_weights(data, type = "queen")
    } else if (identical(weights, "rook")) {
      sp_weights(data, type = "rook")
    } else if (inherits(weights, "nb")) {
      weights
    } else if (is.list(weights) && "type" %in% names(weights)) {
      # Handle list specification like list(type = "knn", k = 6)
      do.call(sp_weights, c(list(data = data), weights))
    } else {
      stop("`weights` must be 'queen', 'rook', an nb object, or a list with type/parameters",
           call. = FALSE)
    }
  })

  # Check connectivity - either error or bridge depending on parameter
  nb <- check_nb_connectivity(nb, bridge_islands = bridge_islands, data = data, call_name = call_name)

  nb
}

# Check neighbor object connectivity
# Returns list with n_components and n_isolates
get_nb_connectivity <- function(nb) {
  n_isolates <- sum(sapply(nb, function(x) length(x) == 1 && x[1] == 0))
  comp <- spdep::n.comp.nb(nb)

  list(
    n_components = comp$nc,
    n_isolates = n_isolates,
    component_labels = comp$comp.id
  )
}

# Bridge disconnected components using KNN
# Adds edges between closest units in different components
bridge_disconnected_nb <- function(nb, data, call_name = NULL) {
  prefix <- if (!is.null(call_name)) paste0(call_name, ": ") else ""

 conn <- get_nb_connectivity(nb)

  if (conn$n_components <= 1 && conn$n_isolates == 0) {
    return(nb)  # Already connected
  }

  # Get centroids for KNN bridging
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(data)))

  # For each disconnected component, find nearest neighbor in different component
  # and add bidirectional edge
  component_labels <- conn$component_labels
  n <- length(nb)

  # Convert nb to mutable list (nb objects can be finicky)
  nb_list <- lapply(nb, function(x) {
    if (length(x) == 1 && x[1] == 0) integer(0) else as.integer(x)
  })

  # Find bridges between components
  # Strategy: for each component except the largest, find the closest point

  # in any other component and add a symmetric edge
  unique_components <- unique(component_labels)

  if (length(unique_components) > 1) {
    bridges_added <- 0

    for (comp in unique_components[-1]) {  # Skip first (usually largest) component
      # Units in this component
      comp_units <- which(component_labels == comp)
      # Units in other components
      other_units <- which(component_labels != comp)

      # Find closest pair between this component and others
      min_dist <- Inf
      best_from <- NA
      best_to <- NA

      for (i in comp_units) {
        for (j in other_units) {
          d <- sqrt(sum((coords[i, ] - coords[j, ])^2))
          if (d < min_dist) {
            min_dist <- d
            best_from <- i
            best_to <- j
          }
        }
      }

      # Add bidirectional edge
      if (!is.na(best_from) && !is.na(best_to)) {
        if (!(best_to %in% nb_list[[best_from]])) {
          nb_list[[best_from]] <- c(nb_list[[best_from]], best_to)
        }
        if (!(best_from %in% nb_list[[best_to]])) {
          nb_list[[best_to]] <- c(nb_list[[best_to]], best_from)
        }
        bridges_added <- bridges_added + 1
      }
    }

    message(
      sprintf("%sBridged %d disconnected component(s) using nearest-neighbor connections.",
              prefix, bridges_added)
    )
  }

  # Handle remaining isolates (units with no neighbors even after bridging)
  for (i in seq_len(n)) {
    if (length(nb_list[[i]]) == 0) {
      # Find nearest neighbor
      dists <- sqrt(rowSums((coords - matrix(coords[i, ], nrow = n, ncol = 2, byrow = TRUE))^2))
      dists[i] <- Inf  # Exclude self
      nearest <- which.min(dists)

      # Add bidirectional edge
      nb_list[[i]] <- c(nb_list[[i]], nearest)
      if (!(i %in% nb_list[[nearest]])) {
        nb_list[[nearest]] <- c(nb_list[[nearest]], i)
      }
    }
  }

  # Convert back to nb object
  # Sort neighbors and ensure proper nb structure
  nb_list <- lapply(nb_list, function(x) {
    if (length(x) == 0) 0L else sort(as.integer(x))
  })

  class(nb_list) <- "nb"
  attr(nb_list, "region.id") <- attr(nb, "region.id")
  attr(nb_list, "call") <- attr(nb, "call")
  attr(nb_list, "type") <- "bridged"
  attr(nb_list, "sym") <- TRUE

  nb_list
}

# Check neighbor object connectivity and error or warn
check_nb_connectivity <- function(nb, bridge_islands = FALSE, data = NULL, call_name = NULL) {
  prefix <- if (!is.null(call_name)) paste0(call_name, ": ") else ""

  conn <- get_nb_connectivity(nb)

  if (conn$n_isolates > 0 || conn$n_components > 1) {
    msg_parts <- character(0)

    if (conn$n_isolates > 0) {
      msg_parts <- c(msg_parts,
        sprintf("%d observation(s) have no neighbors", conn$n_isolates))
    }

    if (conn$n_components > 1) {
      msg_parts <- c(msg_parts,
        sprintf("graph has %d disconnected components", conn$n_components))
    }

    if (bridge_islands && !is.null(data)) {
      # Bridge and return connected nb
      return(bridge_disconnected_nb(nb, data, call_name))
    } else {
      # Error with helpful message
      stop(
        sprintf(
          "%s%s.\nSet `bridge_islands = TRUE` to automatically connect them using nearest-neighbor edges,\nor provide a connected weights object (e.g., KNN weights).",
          prefix, paste(msg_parts, collapse = " and ")
        ),
        call. = FALSE
      )
    }
  }

  nb
}
