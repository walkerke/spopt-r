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
prepare_weights <- function(data, weights) {
  if (is.null(weights) || identical(weights, "queen")) {
    sp_weights(data, type = "queen")
  } else if (identical(weights, "rook")) {
    sp_weights(data, type = "rook")
  } else if (inherits(weights, "nb")) {
    weights
  } else {
    stop("`weights` must be 'queen', 'rook', or an nb object", call. = FALSE)
  }
}
