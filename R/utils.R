#' Attach spopt metadata to result
#'
#' @param result An sf object to attach metadata to.
#' @param metadata A list of algorithm metadata.
#'
#' @return The sf object with "spopt" attribute attached.
#'
#' @keywords internal
attach_spopt_metadata <- function(result, metadata) {
  attr(result, "spopt") <- metadata
  result
}

#' Extract attribute columns from sf
#'
#' @param data An sf object.
#' @param attrs A character vector of column names (e.g., `c("var1", "var2")`).
#'   If NULL, uses all numeric columns.
#'
#' @return A numeric matrix of attribute values.
#'
#' @keywords internal
extract_attrs <- function(data, attrs) {
  if (is.null(attrs)) {
    # Use all numeric columns except geometry
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

  # Check all columns exist
  missing <- setdiff(attr_names, names(data))
  if (length(missing) > 0) {
    stop("Columns not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  # Extract and convert to matrix
  attr_df <- sf::st_drop_geometry(data)[, attr_names, drop = FALSE]

  # Check all numeric
  non_numeric <- !sapply(attr_df, is.numeric)
  if (any(non_numeric)) {
    stop("Non-numeric columns: ", paste(names(non_numeric)[non_numeric], collapse = ", "),
         call. = FALSE)
  }

  as.matrix(attr_df)
}

#' Validate and prepare spatial weights
#'
#' @param data An sf object.
#' @param weights Either an nb object, or one of "queen", "rook".
#'
#' @return An nb object.
#'
#' @keywords internal
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
