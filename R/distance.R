#' Compute distance matrix between sf objects
#'
#' Computes pairwise distances between geometries. For geographic (longlat)
#' coordinate systems, uses great circle distances in meters via [sf::st_distance()].
#' For projected coordinate systems, uses fast Euclidean distance in the CRS units
#' (typically meters).
#'
#' @param x An sf object (demand points for facility location, or areas for regionalization).
#' @param y An sf object (facility locations). If NULL, computes distances within x.
#' @param type Distance type: "euclidean" (default) or "manhattan". Note that for
#'   geographic CRS, only "euclidean" (great circle) distance is available.
#' @param use_centroids Logical. If TRUE (default for polygons), use polygon centroids.
#'
#' @return A numeric matrix of distances. Rows correspond to x, columns to y.
#'   For geographic CRS, distances are in meters. For projected CRS, distances
#'   are in the CRS units (usually meters).
#'
#' @examples
#' \dontrun{
#' library(sf)
#' demand <- st_as_sf(data.frame(x = runif(10), y = runif(10)), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(5), y = runif(5)), coords = c("x", "y"))
#' d <- distance_matrix(demand, facilities)
#' }
#'
#' @export
distance_matrix <- function(x, y = NULL, type = c("euclidean", "manhattan"),
                            use_centroids = NULL) {
  type <- match.arg(type)

  if (!inherits(x, "sf")) {
    stop("`x` must be an sf object", call. = FALSE)
  }
  if (!is.null(y) && !inherits(y, "sf")) {
    stop("`y` must be an sf object or NULL", call. = FALSE)
  }

  # Check if data is in geographic (longlat) CRS
  is_longlat <- sf::st_is_longlat(x)

  # For geographic CRS, use sf::st_distance() which handles great circle distances
  if (isTRUE(is_longlat)) {
    if (type == "manhattan") {
      warning("Manhattan distance not available for geographic CRS; using great circle distance.")
    }
    return(distance_matrix_geographic(x, y, use_centroids))
  }

  # For projected CRS, use fast Euclidean/Manhattan calculation
  # Determine if we should use centroids
  geom_type_x <- unique(sf::st_geometry_type(x))
  is_polygon_x <- any(geom_type_x %in% c("POLYGON", "MULTIPOLYGON"))

  if (is.null(use_centroids)) {
    use_centroids <- is_polygon_x
  }

  # Extract coordinates
  if (use_centroids || is_polygon_x) {
    coords_x <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(x)))
  } else {
    coords_x <- sf::st_coordinates(x)
  }

  if (is.null(y)) {
    coords_y <- coords_x
  } else {
    geom_type_y <- unique(sf::st_geometry_type(y))
    is_polygon_y <- any(geom_type_y %in% c("POLYGON", "MULTIPOLYGON"))

    if (use_centroids || is_polygon_y) {
      coords_y <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(y)))
    } else {
      coords_y <- sf::st_coordinates(y)
    }
  }

  # Compute distances - call Rust for large matrices, pure R for small
  n_x <- nrow(coords_x)
  n_y <- nrow(coords_y)

  # Use Rust for large matrices (>1000 elements), pure R for small
  use_rust <- (n_x * n_y) > 1000

  if (type == "euclidean") {
    if (use_rust) {
      dist_matrix <- rust_distance_matrix_euclidean(
        coords_x[, 1], coords_x[, 2],
        coords_y[, 1], coords_y[, 2]
      )
    } else {
      dist_matrix <- sqrt(outer(coords_x[, 1], coords_y[, 1], "-")^2 +
                          outer(coords_x[, 2], coords_y[, 2], "-")^2)
    }
  } else if (type == "manhattan") {
    if (use_rust) {
      dist_matrix <- rust_distance_matrix_manhattan(
        coords_x[, 1], coords_x[, 2],
        coords_y[, 1], coords_y[, 2]
      )
    } else {
      dist_matrix <- abs(outer(coords_x[, 1], coords_y[, 1], "-")) +
                     abs(outer(coords_x[, 2], coords_y[, 2], "-"))
    }
  }

  dist_matrix
}

#' Compute distance matrix for geographic CRS using sf::st_distance
#'
#' Internal function that uses sf::st_distance() for accurate great circle
#' distances on geographic (longlat) coordinate systems.
#'
#' @param x An sf object
#' @param y An sf object or NULL
#' @param use_centroids Logical. If TRUE, use polygon centroids.
#'
#' @return A numeric matrix of distances in meters.
#'
#' @keywords internal
distance_matrix_geographic <- function(x, y = NULL, use_centroids = NULL) {
  # Determine if we should use centroids
  geom_type_x <- unique(sf::st_geometry_type(x))
  is_polygon_x <- any(geom_type_x %in% c("POLYGON", "MULTIPOLYGON"))

  if (is.null(use_centroids)) {
    use_centroids <- is_polygon_x
  }

  # Get geometries (with centroids if needed)
  if (use_centroids || is_polygon_x) {
    geom_x <- sf::st_centroid(sf::st_geometry(x))
  } else {
    geom_x <- sf::st_geometry(x)
  }

  if (is.null(y)) {
    geom_y <- geom_x
  } else {
    geom_type_y <- unique(sf::st_geometry_type(y))
    is_polygon_y <- any(geom_type_y %in% c("POLYGON", "MULTIPOLYGON"))

    if (use_centroids || is_polygon_y) {
      geom_y <- sf::st_centroid(sf::st_geometry(y))
    } else {
      geom_y <- sf::st_geometry(y)
    }
  }

  # Use sf::st_distance which handles geographic CRS properly (returns meters)
  dist_units <- sf::st_distance(geom_x, geom_y)

  # Convert to plain numeric matrix (dropping units)
  dist_matrix <- as.matrix(dist_units)
  class(dist_matrix) <- "matrix"
  attr(dist_matrix, "units") <- NULL

  dist_matrix
}
