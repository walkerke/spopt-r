#' Compute distance matrix between sf objects
#'
#' Computes pairwise distances between geometries. For points, uses
#' Euclidean distance. For polygons, uses centroid-to-centroid distance
#' by default.
#'
#' @param x An sf object (demand points for facility location, or areas for regionalization).
#' @param y An sf object (facility locations). If NULL, computes distances within x.
#' @param type Distance type: "euclidean" (default) or "manhattan".
#' @param use_centroids Logical. If TRUE (default for polygons), use polygon centroids.
#'
#' @return A numeric matrix of distances. Rows correspond to x, columns to y.
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
