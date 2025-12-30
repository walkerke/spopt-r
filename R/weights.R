#' Create spatial weights from an sf object
#'
#' Constructs spatial weights (neighborhood structure) from sf geometries.
#' Wraps spdep functions with a convenient interface.
#'
#' @param data An sf object with polygon or point geometries.
#' @param type Type of weights. One of:
#'   \itemize{
#'     \item "queen" (default): Polygons sharing any boundary point are neighbors
#'     \item "rook": Polygons sharing an edge are neighbors
#'     \item "knn": K-nearest neighbors based on centroid distance
#'     \item "distance": All units within a distance threshold are neighbors
#'   }
#' @param k Number of nearest neighbors. Required when `type = "knn"`.
#' @param d Distance threshold. Required when `type = "distance"`. Units match
#'   the CRS of the data (e.g., meters for projected CRS).
#' @param ... Additional arguments passed to spdep functions.
#'
#' @return A neighbors list object (class "nb") compatible with spdep.
#'
#' @details
#' **Choosing a weight type:**
#'
#' - Use **queen/rook** for polygon data where physical adjacency matters
#' - Use **knn** when you need guaranteed connectivity (no isolates) or for point data
#' - Use **distance** for point data or when interaction depends on proximity
#'
#' **KNN weights** always produce a connected graph (if k >= 1), making them
#' useful for datasets with islands or disconnected polygons.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Queen contiguity (default)
#' w_queen <- sp_weights(nc, type = "queen")
#'
#' # K-nearest neighbors (guarantees connectivity)
#' w_knn <- sp_weights(nc, type = "knn", k = 6)
#'
#' # Distance-based (e.g., 50km for projected data)
#' nc_proj <- st_transform(nc, 32119)  # NC State Plane
#' w_dist <- sp_weights(nc_proj, type = "distance", d = 50000)
#' }
#'
#' @export
sp_weights <- function(data, type = c("queen", "rook", "knn", "distance"),
                       k = NULL, d = NULL, ...) {
  type <- match.arg(type)

  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  if (type == "knn") {
    if (is.null(k)) {
      stop("`k` must be specified when type = 'knn'", call. = FALSE)
    }
    if (!is.numeric(k) || k < 1) {
      stop("`k` must be a positive integer", call. = FALSE)
    }
    # Get centroids for knn
    coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(data)))
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = as.integer(k)), ...)

  } else if (type == "distance") {
    if (is.null(d)) {
      stop("`d` must be specified when type = 'distance'", call. = FALSE)
    }
    if (!is.numeric(d) || d <= 0) {
      stop("`d` must be a positive number", call. = FALSE)
    }
    # Get centroids for distance-based
    coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(data)))
    nb <- spdep::dnearneigh(coords, d1 = 0, d2 = d, ...)

  } else {
    # Queen or rook contiguity
    queen <- type == "queen"
    nb <- spdep::poly2nb(data, queen = queen, ...)
  }

  nb
}

# Convert neighbor list to sparse adjacency matrix (internal)
nb_to_sparse <- function(nb, style = "B") {
  n <- length(nb)
  lw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)

  # Build sparse matrix from listw
  i <- rep(seq_len(n), times = spdep::card(nb))
  j <- unlist(nb)
  x <- unlist(lw$weights)

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
}
