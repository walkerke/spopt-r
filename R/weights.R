#' Create spatial weights from an sf object
#'
#' Constructs spatial weights (neighborhood structure) from polygon geometries.
#' Wraps spdep functions with a convenient interface.
#'
#' @param data An sf object with polygon geometries.
#' @param type Type of contiguity. One of "queen" (default), "rook", or "knn".
#' @param k Number of nearest neighbors (only used when type = "knn").
#' @param ... Additional arguments passed to spdep functions.
#'
#' @return A neighbors list object (class "nb") compatible with spdep.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#' w <- sp_weights(nc, type = "queen")
#' }
#'
#' @export
sp_weights <- function(data, type = c("queen", "rook", "knn"), k = NULL, ...) {
  type <- match.arg(type)


  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  if (type == "knn") {
    if (is.null(k)) {
      stop("`k` must be specified when type = 'knn'", call. = FALSE)
    }
    # Get centroids for knn
    coords <- sf::st_coordinates(sf::st_centroid(data))
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k), ...)
  } else {
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
