#' Spatially-Encouraged Spectral Clustering (SPENC)
#'
#' Performs spectral clustering with spatial constraints by combining
#' spatial connectivity with attribute similarity using kernel methods.
#' This approach is useful for clustering with highly non-convex clusters
#' or irregular topologies in geographic contexts.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create.
#' @param weights Spatial weights specification. One of "queen" (default),
#'   "rook", or an nb object from spdep.
#' @param gamma Numeric. RBF kernel parameter controlling attribute similarity
#'   decay. Larger values = faster decay = more local similarity. Default is 1.
#'   Can also be "auto" to estimate from data.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "spenc"
#'     \item n_regions: Number of regions created
#'     \item objective: Within-cluster sum of squared distances in embedding space
#'     \item gamma: The gamma parameter used
#'     \item solve_time: Time to solve in seconds
#'   }
#'
#' @details
#' SPENC (Wolf, 2020) extends spectral clustering to incorporate spatial
#' constraints. The algorithm:
#'
#' 1. Computes attribute affinity using an RBF (Gaussian) kernel
#' 2. Multiplies element-wise with spatial weights (only neighbors have affinity)
#' 3. Computes the normalized Laplacian of the combined affinity matrix
#' 4. Extracts the k smallest eigenvectors as a spectral embedding

#' 5. Applies k-means clustering to the embedding
#'
#' Key advantages:
#' \itemize{
#'   \item Can find non-convex cluster shapes
#'   \item Respects spatial connectivity
#'   \item Balances attribute similarity with spatial proximity
#' }
#'
#' The gamma parameter controls how quickly attribute similarity decays with
#' distance in attribute space. Larger values create more localized clusters.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Basic SPENC with 8 regions
#' result <- spenc(nc, attrs = c("SID74", "SID79"), n_regions = 8)
#'
#' # Adjust gamma for different cluster tightness
#' result <- spenc(nc, attrs = c("SID74", "SID79"), n_regions = 8, gamma = 0.5)
#'
#' # View results
#' plot(result[".region"])
#' }
#'
#' @references
#' Wolf, L.J. (2020). Spatially-encouraged spectral clustering: a technique
#' for blending map typologies and regionalization. OSF Preprints.
#'
#' @export
spenc <- function(data,
                  attrs = NULL,
                  n_regions,
                  weights = "queen",
                  gamma = 1.0,
                  scale = TRUE,
                  seed = NULL,
                  verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  if (!is.numeric(n_regions) || n_regions < 2) {
    stop("`n_regions` must be an integer >= 2", call. = FALSE)
  }

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  # Extract attributes
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Handle gamma = "auto"
  if (is.character(gamma) && gamma == "auto") {
    # Estimate gamma as 1 / (2 * median squared distance)
    if (nrow(attr_matrix) > 100) {
      # Sample for efficiency
      idx <- sample(nrow(attr_matrix), 100)
      sample_attrs <- attr_matrix[idx, , drop = FALSE]
    } else {
      sample_attrs <- attr_matrix
    }
    dists <- as.vector(stats::dist(sample_attrs))
    gamma <- 1 / (2 * stats::median(dists)^2)
    if (!is.finite(gamma) || gamma <= 0) {
      gamma <- 1.0
    }
  }

  if (!is.numeric(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive number or 'auto'", call. = FALSE)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights)

  # Convert nb to adjacency indices
  adj <- nb_to_adj_indices(nb)

  if (verbose) {
    message(sprintf(
      "SPENC: n=%d, n_regions=%d, gamma=%.4f, attrs=%d",
      n, n_regions, gamma, ncol(attr_matrix)
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_spenc(
    attrs = attr_matrix,
    n_regions = as.integer(n_regions),
    adj_i = adj$i,
    adj_j = adj$j,
    gamma = as.numeric(gamma),
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  objective <- result_list$objective

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      length(unique(labels)), objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Attach metadata
  metadata <- list(
    algorithm = "spenc",
    n_regions = length(unique(labels)),
    objective = objective,
    gamma = gamma,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale
  )

  attach_spopt_metadata(result, metadata)
}
