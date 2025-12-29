#' SKATER Spatial Clustering
#'
#' Performs spatial clustering using the SKATER algorithm (Spatial 'K'luster
#' Analysis by Tree Edge Removal). The algorithm builds a minimum spanning
#' tree from the spatial contiguity graph, then iteratively removes edges
#' to create spatially contiguous clusters.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create.
#' @param weights Spatial weights specification. One of "queen" (default),
#'   "rook", or an nb object from spdep.
#' @param floor Optional. Column name specifying a floor constraint variable.
#' @param floor_value Numeric. Minimum sum of floor variable required per region.
#'   Only used if `floor` is specified.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Cluster into 5 regions based on SIDS rates
#' result <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 5)
#'
#' # With floor constraint: each region must have at least 100,000 births
#' result <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 5,
#'                  floor = "BIR74", floor_value = 100000)
#'
#' # View results
#' plot(result[".region"])
#' }
#'
#' @references
#' Assuncao, R. M., Neves, M. C., Camara, G., & Freitas, C. da C. (2006).
#' Efficient regionalization techniques for socio-economic geographical units
#' using minimum spanning trees. International Journal of Geographical
#' Information Science, 20(7), 797-811.
#'
#' @export
skater <- function(data,
                   attrs = NULL,
                   n_regions,
                   weights = "queen",
                   floor = NULL,
                   floor_value = 0,
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

  # Determine which columns to check for NAs
  check_cols <- c(
    if (!is.null(attrs)) attrs else character(0),
    if (!is.null(floor)) floor else character(0)
  )

  # Validate data: remove empty geometries, check for NAs
  validated <- validate_regionalization_data(data, check_cols, call_name = "skater")
  data <- validated$data

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  # Extract attributes
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights)

  # Convert nb to sparse matrix indices
  adj <- nb_to_adj_indices(nb)

  # Handle floor constraint
  floor_var <- NULL
  if (!is.null(floor)) {
    if (!floor %in% names(data)) {
      stop(paste0("Floor variable '", floor, "' not found in data"), call. = FALSE)
    }
    floor_var <- as.numeric(data[[floor]])
  }

  # Call Rust implementation
  start_time <- Sys.time()

  labels <- rust_skater(
    attrs = attr_matrix,
    adj_i = adj$i,
    adj_j = adj$j,
    n_regions = as.integer(n_regions),
    floor_var = floor_var,
    floor_value = as.numeric(floor_value),
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

end_time <- Sys.time()

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  # Compute objective (total within-cluster SSD)
  objective <- compute_ssd(attr_matrix, labels)

  # Attach metadata
  metadata <- list(
    algorithm = "skater",
    n_regions = length(unique(labels)),
    objective = objective,
    floor = floor,
    floor_value = floor_value,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale
  )

  attach_spopt_metadata(result, metadata)
}

# Convert nb object to adjacency indices (internal)
nb_to_adj_indices <- function(nb) {
  n <- length(nb)
  i <- integer(0)
  j <- integer(0)

  for (idx in seq_len(n)) {
    neighbors <- nb[[idx]]
    if (length(neighbors) > 0 && neighbors[1] != 0) {
      i <- c(i, rep(idx - 1L, length(neighbors)))  # 0-based for Rust
      j <- c(j, neighbors - 1L)
    }
  }

  list(i = i, j = j)
}

# Compute total sum of squared deviations (internal)
compute_ssd <- function(attrs, labels) {
  unique_labels <- unique(labels)
  total_ssd <- 0

  for (lab in unique_labels) {
    cluster_idx <- which(labels == lab)
    if (length(cluster_idx) > 1) {
      cluster_attrs <- attrs[cluster_idx, , drop = FALSE]
      centroid <- colMeans(cluster_attrs)
      ssd <- sum(sweep(cluster_attrs, 2, centroid)^2)
      total_ssd <- total_ssd + ssd
    }
  }

  total_ssd
}
