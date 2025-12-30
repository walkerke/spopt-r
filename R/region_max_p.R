#' Max-P Regions
#'
#' Perform Max-P regionalization to maximize the number of spatially contiguous
#' regions such that each region satisfies a minimum threshold constraint on a
#' specified attribute. This is useful for creating regions that meet minimum
#' population or sample size requirements.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for computing
#'   within-region dissimilarity (e.g., `c("var1", "var2")`). If NULL,
#'   uses all numeric columns.
#' @param threshold_var Character. Name of the column containing the threshold
#'   variable (e.g., population, income). Each region must have a sum of this
#'   variable >= `threshold`.
#' @param threshold Numeric. Minimum sum of `threshold_var` required per region.
#' @param weights Spatial weights specification. Can be:
#'   \itemize{
#'     \item `"queen"` (default): Polygons sharing any boundary point are neighbors
#'     \item `"rook"`: Polygons sharing an edge are neighbors
#'     \item An `nb` object from spdep or created with [sp_weights()]
#'     \item A list for other weight types: `list(type = "knn", k = 6)` for
#'       k-nearest neighbors, or `list(type = "distance", d = 5000)` for
#'       distance-based weights
#'   }
#'   KNN weights guarantee connectivity (no islands), which can be useful
#'   for datasets with disconnected polygons.
#' @param bridge_islands Logical. If TRUE, automatically connect disconnected
#'   components (e.g., islands) using nearest-neighbor edges. If FALSE (default),
#'   the function will error when the spatial weights graph is disconnected.
#'   This is useful for datasets like LA County with Catalina Islands, or
#'   archipelago data where physical adjacency doesn't exist but regionalization
#'   is still desired.
#' @param compact Logical. If TRUE, optimize for region compactness in addition
#'   to attribute homogeneity. Compact regions have more regular shapes, which
#'   is useful for sales territories, patrol areas, and electoral districts.
#'   Default is FALSE.
#' @param compact_weight Numeric between 0 and 1. Weight for compactness vs
#'   attribute homogeneity when `compact = TRUE`. Higher values prioritize
#'   compact shapes over attribute similarity. Default is 0.5.
#' @param n_iterations Integer. Number of construction phase iterations (default 100).
#'   Higher values explore more random starting solutions.
#' @param n_sa_iterations Integer. Number of simulated annealing iterations (default 100).
#'   Set to 0 to skip the SA refinement phase.
#' @param cooling_rate Numeric. SA cooling rate between 0 and 1 (default 0.99).
#'   Smaller values cool faster, larger values allow more exploration.
#' @param tabu_length Integer. Length of tabu list for SA phase (default 10).
#' @param scale Logical. If TRUE (default), standardize attributes before
#'   computing dissimilarity.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing region assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "max_p"
#'     \item n_regions: Number of regions created (the "p" in max-p)
#'     \item objective: Total within-region sum of squared deviations
#'     \item threshold_var: Name of threshold variable
#'     \item threshold: Threshold value used
#'     \item solve_time: Time to solve in seconds
#'     \item mean_compactness: Mean Polsby-Popper compactness (if `compact = TRUE`)
#'     \item region_compactness: Per-region compactness scores (if `compact = TRUE`)
#'   }
#'
#' @details
#' The Max-P algorithm (Duque, Anselin & Rey, 2012; Wei, Rey & Knaap, 2021)
#' solves the problem of aggregating n geographic areas into the maximum number
#' of homogeneous regions while ensuring:
#'
#' 1. Each region is spatially contiguous (connected)
#' 2. Each region satisfies a minimum threshold on a specified attribute
#'
#' The algorithm has two phases:
#' \enumerate{
#'   \item Construction phase: Builds feasible solutions via randomized greedy
#'     region growing. Multiple random starts are explored in parallel.
#'   \item Simulated annealing phase: Refines solutions by moving border areas
#'     between regions to minimize within-region dissimilarity while respecting
#'     constraints.
#' }
#'
#' When `compact = TRUE`, the algorithm additionally optimizes for compact region
#' shapes based on Feng, Rey, & Wei (2022). Compact regions:
#' \itemize{
#'   \item Minimize travel time within regions (useful for service territories)
#'   \item Reduce gerrymandering potential (electoral districts)
#'   \item Often result in finding MORE regions due to efficient space usage
#' }
#'
#' \strong{Compactness metric}: This implementation uses a centroid dispersion
#' measure during optimization, rather than the Normalized Moment of Inertia (NMI)
#' described in Feng et al. (2022). This design choice provides two advantages:
#' \enumerate{
#'   \item \strong{Point-based regionalization}: The algorithm works with both
#'     polygon and point geometries. For point data, use KNN or distance-based
#'     weights (e.g., `weights = list(type = "knn", k = 6)`).
#'   \item \strong{Computational efficiency}: Centroid dispersion is O(n) per
#'     region versus O(v) for NMI where v = total polygon vertices.
#' }
#' For polygon data, centroids are computed via [sf::st_centroid()]. Users should
#' be aware that centroid-based compactness may be less accurate for highly
#' irregular shapes or large, sparsely-populated areas where the centroid poorly
#' represents the polygon's spatial extent.
#'
#' The reported `mean_compactness` and `region_compactness` in results use
#' Polsby-Popper (4*pi*A/P^2), a standard geometric compactness measure for
#' polygons. For point data, these metrics are not computed.
#'
#' This implementation is optimized for speed using:
#' \itemize{
#'   \item Parallel construction with early termination
#'   \item Efficient articulation point detection for move eligibility
#'   \item Incremental threshold tracking
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Create regions where each has at least 100,000 in BIR74
#' result <- max_p_regions(
#'   nc,
#'   attrs = c("SID74", "SID79"),
#'   threshold_var = "BIR74",
#'   threshold = 100000
#' )
#'
#' # Check number of regions created
#' attr(result, "spopt")$n_regions
#'
#' # With compactness optimization (for sales territories)
#' result_compact <- max_p_regions(
#'   nc,
#'   attrs = c("SID74", "SID79"),
#'   threshold_var = "BIR74",
#'   threshold = 100000,
#'   compact = TRUE,
#'   compact_weight = 0.5
#' )
#'
#' # Check compactness
#' attr(result_compact, "spopt")$mean_compactness
#'
#' # Plot results
#' plot(result[".region"])
#'
#' # Point-based regionalization (e.g., store locations, sensor networks)
#' # Use KNN weights since points don't have polygon contiguity
#' points <- st_as_sf(data.frame(
#'   x = runif(200), y = runif(200),
#'   customers = rpois(200, 100),
#'   avg_income = rnorm(200, 50000, 15000)
#' ), coords = c("x", "y"))
#'
#' result_points <- max_p_regions(
#'   points,
#'   attrs = "avg_income",
#'   threshold_var = "customers",
#'   threshold = 500,
#'   weights = list(type = "knn", k = 6),
#'   compact = TRUE
#' )
#' }
#'
#' @references
#' Duque, J. C., Anselin, L., & Rey, S. J. (2012). The max-p-regions problem.
#' Journal of Regional Science, 52(3), 397-419.
#'
#' Wei, R., Rey, S., & Knaap, E. (2021). Efficient regionalization for spatially
#' explicit neighborhood delineation. International Journal of Geographical
#' Information Science, 35(1), 135-151. \doi{10.1080/13658816.2020.1759806}
#'
#' Feng, X., Rey, S., & Wei, R. (2022). The max-p-compact-regions problem.
#' Transactions in GIS, 26, 717-734. \doi{10.1111/tgis.12874}
#'
#' @export
max_p_regions <- function(data,
                          attrs = NULL,
                          threshold_var,
                          threshold,
                          weights = "queen",
                          bridge_islands = FALSE,
                          compact = FALSE,
                          compact_weight = 0.5,
                          n_iterations = 100L,
                          n_sa_iterations = 100L,
                          cooling_rate = 0.99,
                          tabu_length = 10L,
                          scale = TRUE,
                          seed = NULL,
                          verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  # Validate threshold variable exists before geometry check
  if (missing(threshold_var) || !is.character(threshold_var) || length(threshold_var) != 1) {
    stop("`threshold_var` must be a single column name", call. = FALSE)
  }
  if (!threshold_var %in% names(data)) {
    stop(paste0("Threshold variable '", threshold_var, "' not found in data"), call. = FALSE)
  }

  # Determine which columns to check for NAs
  check_cols <- c(threshold_var, if (!is.null(attrs)) attrs else character(0))

  # Validate data: remove empty geometries, check for NAs

  validated <- validate_regionalization_data(data, check_cols, call_name = "max_p_regions")
  data <- validated$data

  n <- nrow(data)
  if (n < 2) {
    stop("Data must have at least 2 observations", call. = FALSE)
  }

  # Extract threshold values (already validated for NAs)
  threshold_values <- as.numeric(data[[threshold_var]])

  # Validate threshold
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop("`threshold` must be a positive number", call. = FALSE)
  }

  total_threshold <- sum(threshold_values)
  if (total_threshold < threshold) {
    stop(paste0(
      "Total of threshold variable (", round(total_threshold, 2),
      ") is less than threshold (", threshold,
      "). No valid regionalization possible."
    ), call. = FALSE)
  }

  # Validate other parameters
  if (!is.numeric(n_iterations) || n_iterations < 1) {
    stop("`n_iterations` must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_sa_iterations) || n_sa_iterations < 0) {
    stop("`n_sa_iterations` must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1) {
    stop("`cooling_rate` must be between 0 and 1 (exclusive)", call. = FALSE)
  }
  if (!is.numeric(tabu_length) || tabu_length < 0) {
    stop("`tabu_length` must be a non-negative integer", call. = FALSE)
  }
  if (!is.logical(compact) || length(compact) != 1) {
    stop("`compact` must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(compact_weight) || compact_weight < 0 || compact_weight > 1) {
    stop("`compact_weight` must be between 0 and 1", call. = FALSE)
  }
  if (!is.logical(bridge_islands) || length(bridge_islands) != 1) {
    stop("`bridge_islands` must be TRUE or FALSE", call. = FALSE)
  }

  # Extract centroids if compactness is enabled
  centroids_x <- NULL
  centroids_y <- NULL
  if (compact) {
    # Get centroids of all units
    cents <- sf::st_centroid(sf::st_geometry(data))
    coords <- sf::st_coordinates(cents)
    centroids_x <- as.numeric(coords[, 1])
    centroids_y <- as.numeric(coords[, 2])
  }

  # Extract attributes for dissimilarity
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights, bridge_islands = bridge_islands, call_name = "max_p_regions")

  # Convert nb to sparse matrix indices
  adj <- nb_to_adj_indices(nb)

  if (verbose) {
    message(sprintf(
      "Max-P: n=%d, threshold=%.2f (var=%s), attrs=%d, edges=%d",
      n, threshold, threshold_var, ncol(attr_matrix), length(adj$i)
    ))
    message(sprintf(
      "  Construction iterations: %d, SA iterations: %d",
      n_iterations, n_sa_iterations
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_max_p(
    attrs = attr_matrix,
    threshold_var = threshold_values,
    threshold = as.numeric(threshold),
    adj_i = adj$i,
    adj_j = adj$j,
    n_iterations = as.integer(n_iterations),
    n_sa_iterations = as.integer(n_sa_iterations),
    cooling_rate = as.numeric(cooling_rate),
    tabu_length = as.integer(tabu_length),
    seed = if (!is.null(seed)) as.integer(seed) else NULL,
    centroids_x = centroids_x,
    centroids_y = centroids_y,
    compact = compact,
    compact_weight = as.numeric(compact_weight)
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  n_regions <- result_list$n_regions
  objective <- result_list$objective

  # Attach results to sf object
  result <- data
  result$.region <- as.integer(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      n_regions, objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Compute region statistics
  region_stats <- compute_region_stats(result, threshold_var, threshold)

  # Compute Polsby-Popper compactness if compact mode was used
  compactness_metrics <- NULL
  if (compact) {
    compactness_metrics <- compute_polsby_popper(result)
  }

  # Attach metadata
  metadata <- list(
    algorithm = "max_p",
    n_regions = n_regions,
    objective = objective,
    threshold_var = threshold_var,
    threshold = threshold,
    region_stats = region_stats,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale,
    n_iterations = n_iterations,
    n_sa_iterations = n_sa_iterations,
    compact = compact,
    compact_weight = compact_weight,
    mean_compactness = if (!is.null(compactness_metrics)) compactness_metrics$mean else NULL,
    region_compactness = if (!is.null(compactness_metrics)) compactness_metrics$by_region else NULL
  )

  attach_spopt_metadata(result, metadata)
}

# Compute region statistics (internal)
compute_region_stats <- function(result, threshold_var, threshold) {
  # Compute stats for each region
  regions <- unique(result$.region)
  stats <- lapply(regions, function(r) {
    idx <- result$.region == r
    list(
      region = r,
      n_areas = sum(idx),
      threshold_sum = sum(result[[threshold_var]][idx]),
      meets_threshold = sum(result[[threshold_var]][idx]) >= threshold
    )
  })

  do.call(rbind, lapply(stats, as.data.frame))
}

# Compute Polsby-Popper compactness for regions (internal)
# Polsby-Popper = 4 * pi * Area / Perimeter^2 (0-1, 1 = circle)
# Returns NULL for point geometries (compactness not meaningful)
compute_polsby_popper <- function(result) {
  # Check geometry type - Polsby-Popper only makes sense for polygons

  geom_type <- sf::st_geometry_type(result, by_geometry = FALSE)
  if (geom_type %in% c("POINT", "MULTIPOINT")) {
    return(NULL)
  }

  regions <- sort(unique(result$.region))
  by_region <- numeric(length(regions))

  for (i in seq_along(regions)) {
    r <- regions[i]
    region_geom <- sf::st_union(sf::st_geometry(result[result$.region == r, ]))

    # Compute area and perimeter
    area <- as.numeric(sf::st_area(region_geom))

    # Compute perimeter - use boundary length
    boundary <- sf::st_boundary(region_geom)
    if (inherits(boundary, "sfc_MULTILINESTRING")) {
      boundary <- sf::st_cast(boundary, "LINESTRING")
    }
    perimeter <- sum(as.numeric(sf::st_length(boundary)))

    # Polsby-Popper: 4 * pi * A / P^2
    if (perimeter > 0) {
      by_region[i] <- 4 * pi * area / (perimeter^2)
    } else {
      by_region[i] <- 0
    }
  }

  names(by_region) <- as.character(regions)

  list(
    mean = mean(by_region, na.rm = TRUE),
    by_region = by_region
  )
}
