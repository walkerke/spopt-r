#' Automatic Zoning Procedure (AZP)
#'
#' Performs regionalization using the Automatic Zoning Procedure algorithm.
#' AZP uses local search to minimize within-region heterogeneity while
#' maintaining spatial contiguity. Three variants are available: basic
#' (greedy), tabu search, and simulated annealing.
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create. May be
#'   omitted when `initial_regions` is supplied, in which case it is taken
#'   from the number of distinct starting regions.
#' @param weights Spatial weights specification. Can be:
#'   \itemize{
#'     \item `"queen"` (default): Polygons sharing any boundary point are neighbors
#'     \item `"rook"`: Polygons sharing an edge are neighbors
#'     \item An `nb` object from spdep or created with [sp_weights()]
#'     \item A list for other weight types: `list(type = "knn", k = 6)` for
#'       k-nearest neighbors, or `list(type = "distance", d = 5000)` for
#'       distance-based weights
#'   }
#' @param bridge_islands Logical. If TRUE, automatically connect disconnected
#'   components (e.g., islands) using nearest-neighbor edges. If FALSE (default),
#'   the function will error when the spatial weights graph is disconnected.
#' @param method Character. Optimization method: "basic" (greedy local search),
#'   "tabu" (tabu search), or "sa" (simulated annealing). Default is "tabu".
#' @param max_iterations Integer. Iteration budget (default 1000). What one
#'   iteration means depends on `method`: for `"basic"` it is a full pass over
#'   all regions (the search stops earlier once a pass makes no move); for
#'   `"tabu"` it is one move (the search stops earlier after
#'   `max(100, 10 * tabu_length)` consecutive non-improving moves); for `"sa"`
#'   it is one sweep over the current border followed by one cooling step.
#' @param tabu_length Integer. Length of tabu list for tabu method (default 10).
#' @param cooling_rate Numeric. Cooling rate for SA method, between 0 and 1,
#'   applied once per sweep (default 0.99).
#' @param initial_temperature Numeric. Initial temperature for SA method.
#'   If 0 (default), it is calibrated from a sample of uphill moves so that a
#'   typical uphill move is accepted with probability 0.8 at the start.
#' @param initial_regions Optional warm start: a vector of length `nrow(data)`
#'   giving a starting region for every area, or an sf object returned by
#'   another spopt regionalization function (its `.region` column is used).
#'   Every starting region must be spatially contiguous under `weights`.
#'   Starting AZP from a SKATER or Ward solution usually beats a random start.
#'   If NULL (default), a random contiguous partition is grown from seeds.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param seed Optional integer for reproducibility.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "azp"
#'     \item method: The optimization method used
#'     \item n_regions: Number of regions created
#'     \item objective: Total within-region sum of squared deviations
#'     \item initial_objective: The same quantity for the starting partition
#'     \item solve_time: Time to solve in seconds
#'   }
#'
#' @details
#' The Automatic Zoning Procedure (AZP) was introduced by Openshaw (1977) and
#' refined by Openshaw & Rao (1995). It is a local search algorithm that:
#'
#' 1. Starts with an initial partition into n_regions (random, or `initial_regions`)
#' 2. Iteratively moves border areas between regions to reduce heterogeneity
#' 3. Maintains spatial contiguity throughout
#' 4. Terminates when no improving moves are found
#'
#' Three variants are available:
#' \itemize{
#'   \item \strong{basic}: Openshaw's greedy local search. Regions are visited
#'     in random order and border areas are pulled in while any improving move
#'     exists. Converges to a local optimum under single-area moves.
#'   \item \strong{tabu}: Steepest-descent / mildest-ascent search with a tabu
#'     list on reversed moves and an aspiration rule, so it can leave local optima.
#'   \item \strong{sa}: Simulated annealing that accepts worse moves with
#'     decreasing probability as temperature cools, one cooling step per sweep
#'     of the border.
#' }
#'
#' The tabu and SA results are finished with a greedy pass, so every returned
#' partition is a local optimum. Move evaluation uses the exact O(p) update of a
#' within-region sum of squares when one area is added or removed, and the
#' contiguity check is only run for moves that would otherwise be accepted.
#'
#' @examples
#' \donttest{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' # Basic AZP with 8 regions
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8)
#'
#' # Tabu search (often finds better solutions)
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
#'               method = "tabu", tabu_length = 15)
#'
#' # Simulated annealing
#' result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 8,
#'               method = "sa", cooling_rate = 0.95)
#'
#' # Warm start from a SKATER solution (usually the best use of AZP)
#' sk <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 8)
#' result <- azp(nc, attrs = c("SID74", "SID79"), initial_regions = sk, method = "sa")
#'
#' # View results
#' plot(result[".region"])
#' }
#'
#' @references
#' Openshaw, S. (1977). A geographical solution to scale and aggregation
#' problems in region-building, partitioning and spatial modelling.
#' Transactions of the Institute of British Geographers, 2(4), 459-472.
#'
#' Openshaw, S., & Rao, L. (1995). Algorithms for reengineering 1991 Census
#' geography. Environment and Planning A, 27(3), 425-446.
#'
#' @export
azp <- function(data,
                attrs = NULL,
                n_regions,
                weights = "queen",
                bridge_islands = FALSE,
                method = c("tabu", "basic", "sa"),
                max_iterations = 1000L,
                tabu_length = 10L,
                cooling_rate = 0.99,
                initial_temperature = 0,
                initial_regions = NULL,
                scale = TRUE,
                seed = NULL,
                verbose = FALSE) {
  # Input validation
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }

  method <- match.arg(method)

  # Warm start: accept a vector or an sf result with a .region column
  init_labels <- NULL
  if (!is.null(initial_regions)) {
    if (inherits(initial_regions, "sf")) {
      if (!".region" %in% names(initial_regions)) {
        stop("`initial_regions` sf object has no `.region` column", call. = FALSE)
      }
      initial_regions <- initial_regions$.region
    }
    if (!is.atomic(initial_regions) || length(initial_regions) != nrow(data)) {
      stop("`initial_regions` must have one entry per row of `data`", call. = FALSE)
    }
    if (anyNA(initial_regions)) {
      stop("`initial_regions` contains missing values", call. = FALSE)
    }
    init_labels <- as.integer(factor(initial_regions)) - 1L
    k_init <- length(unique(init_labels))
    if (missing(n_regions) || is.null(n_regions)) {
      n_regions <- k_init
    } else if (k_init != n_regions) {
      stop(sprintf("`initial_regions` has %d regions but `n_regions` is %d", k_init, n_regions),
           call. = FALSE)
    }
  }

  if (missing(n_regions) || is.null(n_regions) || !is.numeric(n_regions) || n_regions < 2) {
    stop("`n_regions` must be an integer >= 2", call. = FALSE)
  }

  # Determine which columns to check for NAs
  check_cols <- if (!is.null(attrs)) attrs else character(0)

  # Validate data: remove empty geometries, check for NAs
  validated <- validate_regionalization_data(data, check_cols, call_name = "azp")
  data <- validated$data

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  if (!is.numeric(max_iterations) || max_iterations < 1) {
    stop("`max_iterations` must be a positive integer", call. = FALSE)
  }

  if (method == "tabu" && (!is.numeric(tabu_length) || tabu_length < 1)) {
    stop("`tabu_length` must be a positive integer", call. = FALSE)
  }

  if (method == "sa" && (!is.numeric(cooling_rate) || cooling_rate <= 0 || cooling_rate >= 1)) {
    stop("`cooling_rate` must be between 0 and 1 (exclusive)", call. = FALSE)
  }

  # Extract attributes
  attr_matrix <- extract_attrs(data, attrs)

  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  # Prepare spatial weights
  nb <- prepare_weights(data, weights, bridge_islands = bridge_islands, call_name = "azp")

  # Convert nb to adjacency indices
  adj <- nb_to_adj_indices(nb)

  # Every starting region must be contiguous, or AZP cannot repair it
  if (!is.null(init_labels)) {
    for (r in unique(init_labels)) {
      idx <- which(init_labels == r)
      if (length(idx) > 1 &&
          suppressWarnings(spdep::n.comp.nb(spdep::subset.nb(nb, seq_len(n) %in% idx))$nc) > 1) {
        stop(sprintf("Starting region %s in `initial_regions` is not spatially contiguous",
                     unique(initial_regions)[match(r, unique(init_labels))]), call. = FALSE)
      }
    }
  }

  if (verbose) {
    message(sprintf(
      "AZP: n=%d, n_regions=%d, method=%s, attrs=%d%s",
      n, n_regions, method, ncol(attr_matrix),
      if (is.null(init_labels)) "" else ", warm start"
    ))
  }

  # Call Rust implementation
  start_time <- Sys.time()

  result_list <- rust_azp(
    attrs = attr_matrix,
    n_regions = as.integer(n_regions),
    adj_i = adj$i,
    adj_j = adj$j,
    method = method,
    max_iterations = as.integer(max_iterations),
    tabu_length = as.integer(tabu_length),
    cooling_rate = as.numeric(cooling_rate),
    initial_temperature = as.numeric(initial_temperature),
    initial_regions = init_labels,
    seed = if (!is.null(seed)) as.integer(seed) else NULL
  )

  end_time <- Sys.time()

  # Extract results
  labels <- result_list$labels
  objective <- result_list$objective

  # Attach results to sf object
  result <- data
  result$.region <- as.character(labels)

  if (verbose) {
    message(sprintf(
      "  Result: %d regions, objective=%.4f, time=%.3fs",
      length(unique(labels)), objective,
      as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
  }

  # Attach metadata
  metadata <- list(
    algorithm = "azp",
    method = method,
    n_regions = length(unique(labels)),
    objective = objective,
    initial_objective = result_list$initial_objective,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale,
    max_iterations = max_iterations,
    warm_start = !is.null(init_labels)
  )

  attach_spopt_metadata(result, metadata)
}
