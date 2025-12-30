#' Huff Model for Market Share Analysis
#'
#' Computes probability surfaces to predict market share and sales potential
#' based on distance decay and store attractiveness. The Huff model is widely
#' used in retail site selection to estimate the probability that a consumer
#' at a given location will choose a particular store.
#'
#' @param demand An sf object representing demand points or areas. Can be
#'   customer locations, census block groups, grid cells, etc.
#' @param stores An sf object representing store/facility locations.
#' @param attractiveness_col Character vector. Column name(s) in `stores`
#'   containing attractiveness values (e.g., square footage, parking spaces).
#'   Multiple columns can be specified for composite attractiveness.
#' @param attractiveness_exponent Numeric vector. Exponent(s) for attractiveness
#'   (default 1). Must be same length as `attractiveness_col` or length 1
#'   (recycled). Higher values increase the importance of that variable.
#' @param distance_exponent Numeric. Distance decay exponent (default -1.5).
#'
#'   Should be negative; more negative = faster decay with distance.
#' @param sales_potential_col Optional character. Column name in `demand`
#'   containing sales potential values (e.g., disposable income, population).
#'   If NULL, each demand point is weighted equally.
#' @param cost_matrix Optional. Pre-computed distance/cost matrix (demand x stores).
#'   If NULL, Euclidean distance is computed from geometries.
#' @param distance_metric Distance metric if cost_matrix is NULL:
#'   "euclidean" (default) or "manhattan".
#'
#' @return A list with:
#'   \itemize{
#'     \item `$demand`: Original demand sf with added columns:
#'       \itemize{
#'         \item `.primary_store`: ID of highest-probability store
#'         \item `.entropy`: Competition measure (higher = more competition)
#'         \item `.prob_<store_id>`: Probability columns for each store
#'       }
#'     \item `$stores`: Original stores sf with added columns:
#'       \itemize{
#'         \item `.market_share`: Proportion of total market captured
#'         \item `.expected_sales`: Expected sales (sum of prob × potential)
#'       }
#'     \item `$probability_matrix`: Full probability matrix (n_demand × n_stores)
#'   }
#'   Metadata in "spopt" attribute includes parameters used.
#'
#' @details
#' The Huff model calculates the probability that a consumer at location i
#' will choose store j using:
#'
#' \deqn{P_{ij} = \frac{A_j \times D_{ij}^\beta}{\sum_k A_k \times D_{ik}^\beta}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{A_j} is the composite attractiveness of store j
#'   \item \eqn{D_{ij}} is the distance from i to j
#'   \item \eqn{\beta} is the distance decay exponent (default -1.5)
#' }
#'
#' When multiple attractiveness variables are specified, the composite
#' attractiveness is computed as:
#'
#' \deqn{A_j = \prod_m V_{jm}^{\alpha_m}}
#'
#' Where \eqn{V_{jm}} is the value of attractiveness variable m for store j,
#' and \eqn{\alpha_m} is the corresponding exponent.
#'
#' The distance exponent is typically negative because probability decreases
#' with distance. Common values range from -1 to -3.
#'
#' @section Outputs:
#' **Market Share**: The weighted average probability across all demand points,
#' representing the proportion of total market potential captured by each store.
#'
#' **Expected Sales**: The sum of (probability × sales_potential) for each store,
#' representing the expected sales volume.
#'
#' **Entropy**: A measure of local competition. Higher entropy indicates more
#' competitive areas where multiple stores have similar probabilities.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create demand grid with spending potential
#' demand <- st_as_sf(expand.grid(x = 1:10, y = 1:10), coords = c("x", "y"))
#' demand$spending <- runif(100, 1000, 5000)
#'
#' # Existing stores with varying sizes (attractiveness)
#' stores <- st_as_sf(data.frame(
#'   id = c("Store_A", "Store_B", "Store_C"),
#'   sqft = c(50000, 25000, 75000),
#'   parking = c(200, 100, 300),
#'   x = c(2, 8, 5), y = c(2, 8, 5)
#' ), coords = c("x", "y"))
#'
#' # Single attractiveness variable
#' result <- huff(demand, stores,
#'                attractiveness_col = "sqft",
#'                distance_exponent = -2,
#'                sales_potential_col = "spending")
#'
#' # Multiple attractiveness variables with different exponents
#' # Composite: A = sqft^1.0 * parking^0.5
#' result_multi <- huff(demand, stores,
#'                      attractiveness_col = c("sqft", "parking"),
#'                      attractiveness_exponent = c(1.0, 0.5),
#'                      distance_exponent = -2,
#'                      sales_potential_col = "spending")
#'
#' # View market shares
#' result_multi$stores[, c("id", "sqft", "parking", ".market_share", ".expected_sales")]
#'
#' # Evaluate a new candidate store
#' candidate <- st_as_sf(data.frame(
#'   id = "New_Store", sqft = 40000, parking = 250, x = 3, y = 7
#' ), coords = c("x", "y"))
#'
#' all_stores <- rbind(stores, candidate)
#' result_with_candidate <- huff(demand, all_stores,
#'                               attractiveness_col = c("sqft", "parking"),
#'                               attractiveness_exponent = c(1.0, 0.5),
#'                               distance_exponent = -2,
#'                               sales_potential_col = "spending")
#'
#' # Compare market shares with and without candidate
#' result_with_candidate$stores[, c("id", ".market_share")]
#' }
#'
#' @references
#' Huff, D. L. (1963). A Probabilistic Analysis of Shopping Center Trade Areas.
#' Land Economics, 39(1), 81-90. \doi{10.2307/3144521}
#'
#' Huff, D. L. (1964). Defining and Estimating a Trading Area. Journal of
#' Marketing, 28(3), 34-38. \doi{10.1177/002224296402800307}
#'
#' @export
huff <- function(demand,
                 stores,
                 attractiveness_col,
                 attractiveness_exponent = 1,
                 distance_exponent = -1.5,
                 sales_potential_col = NULL,
                 cost_matrix = NULL,
                 distance_metric = "euclidean") {
  # Input validation
  if (!inherits(demand, "sf")) {
    stop("`demand` must be an sf object", call. = FALSE)
  }
  if (!inherits(stores, "sf")) {
    stop("`stores` must be an sf object", call. = FALSE)
  }

  # Validate attractiveness columns
  attractiveness_col <- as.character(attractiveness_col)
  missing_cols <- attractiveness_col[!attractiveness_col %in% names(stores)]
  if (length(missing_cols) > 0) {
    stop(paste0("Attractiveness column(s) not found in stores: ",
                paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Handle exponent recycling
  if (length(attractiveness_exponent) == 1) {
    attractiveness_exponent <- rep(attractiveness_exponent, length(attractiveness_col))
  } else if (length(attractiveness_exponent) != length(attractiveness_col)) {
    stop(paste0("attractiveness_exponent must be length 1 or same length as ",
                "attractiveness_col (", length(attractiveness_col), ")"), call. = FALSE)
  }

  # Compute composite attractiveness: A_j = prod(V_jm ^ alpha_m)
  attractiveness_adj <- rep(1.0, nrow(stores))
  for (k in seq_along(attractiveness_col)) {
    col_vals <- as.numeric(stores[[attractiveness_col[k]]])
    if (any(is.na(col_vals))) {
      stop(paste0("Attractiveness column '", attractiveness_col[k],
                  "' contains NA values"), call. = FALSE)
    }
    if (any(col_vals <= 0)) {
      stop(paste0("All values in '", attractiveness_col[k],
                  "' must be positive"), call. = FALSE)
    }
    attractiveness_adj <- attractiveness_adj * (col_vals ^ attractiveness_exponent[k])
  }

  # Get sales potential if provided
  sales_potential <- NULL
  if (!is.null(sales_potential_col)) {
    if (!sales_potential_col %in% names(demand)) {
      stop(paste0("Sales potential column '", sales_potential_col,
                  "' not found in demand"), call. = FALSE)
    }
    sales_potential <- as.numeric(demand[[sales_potential_col]])
    if (any(is.na(sales_potential))) {
      stop("Sales potential column contains NA values", call. = FALSE)
    }
  }

  # Compute distance matrix if needed
  if (is.null(cost_matrix)) {
    cost_matrix <- distance_matrix(demand, stores, type = distance_metric)
  }

  # Validate cost matrix
  if (any(is.na(cost_matrix))) {
    warning("cost_matrix contains NA values. Replacing with large value.")
    max_cost <- max(cost_matrix, na.rm = TRUE)
    cost_matrix[is.na(cost_matrix)] <- max_cost * 100
  }

  n_demand <- nrow(demand)
  n_stores <- nrow(stores)

  start_time <- Sys.time()

  # Call Rust implementation
  result <- rust_huff(
    cost_matrix,
    attractiveness_adj,
    distance_exponent,
    sales_potential
  )

  end_time <- Sys.time()

  # Check for errors
  if (!is.null(result$error)) {
    stop(result$error, call. = FALSE)
  }

  # Build probability matrix
  prob_matrix <- matrix(
    result$probabilities,
    nrow = n_demand,
    ncol = n_stores,
    byrow = TRUE
  )

  # Augment demand sf
  demand_result <- demand
  demand_result$.primary_store <- result$primary_store
  demand_result$.entropy <- result$entropy

  # Add probability columns for each store
  store_ids <- if ("id" %in% names(stores)) {
    stores$id
  } else {
    seq_len(n_stores)
  }

  for (j in seq_len(n_stores)) {
    col_name <- paste0(".prob_", store_ids[j])
    demand_result[[col_name]] <- prob_matrix[, j]
  }

  # Augment stores sf
  stores_result <- stores
  stores_result$.market_share <- result$market_shares
  stores_result$.expected_sales <- result$expected_sales

  # Build output
  output <- list(
    demand = demand_result,
    stores = stores_result,
    probability_matrix = prob_matrix
  )

  metadata <- list(
    algorithm = "huff",
    attractiveness_col = attractiveness_col,
    attractiveness_exponent = attractiveness_exponent,
    distance_exponent = distance_exponent,
    sales_potential_col = sales_potential_col,
    total_sales_potential = result$total_potential,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  )

  attr(output, "spopt") <- metadata
  class(output) <- c("spopt_huff", "list")

  output
}
