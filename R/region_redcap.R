#' REDCAP Regionalization
#'
#' Performs regionalization with the REDCAP algorithm (REgionalization with
#' Dynamically Constrained Agglomerative clustering and Partitioning). A
#' contiguity-constrained hierarchical clustering builds a spanning tree over
#' the areas; the tree is then partitioned by removing the edges whose removal
#' lowers within-region heterogeneity the most, as in [skater()].
#'
#' @param data An sf object with polygon or point geometries.
#' @param attrs Character vector of column names to use for clustering
#'   (e.g., `c("var1", "var2")`). If NULL, uses all numeric columns.
#' @param n_regions Integer. Number of regions (clusters) to create.
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
#' @param method Character. Linkage used to build the tree; names follow GeoDa
#'   and rgeoda:
#'   \itemize{
#'     \item `"fullorder-averagelinkage"` (default): mean attribute distance over
#'       all pairs of areas between two clusters
#'     \item `"fullorder-singlelinkage"`: minimum distance over all pairs
#'     \item `"fullorder-completelinkage"`: maximum distance over all pairs
#'     \item `"fullorder-wardlinkage"`: Ward's minimum-variance criterion
#'     \item `"firstorder-singlelinkage"`: shortest contiguity edge only; the
#'       tree is the minimum spanning tree, so this equals [skater()]
#'   }
#'   "Full order" means the cluster distance considers every pair of areas in
#'   the two clusters, not only adjacent pairs; only spatially contiguous
#'   clusters may merge in all variants.
#' @param floor Optional. Column name specifying a floor constraint variable.
#' @param floor_value Numeric. Minimum sum of floor variable required per region.
#'   Only used if `floor` is specified.
#' @param scale Logical. If TRUE (default), standardize attributes before clustering.
#' @param verbose Logical. Print progress messages.
#'
#' @return An sf object with a `.region` column containing cluster assignments.
#'   Metadata is stored in the "spopt" attribute, including:
#'   \itemize{
#'     \item algorithm: "redcap"
#'     \item method: The linkage used
#'     \item n_regions: Number of regions created
#'     \item objective: Total within-region sum of squared deviations
#'     \item solve_time: Time to solve in seconds
#'   }
#'
#' @details
#' REDCAP (Guo 2008) differs from SKATER only in the tree it cuts. SKATER cuts
#' the minimum spanning tree of the contiguity graph, whose edges reflect the
#' similarity of neighbouring pairs of areas. REDCAP first runs a
#' contiguity-constrained agglomerative clustering and records, for each merge,
#' the most similar pair of adjacent areas between the two clusters being
#' joined; those edges form the tree. Because merges are chosen by cluster-level
#' similarity, the tree tends to respect larger homogeneous zones, and the
#' resulting regions are usually more homogeneous than SKATER's at the same k.
#'
#' With `method = "fullorder-wardlinkage"`, the agglomeration is the same as
#' [ward_spatial()]; cutting the recorded tree rather than reading regions off
#' the merge sequence lets each cut re-optimize where a region splits. The
#' full-order single, average and complete linkages keep an n x n matrix of
#' cluster distances (8 n^2 bytes), as GeoDa does; Ward and first-order single
#' linkage run in O(n) memory.
#'
#' Results match `rgeoda::redcap()` for the same method name, with one
#' exception: rgeoda's average linkage runs on squared Euclidean distances (a
#' line in libgeoda's distance-matrix conversion overwrites the square root it
#' has just taken), while this implementation uses Euclidean distances as in
#' Guo (2008). Single and complete linkage are unaffected because minimum and
#' maximum do not depend on the square root. Reported upstream as
#' <https://github.com/GeoDaCenter/geoda/issues/2613>.
#'
#' @examples
#' \donttest{
#' library(sf)
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"))
#'
#' result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 5)
#' attr(result, "spopt")$objective
#'
#' # Ward linkage is often the most homogeneous
#' result_ward <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 5,
#'                       method = "fullorder-wardlinkage")
#'
#' # Compare with SKATER and constrained Ward at the same k
#' attr(skater(nc, attrs = c("SID74", "SID79"), n_regions = 5), "spopt")$objective
#' attr(ward_spatial(nc, attrs = c("SID74", "SID79"), n_regions = 5), "spopt")$objective
#'
#' # With a floor: each region must have at least 100,000 births
#' result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 5,
#'                  floor = "BIR74", floor_value = 100000)
#' plot(result[".region"])
#' }
#'
#' @references
#' Guo, D. (2008). Regionalization with dynamically constrained agglomerative
#' clustering and partitioning (REDCAP). International Journal of Geographical
#' Information Science, 22(7), 801-823. \doi{10.1080/13658810701674970}
#'
#' @seealso [skater()], [ward_spatial()]
#' @export
redcap <- function(data,
                   attrs = NULL,
                   n_regions,
                   weights = "queen",
                   bridge_islands = FALSE,
                   method = c("fullorder-averagelinkage", "fullorder-singlelinkage",
                              "fullorder-completelinkage", "fullorder-wardlinkage",
                              "firstorder-singlelinkage"),
                   floor = NULL,
                   floor_value = 0,
                   scale = TRUE,
                   verbose = FALSE) {
  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object", call. = FALSE)
  }
  if (!is.numeric(n_regions) || n_regions < 2) {
    stop("`n_regions` must be an integer >= 2", call. = FALSE)
  }
  method <- match.arg(method)

  check_cols <- c(
    if (!is.null(attrs)) attrs else character(0),
    if (!is.null(floor)) floor else character(0)
  )
  validated <- validate_regionalization_data(data, check_cols, call_name = "redcap")
  data <- validated$data

  n <- nrow(data)
  if (n_regions >= n) {
    stop("`n_regions` must be less than number of observations", call. = FALSE)
  }

  attr_matrix <- extract_attrs(data, attrs)
  if (scale) {
    attr_matrix <- scale(attr_matrix)
  }

  nb <- prepare_weights(data, weights, bridge_islands = bridge_islands, call_name = "redcap")
  adj <- nb_to_adj_indices(nb)

  floor_var <- NULL
  if (!is.null(floor)) {
    if (!floor %in% names(data)) {
      stop(paste0("Floor variable '", floor, "' not found in data"), call. = FALSE)
    }
    floor_var <- as.numeric(data[[floor]])
  }

  if (verbose) {
    message(sprintf("REDCAP: n=%d, n_regions=%d, method=%s, attrs=%d",
                    n, n_regions, method, ncol(attr_matrix)))
  }

  start_time <- Sys.time()
  result_list <- rust_redcap(
    attrs = attr_matrix,
    n_regions = as.integer(n_regions),
    adj_i = adj$i,
    adj_j = adj$j,
    method = method,
    floor_var = floor_var,
    floor_value = as.numeric(floor_value)
  )
  end_time <- Sys.time()

  labels <- result_list$labels
  result <- data
  result$.region <- as.character(labels)

  actual_n_regions <- result_list$n_regions
  if (actual_n_regions < n_regions) {
    warning(
      "Could only create ", actual_n_regions, " regions (requested ", n_regions, "). ",
      if (!is.null(floor)) "Floor constraint may be too restrictive." else "Graph may be disconnected.",
      call. = FALSE
    )
  }

  if (verbose) {
    message(sprintf("  Result: %d regions, objective=%.4f, time=%.3fs",
                    actual_n_regions, result_list$objective,
                    as.numeric(difftime(end_time, start_time, units = "secs"))))
  }

  metadata <- list(
    algorithm = "redcap",
    method = method,
    n_regions = actual_n_regions,
    objective = result_list$objective,
    floor = floor,
    floor_value = floor_value,
    solve_time = as.numeric(difftime(end_time, start_time, units = "secs")),
    scaled = scale
  )

  attach_spopt_metadata(result, metadata)
}
