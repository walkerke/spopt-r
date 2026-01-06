#' Tarrant County Travel Time Matrix Example Data
#'
#' A dataset containing Census tract polygons, candidate facility points,
#' and a pre-computed driving travel-time matrix for Tarrant County, Texas.
#' Used for demonstrating facility location algorithms with real travel times.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{tracts}{An sf object with Census tract polygons including GEOID,
#'     NAME, population, and geometry columns.}
#'   \item{demand}{An sf object with tract centroid points (demand locations)
#'     in WGS84 coordinates.}
#'   \item{candidates}{An sf object with 30 randomly sampled candidate
#'     facility locations in WGS84 coordinates.}
#'   \item{matrix}{A numeric matrix of driving travel times in minutes.
#'     Rows correspond to demand points, columns to candidate facilities.
#'     Unreachable pairs are set to Inf.}
#' }
#'
#' @details
#' The travel-time matrix was generated using r5r with OpenStreetMap road
#' network data clipped to Tarrant County. Candidate facilities were sampled
#' using `set.seed(1983)` for reproducibility.
#'
#' @source
#' Census tract data from the American Community Survey via tidycensus.
#' Road network from OpenStreetMap via GeoFabrik.
#' Travel times computed with r5r.
#'
#' @examples
#' data(tarrant_travel_times)
#'
#' # Access components
#' tracts <- tarrant_travel_times$tracts
#' demand <- tarrant_travel_times$demand
#' candidates <- tarrant_travel_times$candidates
#' ttm <- tarrant_travel_times$matrix
#'
#' # Use with p_median
#' result <- p_median(
#'   demand = demand,
#'   facilities = candidates,
#'   n_facilities = 5,
#'   weight_col = "population",
#'   cost_matrix = ttm
#' )
"tarrant_travel_times"
