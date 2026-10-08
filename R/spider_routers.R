#' Road-following spider lines with r5r
#'
#' Creates a routing function for the `route_fun` argument of
#' [spider_lines()], using [r5r::detailed_itineraries()] so that spider lines
#' follow the street network instead of running straight. Use the same r5r
#' network you used to build your travel-time matrix.
#'
#' @param r5r_network An r5r network built with [r5r::build_network()].
#' @param mode Transport mode passed to [r5r::detailed_itineraries()].
#'   Defaults to `"CAR"`.
#' @param departure_datetime Departure time passed to
#'   [r5r::detailed_itineraries()]. Defaults to the current time.
#' @param ... Further arguments passed to [r5r::detailed_itineraries()], such
#'   as `max_trip_duration`.
#'
#' @return A function taking `(from, to)` that returns one MULTILINESTRING per
#'   origin-destination pair (EPSG:4326), with empty geometries for pairs r5r
#'   could not route. Pass it to [spider_lines()] as `route_fun`.
#'
#' @details
#' r5r returns each itinerary as one or more legs. The returned function keeps
#' the first itinerary option for each pair, orders its legs, and combines them
#' into a single geometry. Every pair is routed in one batched call, and
#' [spider_lines()] sends each distinct pair only once.
#'
#' @examples
#' \dontrun{
#' library(r5r)
#' options(java.parameters = "-Xmx4G")
#'
#' net <- build_network(data_path = "path/to/osm")
#' result <- p_median(demand, facilities, n_facilities = 5,
#'                    weight_col = "population", cost_matrix = ttm)
#'
#' lines <- spider_lines(
#'   result,
#'   route_fun = r5r_route_fun(net, mode = "CAR",
#'                             departure_datetime = as.POSIXct("2026-10-07 08:00"))
#' )
#' }
#'
#' @seealso [spider_lines()]
#'
#' @export
r5r_route_fun <- function(r5r_network,
                          mode = "CAR",
                          departure_datetime = Sys.time(),
                          ...) {
  if (!requireNamespace("r5r", quietly = TRUE)) {
    stop("The r5r package is required for r5r_route_fun()", call. = FALSE)
  }
  force(r5r_network)
  force(mode)
  force(departure_datetime)
  extra <- list(...)

  function(from, to) {
    args <- c(
      list(
        r5r_network,
        origins = from,
        destinations = to,
        mode = mode,
        departure_datetime = departure_datetime,
        shortest_path = TRUE,
        all_to_all = FALSE,
        progress = FALSE
      ),
      extra
    )
    itineraries <- do.call(r5r::detailed_itineraries, args)
    .assemble_itineraries(itineraries, from$id)
  }
}

# Combine r5r itinerary legs into one MULTILINESTRING per requested pair.
# `ids` gives the requested pair ids in order; pairs with no itinerary get an
# empty geometry.
.assemble_itineraries <- function(itineraries, ids) {
  n <- length(ids)
  empty <- sf::st_multilinestring()
  out <- rep(list(empty), n)

  if (!is.null(itineraries) && nrow(itineraries) > 0) {
    from_id <- as.character(itineraries$from_id)
    option <- itineraries$option
    segment <- itineraries$segment
    geoms <- sf::st_geometry(sf::st_transform(sf::st_zm(itineraries), 4326))

    for (key in unique(from_id)) {
      pos <- match(key, ids)
      if (is.na(pos)) next
      rows <- which(from_id == key)
      rows <- rows[option[rows] == min(option[rows])]
      rows <- rows[order(segment[rows])]
      parts <- lapply(rows, function(r) {
        g <- geoms[[r]]
        if (inherits(g, "MULTILINESTRING")) {
          lapply(unclass(g), function(m) m[, 1:2, drop = FALSE])
        } else {
          list(unclass(g)[, 1:2, drop = FALSE])
        }
      })
      parts <- unlist(parts, recursive = FALSE)
      parts <- parts[vapply(parts, nrow, integer(1)) >= 2]
      if (length(parts) > 0) out[[pos]] <- sf::st_multilinestring(parts)
    }
  }

  sf::st_sfc(out, crs = 4326)
}
