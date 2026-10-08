#' Spider lines for location-allocation results
#'
#' Builds "spider lines" (also called allocation lines or desire lines): one
#' line from each demand location to the facility that serves it. Works with
#' results from [p_median()], [p_center()], [cflp()], [lscp()], [mclp()], and
#' [huff()]. Lines are straight by default; supply `route_fun` to draw lines
#' that follow a road network instead.
#'
#' @param result A result from [p_median()], [p_center()], [cflp()], [lscp()],
#'   [mclp()], or [huff()], with demand and facility rows in the order the
#'   solver returned them.
#' @param allocations Which links to draw. `"primary"` (default) draws at most
#'   one line per demand row: the assigned facility (`.facility`) or the
#'   highest-probability store (`.primary_store`). `"all"` also draws split
#'   allocations from [cflp()] and every store with non-zero probability from
#'   [huff()]; for the other solvers it is the same as `"primary"`.
#' @param min_share Drop links whose `share` is less than or equal to this
#'   value. Applies to primary and additional links alike; shares are never
#'   renormalized. Must be in `[0, 1)`. Most useful for [huff()] results
#'   with `allocations = "all"`.
#' @param cost_matrix Optional demand x facility cost matrix, normally the one
#'   passed to the solver. Not needed for primary links, which take their cost
#'   from the solver's `.cost` column. If supplied, it takes precedence and
#'   provides `cost` for every link, including additional [cflp()] and [huff()]
#'   links.
#' @param anchor How polygon (or multipoint) geometries are reduced to line
#'   endpoints: `"centroid"` (default, matching [distance_matrix()]) or
#'   `"point_on_surface"` (guaranteed to fall inside the polygon).
#' @param route_fun Optional function that returns network-following
#'   geometries. It is called once as `route_fun(from, to)`, where `from` and
#'   `to` are sf POINT data frames in EPSG:4326 with the same number of rows and
#'   an `id` column; row *i* of `from` should be routed to row *i* of `to`.
#'   Identical coordinate pairs are sent only once. It must return an sf or
#'   sfc object with one LINESTRING or MULTILINESTRING per row, in input order,
#'   with a known CRS. Return an empty geometry for a pair that cannot be
#'   routed; those links fall back to straight lines. [r5r_route_fun()]
#'   builds such a function from an r5r network.
#'
#' @return An sf object in the CRS of `result$demand`, one row per link, with
#'   columns:
#'   \itemize{
#'     \item `demand_index`: row in `result$demand`
#'     \item `facility_index`: row in `result$facilities` (or `result$stores`)
#'     \item `share`: 1 for whole assignments; the allocation fraction for
#'       [cflp()]; the choice probability for [huff()]
#'     \item `weight`: demand weight multiplied by `share`, when the solver
#'       recorded a weight column
#'     \item `cost`: link cost from `cost_matrix`, or from the solver's `.cost`
#'       column for primary links; omitted when unavailable
#'     \item `length`: length of the drawn line
#'     \item `routed`: only when `route_fun` is used; `FALSE` where the router
#'       failed and a straight line was drawn instead
#'   }
#'   Geometry is LINESTRING, or MULTILINESTRING when `route_fun` is used.
#'   Counts of unassigned demand, filtered links, and routing failures are
#'   stored in the "spopt" attribute.
#'
#' @details
#' Indices and matrices in spopt results are positional, so do not reorder or
#' filter `result$demand` or `result$facilities` before calling
#' `spider_lines()`. Join attributes onto the output with `demand_index` and
#' `facility_index` instead.
#'
#' Demand that the solver left unassigned (uncovered points in [lscp()] or
#' [mclp()], or an infeasible [p_center()]) produces no lines.
#'
#' For data in longitude/latitude, lines are drawn as two-vertex segments.
#' Use [sf::st_segmentize()] on the output if you need great-circle arcs.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' demand <- st_as_sf(data.frame(
#'   x = runif(100), y = runif(100), population = rpois(100, 500)
#' ), coords = c("x", "y"))
#' facilities <- st_as_sf(data.frame(x = runif(20), y = runif(20)),
#'                        coords = c("x", "y"))
#'
#' result <- p_median(demand, facilities, n_facilities = 5,
#'                    weight_col = "population")
#' lines <- spider_lines(result)
#'
#' plot(st_geometry(lines), col = lines$facility_index)
#' plot(st_geometry(result$facilities[result$facilities$.selected, ]),
#'      pch = 19, cex = 1.5, add = TRUE)
#' }
#'
#' \dontrun{
#' # Lines that follow the road network, using a local OSRM server
#' osrm_lines <- function(from, to) {
#'   geoms <- lapply(seq_len(nrow(from)), function(i) {
#'     r <- tryCatch(osrm::osrmRoute(src = from[i, ], dst = to[i, ]),
#'                   error = function(e) NULL)
#'     if (is.null(r)) sf::st_linestring() else sf::st_geometry(r)[[1]]
#'   })
#'   sf::st_sfc(geoms, crs = 4326)
#' }
#'
#' routed <- spider_lines(result, route_fun = osrm_lines)
#' }
#'
#' @seealso [r5r_route_fun()] for road-following lines with r5r;
#'   [p_median()], [mclp()], [cflp()], [huff()]
#'
#' @export
spider_lines <- function(result,
                         allocations = c("primary", "all"),
                         min_share = 0,
                         cost_matrix = NULL,
                         anchor = c("centroid", "point_on_surface"),
                         route_fun = NULL) {
  allocations <- match.arg(allocations)
  anchor <- match.arg(anchor)

  if (!inherits(result, c("spopt_locate", "spopt_huff"))) {
    stop(
      "`result` must come from p_median(), p_center(), cflp(), lscp(), ",
      "mclp(), or huff()",
      call. = FALSE
    )
  }
  is_huff <- inherits(result, "spopt_huff")
  demand <- result$demand
  facilities <- if (is_huff) result$stores else result$facilities
  meta <- attr(result, "spopt")
  if (is.null(meta)) meta <- list()
  algorithm <- if (!is.null(meta$algorithm)) meta$algorithm else class(result)[1]

  if (algorithm %in% c("frlm", "p_dispersion")) {
    stop("spider_lines() does not support ", algorithm, "() results",
         call. = FALSE)
  }

  if (!is.numeric(min_share) || length(min_share) != 1 || is.na(min_share) ||
      min_share < 0 || min_share >= 1) {
    stop("`min_share` must be a single number in [0, 1)", call. = FALSE)
  }

  n_demand <- nrow(demand)
  n_fac <- nrow(facilities)

  if (!is.null(cost_matrix)) {
    cost_matrix <- as.matrix(cost_matrix)
    if (!identical(dim(cost_matrix), c(n_demand, n_fac))) {
      stop(sprintf(
        "`cost_matrix` must be %d x %d (demand x facilities), not %d x %d",
        n_demand, n_fac, nrow(cost_matrix), ncol(cost_matrix)
      ), call. = FALSE)
    }
  }

  if (!is.null(route_fun) && !is.function(route_fun)) {
    stop("`route_fun` must be a function or NULL", call. = FALSE)
  }

  # CRS handling
  crs <- sf::st_crs(demand)
  crs_fac <- sf::st_crs(facilities)
  if (is.na(crs) != is.na(crs_fac)) {
    stop("`demand` and `facilities` must both have a CRS or both lack one",
         call. = FALSE)
  }
  if (!is.na(crs) && crs != crs_fac) {
    facilities <- sf::st_transform(facilities, crs)
  }
  if (is.na(crs) && !is.null(route_fun)) {
    stop("`route_fun` requires `demand` and `facilities` to have a known CRS",
         call. = FALSE)
  }

  # Demand-facility pairs
  pr <- .spider_pairs(result, demand, n_fac, algorithm, is_huff, meta,
                      allocations, min_share)
  pairs <- pr$pairs

  # Endpoints
  d_pts <- .spider_anchor(demand, anchor)
  f_pts <- .spider_anchor(facilities, anchor)
  d_xy <- .spider_coords(d_pts)
  f_xy <- .spider_coords(f_pts)

  from_xy <- d_xy[pairs$d, , drop = FALSE]
  to_xy <- f_xy[pairs$f, , drop = FALSE]
  # Empty, invalid, and non-finite (NA/Inf) endpoints are all dropped
  bad_geom <- rowSums(!is.finite(from_xy)) > 0 | rowSums(!is.finite(to_xy)) > 0
  n_invalid_geometry <- sum(bad_geom)
  if (n_invalid_geometry > 0) {
    warning(sprintf(
      "Dropped %d link(s) with an empty or invalid demand or facility geometry",
      n_invalid_geometry
    ), call. = FALSE)
    pairs <- pairs[!bad_geom, , drop = FALSE]
    from_xy <- from_xy[!bad_geom, , drop = FALSE]
    to_xy <- to_xy[!bad_geom, , drop = FALSE]
  }

  n_links <- nrow(pairs)
  if (n_links == 0) {
    warning("No valid demand-facility links; returning an empty result",
            call. = FALSE)
  }
  if (n_links > 1e6) {
    warning(sprintf(
      "Building %s spider lines; consider raising `min_share`",
      format(n_links, big.mark = ",")
    ), call. = FALSE)
  }

  # Geometry
  geom <- .spider_straight(from_xy, to_xy, crs)
  routed <- NULL
  n_route_requests <- 0L
  n_route_failures <- 0L
  if (!is.null(route_fun) && n_links > 0) {
    rt <- .spider_route(from_xy, to_xy, geom, crs, route_fun)
    geom <- rt$geom
    routed <- rt$routed
    n_route_requests <- rt$n_requests
    n_route_failures <- rt$n_failures
  } else if (!is.null(route_fun)) {
    geom <- sf::st_cast(geom, "MULTILINESTRING")
    routed <- logical(0)
  }

  # Attributes
  out <- data.frame(
    demand_index = pairs$d,
    facility_index = pairs$f,
    share = pairs$share
  )

  weight_col <- if (is_huff) meta$sales_potential_col else meta$weight_col
  if (!is.null(weight_col) && weight_col %in% names(demand)) {
    w <- as.numeric(sf::st_drop_geometry(demand)[[weight_col]])
    out$weight <- w[pairs$d] * pairs$share
  }

  if (!is.null(cost_matrix)) {
    out$cost <- as.numeric(cost_matrix[cbind(pairs$d, pairs$f)])
  } else if (any(!pairs$primary)) {
    message(
      "Costs for non-primary links need `cost_matrix`; `cost` column omitted"
    )
  } else if (".cost" %in% names(demand)) {
    out$cost <- as.numeric(demand$.cost[pairs$d])
  }

  out$length <- sf::st_length(geom)
  if (!is.null(routed)) out$routed <- routed

  out <- sf::st_sf(out, geometry = geom)

  attr(out, "spopt") <- list(
    algorithm = algorithm,
    allocations = allocations,
    min_share = min_share,
    n_lines = n_links,
    n_unassigned = pr$n_unassigned,
    n_filtered = pr$n_filtered,
    n_invalid_geometry = n_invalid_geometry,
    routed = !is.null(route_fun),
    n_route_requests = n_route_requests,
    n_route_failures = n_route_failures
  )

  out
}

# Build the table of demand-facility links: d, f, share, primary
.spider_pairs <- function(result, demand, n_fac, algorithm, is_huff, meta,
                          allocations, min_share) {
  tol <- 1e-6

  assign <- if (is_huff) demand$.primary_store else demand$.facility
  if (is.null(assign)) {
    stop("`result$demand` has no assignment column (",
         if (is_huff) ".primary_store" else ".facility", ")", call. = FALSE)
  }
  assign <- as.numeric(assign)
  valid <- !is.na(assign) & assign == round(assign) &
    assign >= 1 & assign <= n_fac
  odd <- !valid & !is.na(assign) & assign != 0
  if (any(odd)) {
    warning(sprintf(
      paste0("%d assignment(s) point outside the facility table and were ",
             "skipped; was `result` reordered or filtered?"),
      sum(odd)
    ), call. = FALSE)
  }
  n_unassigned <- sum(!valid)

  share_matrix <- if (is_huff) {
    result$probability_matrix
  } else if (identical(algorithm, "cflp")) {
    meta$allocation_matrix
  } else {
    NULL
  }

  d_primary <- which(valid)
  f_primary <- as.integer(assign[valid])

  # Shares at or below this are not links at all (CFLP solver residue)
  threshold <- if (is_huff) 0 else tol

  if (allocations == "all" && !is.null(share_matrix)) {
    idx <- which(share_matrix > threshold, arr.ind = TRUE)
    in_bounds <- idx[, 1] <= nrow(demand) & idx[, 2] <= n_fac
    if (any(!in_bounds)) {
      warning(sprintf(
        paste0("%d allocation(s) point outside the demand or facility table ",
               "and were skipped; was `result` reordered or filtered?"),
        sum(!in_bounds)
      ), call. = FALSE)
      idx <- idx[in_bounds, , drop = FALSE]
    }
    idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
    d <- as.integer(idx[, 1])
    f <- as.integer(idx[, 2])
    share <- share_matrix[idx]
    primary <- valid[d] & f == assign[d]
  } else {
    if (allocations == "all" && identical(algorithm, "cflp")) {
      warning("No `allocation_matrix` found; drawing primary links only",
              call. = FALSE)
    }
    d <- d_primary
    f <- f_primary
    share <- if (!is.null(share_matrix)) {
      share_matrix[cbind(d, f)]
    } else {
      rep(1, length(d))
    }
    primary <- rep(TRUE, length(d))
    if (!is.null(share_matrix)) {
      real <- share > threshold
      d <- d[real]
      f <- f[real]
      share <- share[real]
      primary <- primary[real]
    }
  }

  keep <- share > min_share
  list(
    pairs = data.frame(
      d = d[keep],
      f = f[keep],
      share = as.numeric(share[keep]),
      primary = primary[keep]
    ),
    n_unassigned = n_unassigned,
    n_filtered = sum(!keep)
  )
}

# Reduce geometries to points for line endpoints
.spider_anchor <- function(x, anchor) {
  g <- sf::st_zm(sf::st_geometry(x))
  types <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
  empty <- sf::st_is_empty(g)
  allowed <- c("POINT", "MULTIPOINT", "POLYGON", "MULTIPOLYGON")
  bad <- !types %in% allowed & !empty
  if (any(bad)) {
    stop("Unsupported geometry type(s) for spider lines: ",
         paste(unique(types[bad]), collapse = ", "), call. = FALSE)
  }
  # Invalid polygons are treated like empty geometries (link dropped)
  poly <- types %in% c("POLYGON", "MULTIPOLYGON") & !empty
  if (any(poly)) {
    ok <- suppressWarnings(sf::st_is_valid(g[poly]))
    bad_poly <- which(poly)[is.na(ok) | !ok]
    if (length(bad_poly) > 0) {
      g[bad_poly] <- sf::st_sfc(sf::st_point(), crs = sf::st_crs(g))
      empty[bad_poly] <- TRUE
    }
  }

  needs <- types != "POINT" & !empty
  if (any(needs)) {
    reduced <- if (anchor == "centroid") {
      sf::st_centroid(g[needs])
    } else {
      sf::st_point_on_surface(g[needs])
    }
    g[needs] <- reduced
  }
  g
}

# XY matrix for a point sfc; NA rows for empty points
.spider_coords <- function(pts) {
  xy <- matrix(NA_real_, nrow = length(pts), ncol = 2)
  ok <- !sf::st_is_empty(pts)
  if (any(ok)) {
    xy[ok, ] <- sf::st_coordinates(pts[ok])[, 1:2, drop = FALSE]
  }
  xy
}

# Two-vertex LINESTRINGs, built directly for speed
.spider_straight <- function(from_xy, to_xy, crs) {
  n <- nrow(from_xy)
  lines <- vector("list", n)
  for (i in seq_len(n)) {
    lines[[i]] <- structure(
      matrix(c(from_xy[i, 1], to_xy[i, 1], from_xy[i, 2], to_xy[i, 2]),
             ncol = 2),
      class = c("XY", "LINESTRING", "sfg")
    )
  }
  if (n == 0) {
    return(sf::st_sfc(sf::st_linestring(), crs = crs)[0])
  }
  sf::st_sfc(lines, crs = crs)
}

# Call the user's router once on unique pairs and map results back
.spider_route <- function(from_xy, to_xy, straight, crs, route_fun) {
  key <- paste(from_xy[, 1], from_xy[, 2], to_xy[, 1], to_xy[, 2])
  uniq <- unique(key)
  first <- match(uniq, key)
  map <- match(key, uniq)
  n_req <- length(uniq)

  to_points <- function(xy) {
    pts <- sf::st_as_sf(as.data.frame(xy), coords = c(1, 2), crs = crs)
    pts <- sf::st_transform(pts, 4326)
    sf::st_sf(id = as.character(seq_len(nrow(xy))),
              geometry = sf::st_geometry(pts))
  }
  from <- to_points(from_xy[first, , drop = FALSE])
  to <- to_points(to_xy[first, , drop = FALSE])

  res <- route_fun(from, to)

  if (inherits(res, "sf")) {
    geom <- sf::st_geometry(res)
  } else if (inherits(res, "sfc")) {
    geom <- res
  } else {
    stop("`route_fun` must return an sf or sfc object", call. = FALSE)
  }
  if (length(geom) != n_req) {
    stop(sprintf(
      "`route_fun` returned %d geometries for %d requested routes",
      length(geom), n_req
    ), call. = FALSE)
  }
  if (is.na(sf::st_crs(geom))) {
    stop("`route_fun` must return geometries with a known CRS", call. = FALSE)
  }
  failed <- sf::st_is_empty(geom)
  types <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
  bad <- !failed & !types %in% c("LINESTRING", "MULTILINESTRING")
  if (any(bad)) {
    stop("`route_fun` must return LINESTRING or MULTILINESTRING geometries, ",
         "not ", paste(unique(types[bad]), collapse = ", "), call. = FALSE)
  }

  link_failed <- failed[map]
  out <- sf::st_cast(straight, "MULTILINESTRING")
  ok <- which(!link_failed)
  if (length(ok) > 0) {
    good <- sf::st_transform(geom[!failed], crs)
    good <- sf::st_cast(sf::st_zm(good), "MULTILINESTRING")
    out[ok] <- good[match(map[ok], which(!failed))]
  }

  n_fail <- sum(failed)
  if (n_fail > 0) {
    warning(sprintf(
      "%d of %d route(s) failed; drew straight lines for %d link(s)",
      n_fail, n_req, sum(link_failed)
    ), call. = FALSE)
  }

  list(geom = out, routed = !link_failed, n_requests = n_req,
       n_failures = n_fail)
}
