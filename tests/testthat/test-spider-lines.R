skip_if_not_installed("sf")
library(sf)

make_spider_data <- function(n_demand = 30, n_fac = 8, seed = 42, crs = NA) {
  set.seed(seed)
  demand <- st_as_sf(
    data.frame(
      x = runif(n_demand, 0, 10),
      y = runif(n_demand, 0, 10),
      pop = rpois(n_demand, 50) + 1L
    ),
    coords = c("x", "y"), crs = crs
  )
  facilities <- st_as_sf(
    data.frame(x = runif(n_fac, 0, 10), y = runif(n_fac, 0, 10),
               sqft = runif(n_fac, 1, 5)),
    coords = c("x", "y"), crs = crs
  )
  list(demand = demand, facilities = facilities)
}

line_ends <- function(lines) {
  xy <- st_coordinates(lines)
  l1 <- xy[, "L1"]
  list(
    start = xy[!duplicated(l1), c("X", "Y"), drop = FALSE],
    end = xy[!duplicated(l1, fromLast = TRUE), c("X", "Y"), drop = FALSE]
  )
}

# ---------------------------------------------------------------------------
# .cost on solver output
# ---------------------------------------------------------------------------

test_that("solvers add .cost from the original cost matrix", {
  d <- make_spider_data()
  cm <- distance_matrix(d$demand, d$facilities)

  pm <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  i <- seq_len(nrow(d$demand))
  expect_equal(pm$demand$.cost, cm[cbind(i, pm$demand$.facility)])
  expect_identical(attr(pm, "spopt")$weight_col, "pop")

  pc <- p_center(d$demand, d$facilities, n_facilities = 3)
  expect_equal(pc$demand$.cost, cm[cbind(i, pc$demand$.facility)])

  mc <- mclp(d$demand, d$facilities, service_radius = 2, n_facilities = 2,
             weight_col = "pop")
  covered <- !is.na(mc$demand$.facility)
  expect_true(any(!covered))
  expect_true(all(is.na(mc$demand$.cost[!covered])))
  expect_equal(mc$demand$.cost[covered],
               cm[cbind(which(covered), mc$demand$.facility[covered])])

  ls <- lscp(d$demand, d$facilities, service_radius = 8)
  ok <- !is.na(ls$demand$.facility)
  expect_equal(ls$demand$.cost[ok],
               cm[cbind(which(ok), ls$demand$.facility[ok])])

  hf <- huff(d$demand, d$facilities, attractiveness_col = "sqft")
  expect_equal(hf$demand$.cost, cm[cbind(i, hf$demand$.primary_store)])
})

test_that(".cost reports original values, not NA/Inf penalties", {
  d <- make_spider_data(n_demand = 6, n_fac = 3)
  cm <- distance_matrix(d$demand, d$facilities)
  cm[1, ] <- NA  # demand 1 cannot reach anything
  pm <- suppressWarnings(
    p_median(d$demand, d$facilities, n_facilities = 2, weight_col = "pop",
             cost_matrix = cm)
  )
  expect_true(is.na(pm$demand$.cost[1]))
  expect_true(all(!is.na(pm$demand$.cost[-1])))
})

# ---------------------------------------------------------------------------
# Basic structure per solver
# ---------------------------------------------------------------------------

test_that("p_median spider lines connect demand to assigned facility", {
  d <- make_spider_data()
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  lines <- spider_lines(res)

  expect_s3_class(lines, "sf")
  expect_equal(nrow(lines), nrow(d$demand))
  expect_equal(lines$facility_index, res$demand$.facility[lines$demand_index])
  expect_true(all(lines$share == 1))
  expect_equal(lines$weight, d$demand$pop[lines$demand_index])
  expect_equal(lines$cost, res$demand$.cost[lines$demand_index])
  expect_true(all(st_geometry_type(lines) == "LINESTRING"))

  ends <- line_ends(lines)
  expect_equal(unname(ends$start),
               unname(st_coordinates(d$demand)[lines$demand_index, ]))
  expect_equal(unname(ends$end),
               unname(st_coordinates(d$facilities)[lines$facility_index, ]))

  meta <- attr(lines, "spopt")
  expect_equal(meta$algorithm, "p_median")
  expect_equal(meta$n_lines, nrow(lines))
  expect_equal(meta$n_unassigned, 0)
})

test_that("p_center and lscp results are supported; no weight column", {
  d <- make_spider_data()
  pc <- spider_lines(p_center(d$demand, d$facilities, n_facilities = 3))
  expect_equal(nrow(pc), nrow(d$demand))
  expect_false("weight" %in% names(pc))

  ls <- spider_lines(lscp(d$demand, d$facilities, service_radius = 8))
  expect_false("weight" %in% names(ls))
  expect_true(nrow(ls) > 0)
})

test_that("uncovered MCLP demand produces no lines", {
  d <- make_spider_data()
  res <- mclp(d$demand, d$facilities, service_radius = 2, n_facilities = 2,
              weight_col = "pop")
  lines <- spider_lines(res)
  n_cov <- sum(res$demand$.covered)
  expect_equal(nrow(lines), n_cov)
  expect_equal(attr(lines, "spopt")$n_unassigned, nrow(d$demand) - n_cov)
})

test_that("invalid assignments (0, NA, out of range) are skipped", {
  d <- make_spider_data()
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  res$demand$.facility[1:3] <- 0L   # e.g. infeasible p_center
  res$demand$.facility[4] <- NA
  expect_silent(lines <- spider_lines(res))
  expect_equal(nrow(lines), nrow(d$demand) - 4)
  expect_equal(attr(lines, "spopt")$n_unassigned, 4)

  res$demand$.facility[5] <- 999L
  expect_warning(lines <- spider_lines(res), "outside the facility table")
  expect_equal(attr(lines, "spopt")$n_unassigned, 5)
})

test_that("no valid assignments returns an empty sf with the full schema", {
  d <- make_spider_data()
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  res$demand$.facility <- 0L
  expect_warning(lines <- spider_lines(res), "No valid")
  expect_equal(nrow(lines), 0)
  expect_true(all(c("demand_index", "facility_index", "share", "weight",
                    "cost", "length") %in% names(lines)))
})

# ---------------------------------------------------------------------------
# CFLP splits
# ---------------------------------------------------------------------------

make_split_data <- function() {
  demand <- st_as_sf(
    data.frame(x = c(0, 10, 5), y = c(0, 0, 3), pop = c(100, 20, 10)),
    coords = c("x", "y")
  )
  facilities <- st_as_sf(
    data.frame(x = c(1, 9), y = c(0, 0), cap = c(70, 70)),
    coords = c("x", "y")
  )
  cflp(demand, facilities, n_facilities = 2, weight_col = "pop",
       capacity_col = "cap")
}

test_that("cflp 'all' draws split allocations with unrenormalized shares", {
  res <- make_split_data()
  expect_gt(attr(res, "spopt")$n_split_demand, 0)

  prim <- spider_lines(res)
  expect_equal(nrow(prim), 3)
  expect_true(any(prim$share < 1))

  expect_message(all_lines <- spider_lines(res, allocations = "all"),
                 "cost_matrix")
  expect_gt(nrow(all_lines), 3)
  expect_false("cost" %in% names(all_lines))
  sums <- tapply(all_lines$share, all_lines$demand_index, sum)
  expect_equal(as.numeric(sums), rep(1, 3), tolerance = 1e-6)
  expect_equal(sum(all_lines$weight), sum(c(100, 20, 10)), tolerance = 1e-6)
})

test_that("cflp residue below tolerance is excluded", {
  res <- make_split_data()
  alloc <- attr(res, "spopt")$allocation_matrix
  alloc[alloc == 0] <- 1e-9
  attr(res, "spopt")$allocation_matrix <- alloc
  lines <- suppressMessages(spider_lines(res, allocations = "all"))
  expect_true(all(lines$share > 1e-6))
})

test_that("min_share filters primary links and counts them separately", {
  res <- make_split_data()
  prim_share <- spider_lines(res)$share
  cutoff <- min(prim_share[prim_share < 1])
  lines <- spider_lines(res, min_share = cutoff)
  meta <- attr(lines, "spopt")
  expect_equal(meta$n_filtered, sum(prim_share <= cutoff))
  expect_equal(meta$n_unassigned, 0)
  expect_error(spider_lines(res, min_share = 1), "min_share")
  expect_error(spider_lines(res, min_share = -0.1), "min_share")
})

# ---------------------------------------------------------------------------
# Huff
# ---------------------------------------------------------------------------

test_that("huff primary and all links use matrix probabilities", {
  d <- make_spider_data(n_demand = 20, n_fac = 4)
  d$demand$spend <- d$demand$pop * 10
  cm <- distance_matrix(d$demand, d$facilities)
  cm[, 4] <- Inf  # store 4 is unreachable -> zero probability
  res <- huff(d$demand, d$facilities, attractiveness_col = "sqft",
              sales_potential_col = "spend", cost_matrix = cm)
  pm <- res$probability_matrix

  prim <- spider_lines(res)
  expect_equal(prim$share, pm[cbind(prim$demand_index, prim$facility_index)])
  expect_equal(prim$weight, d$demand$spend[prim$demand_index] * prim$share)

  all_lines <- suppressMessages(spider_lines(res, allocations = "all"))
  expect_equal(nrow(all_lines), sum(pm > 0))
  expect_false(4 %in% all_lines$facility_index)

  filt <- suppressMessages(
    spider_lines(res, allocations = "all", min_share = 0.2)
  )
  expect_equal(nrow(filt), sum(pm > 0.2))
  expect_equal(attr(filt, "spopt")$n_filtered, sum(pm > 0) - sum(pm > 0.2))

  with_cost <- spider_lines(res, allocations = "all", cost_matrix = cm)
  expect_equal(with_cost$cost,
               cm[cbind(with_cost$demand_index, with_cost$facility_index)])
})

# ---------------------------------------------------------------------------
# Cost resolution
# ---------------------------------------------------------------------------

test_that("automatic and explicit costs agree; cost_matrix is validated", {
  d <- make_spider_data()
  cm <- distance_matrix(d$demand, d$facilities) * 2
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop",
                  cost_matrix = cm)
  auto <- spider_lines(res)
  explicit <- spider_lines(res, cost_matrix = cm)
  expect_equal(auto$cost, explicit$cost)
  expect_error(spider_lines(res, cost_matrix = cm[-1, ]), "must be")
})

test_that("results without .cost or weight_col still work", {
  d <- make_spider_data()
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  res$demand$.cost <- NULL
  attr(res, "spopt")$weight_col <- NULL
  lines <- spider_lines(res)
  expect_false(any(c("cost", "weight") %in% names(lines)))
  expect_equal(nrow(lines), nrow(d$demand))
})

# ---------------------------------------------------------------------------
# Geometry and CRS
# ---------------------------------------------------------------------------

test_that("polygon demand uses centroid or point-on-surface anchors", {
  # L-shaped polygon whose centroid falls outside the shape
  ell <- st_polygon(list(rbind(c(0, 0), c(4, 0), c(4, 1), c(1, 1), c(1, 4),
                               c(0, 4), c(0, 0))))
  sq <- st_polygon(list(rbind(c(6, 6), c(8, 6), c(8, 8), c(6, 8), c(6, 6))))
  demand <- st_sf(pop = c(10, 20), geometry = st_sfc(ell, sq))
  facilities <- st_as_sf(data.frame(x = c(2, 7), y = c(2, 7)),
                         coords = c("x", "y"))
  res <- p_median(demand, facilities, n_facilities = 2, weight_col = "pop")

  cen <- spider_lines(res)
  expect_equal(unname(line_ends(cen)$start),
               unname(st_coordinates(st_centroid(st_geometry(demand)))))

  pos <- spider_lines(res, anchor = "point_on_surface")
  starts <- st_as_sf(as.data.frame(line_ends(pos)$start), coords = c(1, 2))
  expect_true(all(lengths(st_within(starts, demand)) == 1))
})

test_that("coincident endpoints give zero-length lines", {
  pts <- st_as_sf(data.frame(x = c(0, 5), y = c(0, 5), pop = c(1, 1)),
                  coords = c("x", "y"))
  res <- p_median(pts, pts, n_facilities = 2, weight_col = "pop")
  lines <- spider_lines(res)
  expect_equal(nrow(lines), 2)
  expect_equal(as.numeric(lines$length), c(0, 0))
})

test_that("CRS rules: transform, missing on one side, both missing", {
  d <- make_spider_data(crs = 32613)
  d$demand <- st_transform(d$demand, 4326)
  cm <- distance_matrix(st_transform(d$demand, 32613), d$facilities)
  res <- p_median(d$demand, d$facilities, n_facilities = 3,
                  weight_col = "pop", cost_matrix = cm)
  lines <- spider_lines(res)
  expect_equal(st_crs(lines), st_crs(4326))
  expect_s3_class(lines$length, "units")

  res_bad <- res
  st_crs(res_bad$facilities) <- NA
  expect_error(spider_lines(res_bad), "both have a CRS")

  d0 <- make_spider_data()
  res0 <- p_median(d0$demand, d0$facilities, n_facilities = 3,
                   weight_col = "pop")
  lines0 <- spider_lines(res0)
  expect_true(is.na(st_crs(lines0)))
  expect_type(lines0$length, "double")
  expect_error(spider_lines(res0, route_fun = function(from, to) NULL),
               "known CRS")
})

test_that("unsupported results error clearly", {
  expect_error(spider_lines(list(a = 1)), "must come from")
})

# ---------------------------------------------------------------------------
# route_fun
# ---------------------------------------------------------------------------

make_routable <- function() {
  d <- make_spider_data(crs = 32613)
  # Duplicate demand locations so some routes are requested twice
  d$demand <- rbind(d$demand, d$demand[1:5, ])
  p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
}

test_that("route_fun is called once on unique pairs, with fallbacks", {
  res <- make_routable()
  calls <- 0L
  n_seen <- NA_integer_
  dogleg <- function(from, to) {
    calls <<- calls + 1L
    n_seen <<- nrow(from)
    expect_equal(st_crs(from), st_crs(4326))
    expect_true("id" %in% names(from))
    a <- st_coordinates(from)
    b <- st_coordinates(to)
    geoms <- lapply(seq_len(nrow(a)), function(i) {
      if (i == 1) return(st_linestring())  # first pair fails
      st_linestring(rbind(a[i, ], c(a[i, 1], b[i, 2]), b[i, ]))
    })
    st_sfc(geoms, crs = 4326)
  }

  expect_warning(lines <- spider_lines(res, route_fun = dogleg),
                 "route\\(s\\) failed")
  expect_equal(calls, 1L)
  meta <- attr(lines, "spopt")
  expect_equal(meta$n_route_requests, n_seen)
  expect_lt(n_seen, nrow(lines))
  expect_equal(meta$n_route_failures, 1)
  expect_true(all(st_geometry_type(lines) == "MULTILINESTRING"))
  expect_true(any(!lines$routed))
  expect_true(any(lines$routed))
  # Routed lines have 3 vertices; fallbacks have 2
  nv <- table(st_coordinates(lines)[, "L2"])
  expect_true(all(nv[lines$routed] == 3))
  expect_true(all(nv[!lines$routed] == 2))
})

test_that("malformed route_fun returns error", {
  res <- make_routable()
  short <- function(from, to) st_sfc(st_linestring(), crs = 4326)
  expect_error(spider_lines(res, route_fun = short), "returned 1 geometries")

  no_crs <- function(from, to) {
    st_sfc(lapply(seq_len(nrow(from)), function(i) st_linestring()))
  }
  expect_error(spider_lines(res, route_fun = no_crs), "known CRS")

  pts <- function(from, to) st_geometry(st_transform(from, 4326))
  expect_error(spider_lines(res, route_fun = pts), "LINESTRING")

  expect_error(spider_lines(res, route_fun = function(from, to) "x"),
               "sf or sfc")
})

# ---------------------------------------------------------------------------
# Review follow-ups
# ---------------------------------------------------------------------------

test_that("cflp .cost is the cost to the primary facility", {
  res <- make_split_data()
  cm <- distance_matrix(res$demand, res$facilities)
  i <- seq_len(nrow(res$demand))
  expect_equal(res$demand$.cost, cm[cbind(i, res$demand$.facility)])
})

test_that(".cost keeps an original Inf rather than the penalty", {
  d <- make_spider_data(n_demand = 6, n_fac = 3)
  cm <- distance_matrix(d$demand, d$facilities)
  cm[1, ] <- Inf
  pm <- suppressWarnings(
    p_median(d$demand, d$facilities, n_facilities = 2, weight_col = "pop",
             cost_matrix = cm)
  )
  expect_identical(pm$demand$.cost[1], Inf)
})

test_that("explicit cost_matrix overrides .cost", {
  d <- make_spider_data()
  res <- p_median(d$demand, d$facilities, n_facilities = 3, weight_col = "pop")
  other <- distance_matrix(d$demand, d$facilities) * 10 + 1
  lines <- spider_lines(res, cost_matrix = other)
  expect_equal(lines$cost,
               other[cbind(lines$demand_index, lines$facility_index)])
  expect_false(isTRUE(all.equal(lines$cost, res$demand$.cost)))
})

test_that("MULTIPOINT and Z geometries are accepted", {
  mp <- st_sfc(st_multipoint(rbind(c(0, 0), c(2, 0))),
               st_multipoint(rbind(c(5, 5), c(7, 7))))
  demand <- st_sf(pop = c(1, 1), geometry = mp)
  facilities <- st_as_sf(data.frame(x = c(1, 6), y = c(0, 6), z = c(9, 9)),
                         coords = c("x", "y", "z"))
  res <- p_median(demand, facilities, n_facilities = 2, weight_col = "pop")
  lines <- spider_lines(res)
  expect_equal(nrow(lines), 2)
  expect_equal(unname(line_ends(lines)$start),
               unname(st_coordinates(st_centroid(mp))))
  expect_equal(class(st_geometry(lines)[[1]])[1], "XY")
})

test_that("invalid polygons are dropped and counted", {
  bowtie <- st_polygon(list(rbind(c(0, 0), c(2, 2), c(2, 0), c(0, 2),
                                  c(0, 0))))
  sq <- st_polygon(list(rbind(c(5, 5), c(6, 5), c(6, 6), c(5, 6), c(5, 5))))
  demand <- st_sf(pop = c(1, 1), geometry = st_sfc(bowtie, sq))
  facilities <- st_as_sf(data.frame(x = c(1, 5), y = c(1, 5)),
                         coords = c("x", "y"))
  cm <- matrix(c(1, 9, 9, 1), 2)
  res <- p_median(demand, facilities, n_facilities = 2, weight_col = "pop",
                  cost_matrix = cm)
  expect_warning(lines <- spider_lines(res), "invalid")
  expect_equal(nrow(lines), 1)
  expect_equal(attr(lines, "spopt")$n_invalid_geometry, 1)
})

test_that("'all' links outside a filtered table warn instead of failing", {
  d <- make_spider_data(n_demand = 10, n_fac = 4)
  res <- huff(d$demand, d$facilities, attractiveness_col = "sqft")
  res$stores <- res$stores[1:3, ]
  msgs <- character()
  lines <- withCallingHandlers(
    suppressMessages(spider_lines(res, allocations = "all")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("allocation\\(s\\) point outside", msgs)))
  expect_true(all(lines$facility_index <= 3))
})

test_that("cflp primary links below the residue tolerance are dropped", {
  res <- make_split_data()
  alloc <- attr(res, "spopt")$allocation_matrix
  alloc[2, res$demand$.facility[2]] <- 1e-9
  attr(res, "spopt")$allocation_matrix <- alloc
  lines <- spider_lines(res)
  expect_false(2 %in% lines$demand_index)
})

test_that("non-finite endpoint coordinates are dropped and counted", {
  demand <- st_sf(pop = c(1, 1, 1),
                  geometry = st_sfc(st_point(c(Inf, 0)), st_point(c(1, 1)),
                                    st_point(c(2, 2))))
  facilities <- st_sf(geometry = st_sfc(st_point(c(0, 0)),
                                        st_point(c(-Inf, 5))))
  cm <- matrix(c(1, 1, 1, 9, 9, 9), 3)
  cm[3, ] <- c(9, 1)  # demand 3 goes to the non-finite facility
  res <- p_median(demand, facilities, n_facilities = 2, weight_col = "pop",
                  cost_matrix = cm)
  expect_equal(res$demand$.facility, c(1L, 1L, 2L))
  expect_warning(lines <- spider_lines(res), "invalid")
  expect_equal(lines$demand_index, 2L)
  expect_true(all(is.finite(as.numeric(lines$length))))
  expect_equal(attr(lines, "spopt")$n_invalid_geometry, 2)
})

# ---------------------------------------------------------------------------
# r5r itinerary assembly (no Java needed)
# ---------------------------------------------------------------------------

test_that(".assemble_itineraries orders legs, keeps option 1, fills gaps", {
  leg <- function(...) st_linestring(rbind(...))
  it <- st_sf(
    from_id = c("1", "1", "1", "3"),
    to_id = c("1", "1", "1", "3"),
    option = c(1L, 1L, 2L, 1L),
    segment = c(2L, 1L, 1L, 1L),
    geometry = st_sfc(
      leg(c(1, 1), c(2, 2)),          # pair 1, option 1, segment 2
      leg(c(0, 0), c(1, 1)),          # pair 1, option 1, segment 1
      leg(c(0, 0), c(5, 5)),          # pair 1, option 2 (ignored)
      leg(c(3, 3), c(4, 4), c(5, 3)), # pair 3
      crs = 4326
    )
  )
  out <- spopt:::.assemble_itineraries(it, ids = c("1", "2", "3"))

  expect_s3_class(out, "sfc")
  expect_equal(length(out), 3)
  expect_equal(st_crs(out), st_crs(4326))
  expect_true(st_is_empty(out[2]))  # pair 2 had no itinerary

  p1 <- unclass(out[[1]])
  expect_length(p1, 2)
  expect_equal(p1[[1]][1, ], c(0, 0))  # segment 1 first
  expect_equal(p1[[2]][2, ], c(2, 2))
  expect_equal(nrow(unclass(out[[3]])[[1]]), 3)

  none <- spopt:::.assemble_itineraries(it[0, ], ids = c("1", "2"))
  expect_true(all(st_is_empty(none)))
})
