skip_if_not_installed("sf")
library(sf)

test_that("MULTIPOINT geometries use one centroid per feature", {
  mp <- st_sf(geometry = st_sfc(
    st_multipoint(rbind(c(0, 0), c(2, 0))),
    st_multipoint(rbind(c(4, 4), c(6, 4), c(5, 7)))
  ))
  pts <- st_as_sf(data.frame(x = c(1, 5), y = c(0, 5)), coords = c("x", "y"))

  d <- distance_matrix(mp, pts)
  expect_equal(dim(d), c(2L, 2L))
  cen <- st_coordinates(st_centroid(st_geometry(mp)))
  expected <- sqrt(outer(cen[, 1], c(1, 5), "-")^2 +
                     outer(cen[, 2], c(0, 5), "-")^2)
  expect_equal(unname(d), unname(expected))

  # Reversed roles and manhattan distances
  expect_equal(dim(distance_matrix(pts, mp, type = "manhattan")), c(2L, 2L))
})

test_that("MULTIPOINT works in geographic CRS", {
  mp <- st_sf(geometry = st_sfc(
    st_multipoint(rbind(c(-97.3, 32.7), c(-97.1, 32.7))),
    st_multipoint(rbind(c(-96.8, 32.8), c(-96.7, 32.9))),
    crs = 4326
  ))
  pts <- st_as_sf(data.frame(x = c(-97.2, -96.75), y = c(32.7, 32.85)),
                  coords = c("x", "y"), crs = 4326)
  d <- distance_matrix(mp, pts)
  expect_equal(dim(d), c(2L, 2L))
  expect_equal(
    as.numeric(d),
    as.numeric(st_distance(st_centroid(st_geometry(mp)), pts)),
    tolerance = 1e-6
  )
})

test_that("solvers accept MULTIPOINT demand", {
  mp <- st_sf(pop = c(10, 20, 30), geometry = st_sfc(
    st_multipoint(rbind(c(0, 0), c(1, 0))),
    st_multipoint(rbind(c(5, 5), c(6, 5))),
    st_multipoint(rbind(c(9, 0), c(9, 1)))
  ))
  fac <- st_as_sf(data.frame(x = c(0, 5, 9), y = c(0, 5, 0)),
                  coords = c("x", "y"))
  res <- p_median(mp, fac, n_facilities = 3, weight_col = "pop")
  expect_equal(res$demand$.facility, 1:3)
})
