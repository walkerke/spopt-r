test_that("lscp returns correct structure", {
  skip_if_not_installed("sf")
  skip("Rust compilation required")

  # Create simple test data
  set.seed(42)
  demand <- sf::st_as_sf(
    data.frame(x = runif(20), y = runif(20)),
    coords = c("x", "y")
  )
  facilities <- sf::st_as_sf(
    data.frame(x = runif(8), y = runif(8)),
    coords = c("x", "y")
  )

  result <- lscp(demand, facilities, service_radius = 0.5)

  expect_type(result, "list")
  expect_s3_class(result$demand, "sf")
  expect_s3_class(result$facilities, "sf")
  expect_true(".covered" %in% names(result$demand))
  expect_true(".selected" %in% names(result$facilities))
})

test_that("lscp covers all demand when feasible", {
  skip_if_not_installed("sf")
  skip("Rust compilation required")

  # Dense facilities should cover everything
  set.seed(42)
  demand <- sf::st_as_sf(
    data.frame(x = runif(10), y = runif(10)),
    coords = c("x", "y")
  )
  facilities <- sf::st_as_sf(
    data.frame(
      x = rep(seq(0, 1, 0.2), each = 6),
      y = rep(seq(0, 1, 0.2), times = 6)
    ),
    coords = c("x", "y")
  )

  result <- lscp(demand, facilities, service_radius = 0.3)

  expect_equal(attr(result, "spopt")$coverage_pct, 100)
})
