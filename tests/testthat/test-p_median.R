test_that("p_median returns correct structure", {
  skip_if_not_installed("sf")
  skip("Rust compilation required")

  set.seed(42)
  demand <- sf::st_as_sf(
    data.frame(x = runif(30), y = runif(30), pop = rpois(30, 100)),
    coords = c("x", "y")
  )
  facilities <- sf::st_as_sf(
    data.frame(x = runif(10), y = runif(10)),
    coords = c("x", "y")
  )

  result <- p_median(demand, facilities, n_facilities = 3, weight_col = "pop")

  expect_type(result, "list")
  expect_s3_class(result$demand, "sf")
  expect_s3_class(result$facilities, "sf")
  expect_true(".facility" %in% names(result$demand))
  expect_true(".selected" %in% names(result$facilities))
  expect_equal(sum(result$facilities$.selected), 3)
})

test_that("p_median assigns all demand to selected facilities", {
  skip_if_not_installed("sf")
  skip("Rust compilation required")

  set.seed(42)
  demand <- sf::st_as_sf(
    data.frame(x = runif(20), y = runif(20), pop = rep(1, 20)),
    coords = c("x", "y")
  )
  facilities <- sf::st_as_sf(
    data.frame(x = runif(8), y = runif(8)),
    coords = c("x", "y")
  )

  result <- p_median(demand, facilities, n_facilities = 4, weight_col = "pop")

  # All assignments should be to selected facilities
  selected_ids <- which(result$facilities$.selected)
  expect_true(all(result$demand$.facility %in% selected_ids))
})
