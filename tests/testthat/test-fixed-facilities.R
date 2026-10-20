# Shared test data
make_test_data <- function(fixed_logical = NULL, fixed_status = NULL) {
  set.seed(42)
  demand <- sf::st_as_sf(
    data.frame(x = runif(20), y = runif(20), pop = rpois(20, 100)),
    coords = c("x", "y")
  )
  fac_df <- data.frame(x = runif(8), y = runif(8))
  if (!is.null(fixed_logical)) fac_df$fixed <- fixed_logical
  if (!is.null(fixed_status)) fac_df$status <- fixed_status
  facilities <- sf::st_as_sf(fac_df, coords = c("x", "y"))
  list(demand = demand, facilities = facilities)
}

test_that("mclp with fixed_col (logical) includes fixed in selected", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, TRUE, rep(FALSE, 6)))
  result <- mclp(d$demand, d$facilities, service_radius = 0.5,
                  n_facilities = 3, weight_col = "pop",
                  fixed_col = "fixed")

  meta <- attr(result, "spopt")
  selected <- which(result$facilities$.selected)

  expect_true(1 %in% selected)
  expect_true(2 %in% selected)
  expect_equal(meta$n_selected, 3)
  expect_equal(meta$n_fixed, 2)
  expect_true(result$facilities$.fixed[1])
  expect_true(result$facilities$.fixed[2])
  expect_false(result$facilities$.fixed[3])
})

test_that("mclp with fixed_col (required/candidate) works", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_status = c("required", "required",
                                        rep("candidate", 6)))
  result <- mclp(d$demand, d$facilities, service_radius = 0.5,
                  n_facilities = 3, weight_col = "pop",
                  fixed_col = "status")

  selected <- which(result$facilities$.selected)
  expect_true(1 %in% selected)
  expect_true(2 %in% selected)
  expect_equal(length(selected), 3)
})

test_that("p_median with fixed_col includes fixed in selected", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, rep(FALSE, 7)))
  result <- p_median(d$demand, d$facilities, n_facilities = 3,
                      weight_col = "pop", fixed_col = "fixed")

  meta <- attr(result, "spopt")
  selected <- which(result$facilities$.selected)

  expect_true(1 %in% selected)
  expect_equal(length(selected), 3)
  expect_equal(meta$n_fixed, 1)
})

test_that("p_center binary_search with fixed_col includes fixed in selected", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, rep(FALSE, 7)))
  result <- p_center(d$demand, d$facilities, n_facilities = 3,
                      method = "binary_search", fixed_col = "fixed")

  selected <- which(result$facilities$.selected)
  expect_true(1 %in% selected)
  expect_equal(length(selected), 3)
  expect_true(result$facilities$.fixed[1])
})

test_that("p_center mip with fixed_col includes fixed in selected", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, rep(FALSE, 7)))
  result <- p_center(d$demand, d$facilities, n_facilities = 3,
                      method = "mip", fixed_col = "fixed")

  selected <- which(result$facilities$.selected)
  expect_true(1 %in% selected)
  expect_equal(length(selected), 3)
  expect_true(result$facilities$.fixed[1])
})

test_that("all facilities fixed (degenerate case) works", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, TRUE, TRUE, rep(FALSE, 5)))
  result <- p_median(d$demand, d$facilities, n_facilities = 3,
                      weight_col = "pop", fixed_col = "fixed")

  meta <- attr(result, "spopt")
  selected <- which(result$facilities$.selected)

  expect_equal(sort(selected), c(1, 2, 3))
  expect_equal(meta$n_fixed, 3)
  expect_true(all(result$facilities$.fixed[1:3]))
})

test_that("fixed_col validation catches errors", {
  skip_if_not_installed("sf")

  d <- make_test_data()

  # Column not found
  expect_error(
    p_median(d$demand, d$facilities, n_facilities = 3,
             weight_col = "pop", fixed_col = "nonexistent"),
    "not found"
  )

  # NA values in column
  d2 <- make_test_data(fixed_logical = c(TRUE, NA, rep(FALSE, 6)))
  expect_error(
    p_median(d2$demand, d2$facilities, n_facilities = 3,
             weight_col = "pop", fixed_col = "fixed"),
    "NA"
  )

  # More fixed than requested
  d3 <- make_test_data(fixed_logical = c(TRUE, TRUE, TRUE, rep(FALSE, 5)))
  expect_error(
    p_median(d3$demand, d3$facilities, n_facilities = 2,
             weight_col = "pop", fixed_col = "fixed"),
    "more fixed"
  )

  # Bad character values
  d4 <- make_test_data(fixed_status = c("required", "open", rep("candidate", 6)))
  expect_error(
    p_median(d4$demand, d4$facilities, n_facilities = 3,
             weight_col = "pop", fixed_col = "status"),
    "required.*candidate"
  )
})

test_that("NULL fixed_col preserves backward compatibility", {
  skip_if_not_installed("sf")

  d <- make_test_data()
  result <- p_median(d$demand, d$facilities, n_facilities = 3,
                      weight_col = "pop", fixed_col = NULL)

  meta <- attr(result, "spopt")
  expect_equal(meta$n_fixed, 0)
  expect_true(all(!result$facilities$.fixed))
})

test_that("cflp with fixed_col opens fixed facilities", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(FALSE, FALSE, FALSE, TRUE, rep(FALSE, 4)))
  d$facilities$cap <- 1000
  free <- cflp(d$demand, d$facilities, n_facilities = 3,
               weight_col = "pop", capacity_col = "cap")
  result <- cflp(d$demand, d$facilities, n_facilities = 3,
                 weight_col = "pop", capacity_col = "cap", fixed_col = "fixed")

  meta <- attr(result, "spopt")
  selected <- which(result$facilities$.selected)

  expect_true(4 %in% selected)
  expect_equal(length(selected), 3)
  expect_equal(meta$n_fixed, 1)
  expect_equal(meta$fixed_facilities, 4L)
  expect_equal(result$facilities$.fixed, c(FALSE, FALSE, FALSE, TRUE, rep(FALSE, 4)))
  # Forcing a site open can't beat the unconstrained optimum
  expect_gte(meta$objective, attr(free, "spopt")$objective - 1e-6)
  # Capacities still hold
  expect_true(all(result$facilities$.utilization <= 1 + 1e-9))
})

test_that("cflp fixed_col works in cost-based mode", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_status = c("required", "required", "required",
                                       rep("candidate", 5)))
  d$facilities$cap <- 2000
  d$facilities$fixed_cost <- 1e6  # high enough that the solver opens one site
  result <- cflp(d$demand, d$facilities, n_facilities = 0,
                 weight_col = "pop", capacity_col = "cap",
                 facility_cost_col = "fixed_cost", fixed_col = "status")

  selected <- which(result$facilities$.selected)
  expect_true(all(1:3 %in% selected))
  expect_equal(attr(result, "spopt")$n_fixed, 3)
})

test_that("cflp fixed_col validation and NULL default", {
  skip_if_not_installed("sf")

  d <- make_test_data(fixed_logical = c(TRUE, TRUE, TRUE, rep(FALSE, 5)))
  d$facilities$cap <- 1000
  expect_error(
    cflp(d$demand, d$facilities, n_facilities = 2,
         weight_col = "pop", capacity_col = "cap", fixed_col = "fixed"),
    "more fixed"
  )

  result <- cflp(d$demand, d$facilities, n_facilities = 3,
                 weight_col = "pop", capacity_col = "cap")
  expect_equal(attr(result, "spopt")$n_fixed, 0)
  expect_true(all(!result$facilities$.fixed))
})

test_that("cflp keeps positional compatibility for verbose", {
  skip_if_not_installed("sf")

  d <- make_test_data()
  d$facilities$cap <- 1000
  # Calls written before fixed_col existed pass verbose in eighth position
  result <- cflp(d$demand, d$facilities, 3, "pop", "cap", NULL, NULL,
                 "euclidean", FALSE)
  expect_equal(sum(result$facilities$.selected), 3)
  expect_equal(attr(result, "spopt")$n_fixed, 0)
})
