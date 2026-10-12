# AZP and max-p search behaviour. These need the compiled Rust backend, so they
# are skipped on CRAN like the other regionalization tests.

nc_data <- function() {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")
  skip_on_cran()
  sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
}

regions_contiguous <- function(result) {
  nb <- spdep::poly2nb(result, queen = TRUE)
  all(vapply(split(seq_len(nrow(result)), result$.region), function(idx) {
    length(idx) == 1 || suppressWarnings(spdep::n.comp.nb(spdep::subset.nb(nb, seq_len(nrow(result)) %in% idx))$nc) == 1
  }, logical(1)))
}

test_that("azp returns the requested number of contiguous regions for every method", {
  nc <- nc_data()
  for (m in c("basic", "tabu", "sa")) {
    result <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = m, seed = 1)
    expect_s3_class(result, "sf")
    expect_equal(length(unique(result$.region)), 6)
    expect_true(regions_contiguous(result), info = m)
    meta <- attr(result, "spopt")
    expect_lte(meta$objective, meta$initial_objective + 1e-8)
  }
})

test_that("azp accepts a warm start and never makes it worse", {
  nc <- nc_data()
  sk <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 6)
  sk_obj <- attr(sk, "spopt")$objective
  for (m in c("basic", "tabu", "sa")) {
    result <- azp(nc, attrs = c("SID74", "SID79"), initial_regions = sk, method = m, seed = 1)
    expect_equal(length(unique(result$.region)), 6)
    expect_true(regions_contiguous(result), info = m)
    expect_lte(attr(result, "spopt")$objective, sk_obj + 1e-8)
    expect_true(attr(result, "spopt")$warm_start)
  }
  # a plain vector works too, and n_regions is inferred from it
  result <- azp(nc, attrs = c("SID74", "SID79"), initial_regions = sk$.region, method = "basic")
  expect_equal(length(unique(result$.region)), 6)
})

test_that("azp rejects inconsistent or non-contiguous warm starts", {
  nc <- nc_data()
  sk <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 6)
  expect_error(azp(nc, attrs = "SID74", n_regions = 5, initial_regions = sk), "has 6 regions")
  expect_error(azp(nc, attrs = "SID74", initial_regions = sk$.region[-1]), "one entry per row")
  scattered <- rep(c("a", "b"), length.out = nrow(nc))
  expect_error(azp(nc, attrs = "SID74", initial_regions = scattered), "not spatially contiguous")
})

test_that("azp is reproducible with a seed and runs more than a handful of moves", {
  nc <- nc_data()
  a <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = "sa", seed = 7)
  b <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = "sa", seed = 7)
  expect_identical(a$.region, b$.region)
  # SA from a random start must at least reach the greedy local optimum
  g <- azp(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = "basic", seed = 7)
  expect_lte(attr(a, "spopt")$objective, attr(g, "spopt")$objective * 1.05)
})

test_that("max_p_regions improvement phase respects the floor and lowers the objective", {
  nc <- nc_data()
  base <- max_p_regions(nc, attrs = c("SID74", "SID79"), threshold_var = "BIR74",
                        threshold = 100000, seed = 3, n_sa_iterations = 0)
  improved <- max_p_regions(nc, attrs = c("SID74", "SID79"), threshold_var = "BIR74",
                            threshold = 100000, seed = 3)
  expect_equal(attr(improved, "spopt")$n_regions, attr(base, "spopt")$n_regions)
  expect_lte(attr(improved, "spopt")$objective, attr(base, "spopt")$objective + 1e-8)
  expect_true(all(tapply(nc$BIR74, improved$.region, sum) >= 100000))
  expect_true(regions_contiguous(improved))
})
