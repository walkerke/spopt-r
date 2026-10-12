# REDCAP, all five linkages. Needs the compiled Rust backend, so skipped on
# CRAN like the other regionalization tests.

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

methods <- c("fullorder-averagelinkage", "fullorder-singlelinkage",
             "fullorder-completelinkage", "fullorder-wardlinkage",
             "firstorder-singlelinkage")

test_that("every linkage returns the requested number of contiguous regions", {
  nc <- nc_data()
  for (m in methods) {
    result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = m)
    expect_s3_class(result, "sf")
    expect_equal(length(unique(result$.region)), 6, info = m)
    expect_true(regions_contiguous(result), info = m)
    meta <- attr(result, "spopt")
    expect_equal(meta$algorithm, "redcap")
    expect_equal(meta$method, m)
    expect_gt(meta$objective, 0)
  }
})

test_that("the default method is full-order average linkage, as in GeoDa and rgeoda", {
  nc <- nc_data()
  a <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6)
  b <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = "fullorder-averagelinkage")
  expect_identical(a$.region, b$.region)
  expect_equal(attr(a, "spopt")$method, "fullorder-averagelinkage")
  expect_error(redcap(nc, attrs = "SID74", n_regions = 4, method = "single"), "should be one of")
})

test_that("every linkage respects a floor constraint", {
  nc <- nc_data()
  for (m in methods) {
    result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 4, method = m,
                     floor = "BIR74", floor_value = 50000)
    expect_true(all(tapply(nc$BIR74, result$.region, sum) >= 50000), info = m)
    expect_true(regions_contiguous(result), info = m)
  }
})

test_that("first-order single linkage is SKATER's tree", {
  nc <- nc_data()
  r <- redcap(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = 8, method = "firstorder-singlelinkage")
  s <- skater(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = 8)
  expect_equal(attr(r, "spopt")$objective, attr(s, "spopt")$objective)
  tb <- table(r$.region, s$.region)
  expect_true(all(rowSums(tb > 0) == 1) && all(colSums(tb > 0) == 1))
})

test_that("Ward linkage stays close to the constrained Ward partition at the same k", {
  # Same merge sequence as ward_spatial(); each cut re-optimizes the split, so
  # the two usually land within a few percent (REDCAP is typically, not
  # always, the more homogeneous). Guards against a broken tree only.
  nc <- nc_data()
  for (k in c(4, 6, 10)) {
    r <- redcap(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = k, method = "fullorder-wardlinkage")
    w <- ward_spatial(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = k)
    expect_lte(attr(r, "spopt")$objective, attr(w, "spopt")$objective * 1.10)
  }
})
