test_that("skater returns sf with .region column",
{
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")

  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

  result <- skater(nc, attrs = c("SID74", "SID79"), n_regions = 5)

  expect_s3_class(result, "sf")
  expect_true(".region" %in% names(result))
  expect_equal(length(unique(result$.region)), 5)
})

test_that("skater respects floor constraint", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")

  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

  result <- skater(
    nc,
    attrs = c("SID74", "SID79"),
    n_regions = 5,
    floor = "BIR74",
    floor_value = 40000
  )

  # Check each region meets floor
  region_births <- tapply(nc$BIR74, result$.region, sum)
  expect_equal(length(region_births), 5)
  expect_true(all(region_births >= 40000))
})
