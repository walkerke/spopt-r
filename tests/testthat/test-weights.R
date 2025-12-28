test_that("sp_weights creates queen contiguity", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")

  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

  w <- sp_weights(nc, type = "queen")

  expect_s3_class(w, "nb")
  expect_equal(length(w), nrow(nc))
})

test_that("sp_weights creates rook contiguity", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")

  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

  w <- sp_weights(nc, type = "rook")

  expect_s3_class(w, "nb")
  # Rook should have <= neighbors compared to queen
})

test_that("sp_weights errors on non-sf input", {
  expect_error(sp_weights(data.frame(x = 1:5)), "sf object")
})
