# bridge_islands = TRUE must always return a connected neighbour graph.
# Pure R (spdep), so no Rust build is needed.

three_islands <- function() {
  # A mainland (x < 10) and two small islands (12 < x < 14, 15 < x < 17) that
  # are closer to each other (gap 1) than to the mainland (gap 2). Linking
  # every island only to its nearest unit in another component joins the two
  # islands to each other and never to the mainland.
  g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 17, ymax = 4))),
                        cellsize = 1, square = TRUE)
  cx <- sf::st_coordinates(sf::st_centroid(g))[, 1]
  keep <- cx < 10 | (cx > 12 & cx < 14) | (cx > 15 & cx < 17)
  sf::st_sf(v = seq_along(g[keep]), geometry = g[keep])
}

test_that("bridging connects every component, not just each island to its nearest island", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")
  x <- three_islands()
  nb0 <- suppressWarnings(spdep::poly2nb(x, queen = TRUE))
  expect_gt(spdep::n.comp.nb(nb0)$nc, 1)
  expect_message(prepare_weights(x, "queen", bridge_islands = TRUE), "Bridged")
  nb <- suppressMessages(prepare_weights(x, "queen", bridge_islands = TRUE))
  expect_s3_class(nb, "nb")
  expect_equal(spdep::n.comp.nb(nb)$nc, 1)
  # symmetric and no self-neighbours
  edges <- data.frame(i = rep(seq_along(nb), lengths(nb)), j = unlist(nb))
  expect_false(any(edges$i == edges$j))
  expect_true(all(paste(edges$j, edges$i) %in% paste(edges$i, edges$j)))
})

test_that("a graph with an isolate is also bridged to one component", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")
  x <- three_islands()
  # Detach one square completely by moving it away
  g <- sf::st_geometry(x); g[1] <- g[1] + c(0, 40); sf::st_geometry(x) <- g
  nb <- suppressMessages(prepare_weights(x, "queen", bridge_islands = TRUE))
  expect_equal(spdep::n.comp.nb(nb)$nc, 1)
  expect_true(all(spdep::card(nb) > 0))
})
