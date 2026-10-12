# Synthetic regionalization fixtures: a hexagonal grid of n polygons with p
# spatially smooth attributes (sums of Gaussian bumps plus noise) and a
# population column. Deterministic given the seed. `islands` cuts the grid into
# that many disconnected pieces by deleting cells along vertical bands.
make_grid_fixture <- function(n = 2500, p = 4, seed = 1, islands = 0) {
  set.seed(seed)
  side <- ceiling(sqrt(n))
  g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = side, ymax = side))),
                        cellsize = 1, square = FALSE)
  g <- g[seq_len(min(n, length(g)))]
  if (islands > 0) {
    cx <- sf::st_coordinates(sf::st_centroid(g))[, 1]
    cuts <- quantile(cx, probs = seq_len(islands) / (islands + 1))
    keep <- !Reduce(`|`, lapply(cuts, function(c) abs(cx - c) < 0.6))
    g <- g[keep]
  }
  m <- length(g); xy <- sf::st_coordinates(sf::st_centroid(g))
  bump <- function() { f <- numeric(m); for (b in 1:6) { cx <- runif(1, 0, side); cy <- runif(1, 0, side); s <- runif(1, side / 10, side / 3)
    f <- f + rnorm(1) * exp(-((xy[, 1] - cx)^2 + (xy[, 2] - cy)^2) / (2 * s^2)) }; as.numeric(scale(f)) * 0.9 + rnorm(m) * 0.45 }
  df <- as.data.frame(replicate(p, bump())); names(df) <- paste0("v", seq_len(p))
  df$pop <- round(runif(m, 1000, 8000))
  sf::st_sf(df, geometry = g)
}
