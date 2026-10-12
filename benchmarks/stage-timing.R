# Stage timing of the regionalization workflows on synthetic grids.
# Usage: Rscript benchmarks/stage-timing.R [libpath] [label] [sizes]
# Splits each public call into: weights (poly2nb + connectivity), Rust kernel,
# and the rest (validation, scaling, result assembly). Writes a CSV next to it.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
label <- if (length(args) >= 2) args[2] else "current"
sizes <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(500L, 2500L, 10000L)
suppressPackageStartupMessages({library(sf); library(spopt)})
cat(label, "spopt from", find.package("spopt"), "\n")
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "fixtures.R"))
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), sprintf("stage-timing-%s.csv", label))
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- force(expr); c(secs = proc.time()[["elapsed"]] - t0) }
rows <- list()
rec <- function(n, stage, secs, note = "") { rows[[length(rows) + 1]] <<- data.frame(label, n, stage, secs = round(secs, 3), note); cat(sprintf("%-8s n=%6d  %-34s %8.3fs %s\n", label, n, stage, secs, note)) }
for (n in sizes) {
  fx <- make_grid_fixture(n, p = 4, seed = 1); attrs <- paste0("v", 1:4); n_act <- nrow(fx)
  rec(n_act, "weights: poly2nb queen", tm(nb <- spdep::poly2nb(fx, queen = TRUE)))
  rec(n_act, "weights: prepare_weights()", tm(spopt:::prepare_weights(fx, "queen")))
  rec(n_act, "skater, nb given", tm(sk <- skater(fx, attrs, n_regions = 10, weights = nb)))
  rec(n_act, "skater, floor, nb given", tm(skater(fx, attrs, n_regions = 10, weights = nb, floor = "pop", floor_value = 20000)))
  rec(n_act, "ward_spatial, nb given", tm(ward_spatial(fx, attrs, n_regions = 10, weights = nb)))
  if ("initial_regions" %in% names(formals(azp))) {
    rec(n_act, "azp basic from skater, nb given", tm(azp(fx, attrs, weights = nb, method = "basic", initial_regions = sk, seed = 1)))
    rec(n_act, "azp sa from skater, nb given", tm(azp(fx, attrs, weights = nb, method = "sa", initial_regions = sk, seed = 1)))
  }
  if ("redcap" %in% getNamespaceExports("spopt")) rec(n_act, "redcap fullorder-ward, nb given", tm(redcap(fx, attrs, n_regions = 10, weights = nb)))
  rec(n_act, "max_p 20k floor, nb given", tm(mp <- max_p_regions(fx, attrs, threshold_var = "pop", threshold = 20000, weights = nb, seed = 1)), sprintf("p=%d", attr(mp, "spopt")$n_regions))
  rec(n_act, "max_p construction only", tm(max_p_regions(fx, attrs, threshold_var = "pop", threshold = 20000, weights = nb, seed = 1, n_sa_iterations = 0)))
  if (n_act <= 2500) rec(n_act, "spenc, nb given", tm(spenc(fx, attrs, n_regions = 10, weights = nb, seed = 1)))
  if (n_act <= 2500) { fi <- make_grid_fixture(n, p = 4, seed = 1, islands = 3)
    rec(nrow(fi), "weights: prepare_weights(bridge)", tm(suppressMessages(spopt:::prepare_weights(fi, "queen", bridge_islands = TRUE))), "3 islands") }
}
res <- do.call(rbind, rows); write.csv(res, out, row.names = FALSE); cat("wrote", out, "\n")
