# Capture or compare complete regionalization outputs across package builds.
# Usage: Rscript benchmarks/equivalence.R capture <libpath> <dir>
#        Rscript benchmarks/equivalence.R compare <libpath> <dir>
# Partitions are compared up to relabelling (every region in A maps to exactly
# one region in B), objectives to 1e-9 relative, and metadata fields by value.
args <- commandArgs(trailingOnly = TRUE); mode <- args[1]; .libPaths(c(args[2], .libPaths())); dir <- args[3]
suppressPackageStartupMessages({library(sf); library(spdep); library(spopt)})
cat(mode, "with spopt from", find.package("spopt"), "\n")
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))); source(file.path(here, "fixtures.R"))
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
# Real-data fixtures live next to this script (see fetch-tracts.R to regenerate)
meck <- st_read(file.path(here, "data", "meck.gpkg"), quiet = TRUE)
van <- st_read(file.path(here, "data", "vancouver.gpkg"), quiet = TRUE); cmp <- n.comp.nb(poly2nb(van)); van <- van[cmp$comp.id == which.max(tabulate(cmp$comp.id)), ]
g500 <- make_grid_fixture(500, 4, seed = 1); g2500 <- make_grid_fixture(2500, 4, seed = 2)
gint <- g500; for (v in paste0("v", 1:4)) gint[[v]] <- round(gint[[v]] * 2) / 2   # exact ties
gisl <- make_grid_fixture(600, 4, seed = 3, islands = 3)
va <- c("mgmt","business","sciences","health","education","arts","sales","trades","resources","manufacturing")
cases <- list(
  nc_skater      = quote(skater(nc, c("SID74","SID79"), n_regions = 6)),
  nc_skater_fl   = quote(skater(nc, c("SID74","SID79"), n_regions = 6, floor = "BIR74", floor_value = 50000)),
  nc_ward        = quote(ward_spatial(nc, c("SID74","SID79"), n_regions = 6)),
  nc_azp_basic   = quote(azp(nc, c("SID74","SID79"), n_regions = 6, method = "basic", seed = 1)),
  nc_azp_sa      = quote(azp(nc, c("SID74","SID79"), n_regions = 6, method = "sa", seed = 1)),
  nc_azp_tabu    = quote(azp(nc, c("SID74","SID79"), n_regions = 6, method = "tabu", seed = 1)),
  nc_maxp        = quote(max_p_regions(nc, c("SID74","SID79"), threshold_var = "BIR74", threshold = 100000, seed = 3)),
  meck_skater    = quote(skater(meck, c("medinc","pct_bach","rent","age"), n_regions = 10)),
  meck_skater_fl = quote(skater(meck, c("medinc","pct_bach","rent","age"), n_regions = 10, floor = "pop", floor_value = 50000)),
  meck_ward      = quote(ward_spatial(meck, c("medinc","pct_bach","rent","age"), n_regions = 10)),
  meck_azp_sa    = quote(azp(meck, c("medinc","pct_bach","rent","age"), n_regions = 10, method = "sa", seed = 2)),
  meck_maxp      = quote(max_p_regions(meck, c("medinc","pct_bach","rent","age"), threshold_var = "pop", threshold = 50000, seed = 1)),
  van_skater     = quote(skater(van, va, n_regions = 10)),
  van_skater_fl  = quote(skater(van, va, n_regions = 12, floor = "pop", floor_value = 30000)),
  van_ward       = quote(ward_spatial(van, va, n_regions = 10)),
  van_azp_tabu   = quote(azp(van, va, n_regions = 10, method = "tabu", seed = 1)),
  g500_skater    = quote(skater(g500, paste0("v",1:4), n_regions = 12)),
  g500_skater_fl = quote(skater(g500, paste0("v",1:4), n_regions = 12, floor = "pop", floor_value = 20000)),
  g500_ward      = quote(ward_spatial(g500, paste0("v",1:4), n_regions = 12)),
  g500_azp_sa    = quote(azp(g500, paste0("v",1:4), n_regions = 12, method = "sa", seed = 5)),
  g500_azp_basic = quote(azp(g500, paste0("v",1:4), n_regions = 12, method = "basic", seed = 5)),
  g500_maxp      = quote(max_p_regions(g500, paste0("v",1:4), threshold_var = "pop", threshold = 20000, seed = 2)),
  g2500_skater   = quote(skater(g2500, paste0("v",1:4), n_regions = 20)),
  g2500_skater_fl= quote(skater(g2500, paste0("v",1:4), n_regions = 20, floor = "pop", floor_value = 40000)),
  g2500_azp_sa   = quote(azp(g2500, paste0("v",1:4), n_regions = 20, method = "sa", seed = 1, max_iterations = 200)),
  gint_skater    = quote(skater(gint, paste0("v",1:4), n_regions = 12)),
  gint_ward      = quote(ward_spatial(gint, paste0("v",1:4), n_regions = 12)),
  gisl_skater_br = quote(skater(gisl, paste0("v",1:4), n_regions = 8, bridge_islands = TRUE)),
  gisl_nb_bridge = quote(suppressMessages(spopt:::prepare_weights(gisl, "queen", bridge_islands = TRUE)))
)
same_partition <- function(a, b) { if (length(a) != length(b)) return(FALSE); tb <- table(a, b); all(rowSums(tb > 0) == 1) && all(colSums(tb > 0) == 1) }
# Feature guards, so the harness runs on any subset of the regionalization PRs:
# AZP/max-p cases need the rewritten search engine (marked by azp()'s
# initial_regions argument); REDCAP cases need redcap() to exist.
has_engine <- "initial_regions" %in% names(formals(spopt::azp))
has_redcap <- exists("redcap", envir = asNamespace("spopt"), inherits = FALSE)
needs_ok <- function(nm) {
  if (grepl("azp|maxp", nm) && !has_engine) return(FALSE)
  if (grepl("redcap", nm) && !has_redcap) return(FALSE)
  TRUE
}
dir.create(dir, showWarnings = FALSE, recursive = TRUE); fails <- 0; skipped <- 0
for (nm in names(cases)) {
  f <- file.path(dir, paste0(nm, ".rds"))
  if (!needs_ok(nm) || (mode == "compare" && !file.exists(f))) { skipped <- skipped + 1; cat(sprintf("%-16s skipped (%s)\n", nm, if (!needs_ok(nm)) "function or argument not in this build" else "no baseline file")); next }
  t0 <- proc.time()[["elapsed"]]; res <- suppressMessages(eval(cases[[nm]])); secs <- proc.time()[["elapsed"]] - t0
  if (mode == "capture") { saveRDS(res, f); cat(sprintf("captured %-16s %6.2fs\n", nm, secs)) } else {
    ref <- readRDS(f)
    if (inherits(res, "nb")) { ok <- identical(unclass(res), unclass(ref)); msg <- if (ok) "identical nb" else "nb DIFFERS" } else {
      ma <- attr(ref, "spopt"); mb <- attr(res, "spopt")
      ok_part <- same_partition(ref$.region, res$.region)
      ok_obj <- isTRUE(all.equal(ma$objective, mb$objective, tolerance = 1e-9))
      ok_k <- identical(ma$n_regions, mb$n_regions)
      ok <- ok_part && ok_obj && ok_k
      msg <- sprintf("partition %s, objective %s (%.6g vs %.6g), k %s", if (ok_part) "same" else "DIFFERS", if (ok_obj) "same" else "DIFFERS", ma$objective, mb$objective, if (ok_k) "same" else "DIFFERS") }
    if (!ok) fails <- fails + 1
    cat(sprintf("%-16s %6.2fs  %s\n", nm, secs, msg)) }
}
if (mode == "compare") cat(if (fails == 0) sprintf("ALL EQUIVALENT (%d skipped)\n", skipped) else sprintf("%d CASE(S) DIFFER (%d skipped)\n", fails, skipped))
# Note: gisl_nb_bridge and gisl_skater_br were re-captured after the
# bridge_islands connectivity fix (commit "Fix bridge_islands leaving the graph
# disconnected"); the pre-fix build left that fixture in two components.
