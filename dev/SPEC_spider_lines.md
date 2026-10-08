# Spec: Spider Lines for Location-Allocation Results

**Author**: Kyle Walker
**Date**: 2026-10-08
**Package**: spopt-r
**Status**: Implemented on branch `spider-lines` (v3; Codex spec + implementation reviews 2026-10-08). Deferred: vignette "Mapping allocations" section and the r5r adapter helper.
**Depends on**: existing `p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, `huff()` result objects

---

## Motivation

Location-allocation results are easiest to read as *spider lines*: one line from each demand location to the facility that serves it. ArcGIS Network Analyst's location-allocation layer produces these as "Lines", and the retired `tbart` package drew them as star diagrams. spopt users currently build them by hand from the `.facility` column, which takes about 15 lines of boilerplate and mishandles polygon demand, split allocations (`cflp()`), and probabilistic assignment (`huff()`).

Spider lines are also the most effective visual for demos and teaching: they show at a glance which facility serves which area and how far each service distance stretches.

Goal: one function that turns any spopt location-allocation result into an `sf` line layer, drawn as straight lines by default. An optional hook lets users plug in a router (OSRM, r5r, dodgr, ...) so lines follow the network that built their travel-time matrix.

---

## Current state (what the solvers already return)

All of the following return `list(demand = <sf>, facilities = <sf>)` with `attr(x, "spopt")` metadata:

| Function | Class | Assignment on `demand` | Extra |
|---|---|---|---|
| `p_median()` | `spopt_pmedian`, `spopt_locate` | `.facility` (1-based) | |
| `p_center()` | `spopt_pcenter`, `spopt_locate` | `.facility` (1-based). With `method = "binary_search"`, **`0L` for every row when no feasible solution is found** (`R/highs-backend.R:328`); the MIP method errors on non-optimal status instead (`R/highs-backend.R:453`) | no weight_col |
| `cflp()` | `spopt_cflp`, `spopt_locate` | `.facility` = primary assignment, `.split` | `attr(.,"spopt")$allocation_matrix`: dense n_demand x n_fac raw solver fractions, which may hold residue below 1e-6 on unselected facilities (`tests/testthat/test-solvers-regression.R:182`) |
| `lscp()` | `spopt_lscp`, `spopt_locate` | `.facility` (nearest selected; NA if uncovered), `.covered` | no weight_col |
| `mclp()` | `spopt_mclp`, `spopt_locate` | `.facility` (nearest selected; NA if uncovered), `.covered` | |

`huff()` returns `list(demand, stores, probability_matrix)` with class `spopt_huff` (not `spopt_locate`). `demand$.primary_store` is a **1-based index** into `stores`. Its roxygen calls this an "ID", and the `.prob_<id>` columns use `stores$id` when present. Probabilities can be exactly zero, for example for stores with `Inf` cost (`R/huff.R:201`).

Other facts that matter here:
- **Demand and facilities can have different CRSs.** When `cost_matrix` is supplied, solvers skip distance computation and return both inputs unchanged, with no CRS check (`R/locate_p_median.R:105`).
- **Solvers coerce weights with `as.numeric()`.** The location solvers replace NA/Inf costs internally; `huff()` replaces NA only and keeps Inf (`R/huff.R:201`). So a cost matrix the user passes in may differ from what the optimizer actually used.
- No result currently records `weight_col`, the cost matrix, or the cost of each assignment.

Not in scope: `frlm()` (flow-based) and `p_dispersion()` (no demand).

---

## Proposed API

```r
spider_lines(
  result,
  allocations = c("primary", "all"),
  min_share = 0,
  cost_matrix = NULL,
  anchor = c("centroid", "point_on_surface"),
  route_fun = NULL
)
```

**Parameters**

- `result`: an object of class `spopt_locate` or `spopt_huff`, **with demand and facility rows in the order the solver returned them**. Indices and matrices are positional, so this invariant is documented and bounds-checked, but reordering cannot be detected.
- `allocations`: which demand–facility links to draw.
  - `"primary"` (default): at most one link per demand row, using `.facility` (locate) or `.primary_store` (huff).
  - `"all"`: every allocation above tolerance. For `cflp()` that is each cell of `allocation_matrix` above `tol` (below). For `huff()` it is each store with probability > 0. For the other solvers it is identical to `"primary"`.
- `min_share`: drop any link whose `share` is <= `min_share`. This applies the same way to primary and `"all"` links. Retained shares are **never renormalized**. Must be in `[0, 1)`.
- `cost_matrix`: optional demand x facility matrix, normally the one passed to the solver. Not needed for primary links, which take their cost from the solver's `.cost` column (see "Cost" below). If supplied, it **takes precedence** and supplies `cost` for every emitted link, including additional CFLP/Huff links. Values are reported as-is. Dimensions must match `nrow(demand)` x `nrow(facilities)`.
- `anchor`: how non-point geometries become line endpoints. The default `"centroid"` matches `distance_matrix()` (`st_centroid()`). `"point_on_surface"` keeps endpoints inside concave polygons.
- `route_fun`: `NULL` (straight lines) or a router function (see below).

**Internal tolerance**: `tol = 1e-6` (not user-facing) separates real CFLP allocations from solver numerical residue. A CFLP cell is a link candidate only if it exceeds `tol`; `min_share` then filters on top of that. Huff candidates are all stores with probability > 0 (no extra tolerance), and `min_share` then filters on top of that.

**Share semantics**

| Solver | `share` of a primary link | `share` of `"all"` links |
|---|---|---|
| p_median, p_center, lscp, mclp | 1 | 1 |
| cflp | `allocation_matrix[i, .facility[i]]` (may be < 1 when split) | the matrix cell |
| huff | `probability_matrix[i, .primary_store[i]]` | the matrix cell |

**Valid assignment**: an index counts only if it is a finite integer in `1:nrow(facilities)`. `NA`, `0`, and out-of-range values count as "unassigned". This covers infeasible `p_center()` results and uncovered `lscp`/`mclp` demand.

### Return value

An `sf` object in the demand CRS, one row per link:

| Column | Description |
|---|---|
| `demand_index` | row in `result$demand` |
| `facility_index` | row in `result$facilities` / `result$stores` |
| `share` | see the share table above |
| `weight` | `as.numeric(demand[[weight_col]]) * share`, when the weight column is known; otherwise omitted |
| `cost` | see "Cost" below |
| `length` | `sf::st_length()`, which returns units when the CRS is known and a plain numeric when it is missing |
| `routed` | only when `route_fun` is used: `TRUE` if the router returned a geometry, `FALSE` for a straight-line fallback |

**Geometry type**: LINESTRING for straight output. When `route_fun` is used, all rows (fallbacks included) are cast to MULTILINESTRING so the column has one type. Coincident endpoints are kept as zero-length lines.

The output is a plain `sf` with no extra class, so it works directly in ggplot2, mapgl, and `sf` verbs. Users join facility or demand attributes through the index columns.

**Metadata** `attr(out, "spopt")`:

| Field | Definition |
|---|---|
| `algorithm` | copied from `result` |
| `allocations`, `min_share` | as called |
| `n_lines` | rows returned |
| `n_unassigned` | demand rows with no valid solver assignment (solver outcome, *before* any filtering) |
| `n_filtered` | links removed by `min_share` |
| `n_invalid_geometry` | links dropped because an endpoint geometry was empty or invalid |
| `routed` | logical: was `route_fun` used |
| `n_route_requests` | unique coordinate pairs sent to the router |
| `n_route_failures` | of those unique requests, how many failed (expanded output rows are counted by `routed == FALSE`) |

### Cost

Resolution order for the `cost` column:

1. If `cost_matrix` is supplied, use `cost_matrix[demand_index, facility_index]` for every link.
2. Otherwise, if every emitted link is a primary link, use `result$demand$.cost[demand_index]`.
3. Otherwise (`allocations = "all"` produced non-primary CFLP/Huff links, and there is no matrix): **omit the `cost` column** and emit a message telling the user to pass `cost_matrix`. Never recycle a primary cost onto other links, and never return a partly filled `cost` column.
4. Results made before `.cost` existed have no `.cost` column. Without a matrix, `cost` is omitted.

### Solver change: `.cost` on demand

`p_median()`, `p_center()`, `cflp()`, `lscp()`, `mclp()`, and `huff()` add a `.cost` column to `result$demand`:

- **Definition**: the pair cost between the demand row and its assigned facility (`.facility`) or primary store (`.primary_store`), taken from the **original** supplied or internally computed matrix. That means *before* any NA/Inf replacement, so penalty values never show up as reported travel times. It is not multiplied by demand weight or allocation share.
- **Unassigned demand** (uncovered LSCP/MCLP, infeasible p_center) gets `NA`.
- **CFLP**: the cost of the primary assignment, consistent with `.facility` being the primary.
- Implementation: keep a reference to the cost matrix before the solver's NA/Inf handling, and index it with the final assignments.

This makes assignment cost a solver output in its own right. "Which block groups have the longest assigned trips?" can be answered from `result$demand` without building lines, including when the solver computed the matrix internally and the user never saw it.

### Solver change: record `weight_col`

To fill `weight` without asking the user again, add `weight_col` to the metadata of `p_median()`, `cflp()`, and `mclp()`; `huff()` already records `sales_potential_col`. This is additive. Results created before this change lack the field, and `spider_lines()` then omits `weight`. That is tested.

---

## Geometry and CRS handling

- **Demand and facility CRSs differ** (possible via `cost_matrix`): transform facilities to the demand CRS. This is required handling, not just defensive.
- **One CRS missing and the other set**: error.
- **Both CRSs missing** (as in the existing test fixtures): straight lines are allowed, and `length` is a plain numeric. `route_fun` errors with a clear missing-CRS message.
- **Supported endpoint geometries**: POINT, MULTIPOINT, POLYGON, MULTIPOLYGON. Anything other than POINT is reduced by `anchor`. GEOMETRYCOLLECTION and line geometries: error. Z/M dimensions are dropped (`st_zm()`) before building lines.
- **Lon/lat**: lines stay as two-vertex lines. Great-circle densification (`st_segmentize()`) is out of scope and mentioned in the docs.

---

## Network-following lines (`route_fun`)

spopt never sees road geometry: a travel-time matrix from r5r/OSRM/dodgr is just numbers. Rather than depend on one router, `spider_lines()` accepts a user function.

**Contract**

```r
route_fun(from, to)
```

- Inputs: `from` and `to` are `sf` POINT data frames with the same number of rows *n*, in **EPSG:4326**, each with an `id` column (`"1"`, ..., `"n"`). Row *i* of `from` is routed to row *i* of `to`.
- Return value: an `sf` (with `nrow() == n`) or an `sfc` (with `length() == n`), with a **known CRS** and LINESTRING/MULTILINESTRING geometry, in input order.
- Failures: a failed route is an **empty geometry**. (An sfc cannot hold `NA`, so empty geometry is the only failure representation.) Failed routes fall back to straight lines with `routed = FALSE`, and spopt emits one warning with the failure count.
- Malformed returns are errors, not fallbacks: wrong row count, unknown CRS, or a non-line geometry type.
- Batching: spopt **deduplicates by coordinate pair** and calls `route_fun` **once** with all unique pairs, then maps results (and failures) back to every link that shares the pair.
- Router errors: if `route_fun` itself throws, the error propagates.

**Documented recipes** (in `\dontrun{}` and the vignette; exact arguments to be verified against current router versions at implementation time):

```r
# OSRM (local server): one request per pair
osrm_lines <- function(from, to) {
  geoms <- lapply(seq_len(nrow(from)), function(i) {
    r <- tryCatch(osrm::osrmRoute(src = from[i, ], dst = to[i, ]),
                  error = function(e) NULL)
    if (is.null(r)) sf::st_linestring() else sf::st_geometry(r)[[1]]
  })
  sf::st_sfc(geoms, crs = 4326)
}
```

The r5r recipe needs real assembly logic. `detailed_itineraries()` returns one row per itinerary *leg*, so the adapter must order the legs, union them into one geometry per O-D pair, and restore empty geometries for pairs with no itinerary. Provide it as a tested helper snippet, not a one-liner.

Expected performance: straight lines for 50k links in under about 1 s. Routed performance is bounded by the router.

---

## Internal design

```
spider_lines()
  ├─ .spider_validate(result, cost_matrix, min_share)   class, dims, CRS rules
  ├─ .spider_pairs(result, allocations, min_share)      -> data.frame(d, f, share) + counters
  ├─ .spider_anchor(geom, anchor)                       -> sfc POINT (st_zm, reduce non-points)
  ├─ .spider_straight(from_xy, to_xy, crs)              -> sfc LINESTRING (vectorised)
  ├─ .spider_route(from, to, route_fun, crs)            -> dedupe, single call, validate, fallback
  └─ assemble sf + metadata
```

- `.spider_pairs()` dispatches on class:
  - **locate primary**: rows with valid `.facility`; share from the share table.
  - **cflp all**: `which(allocation_matrix > tol, arr.ind = TRUE)`, ordered by demand then facility.
  - **huff primary / all**: `.primary_store` / `which(probability_matrix > 0, arr.ind = TRUE)`.
  - Candidates are enumerated first. `min_share` is applied last to every path, and the removals are counted in `n_filtered`.
- `.spider_straight()`: benchmark a per-row `st_linestring()` loop at 50k links. If it misses the target, build the coordinate matrix and create the sfc in one call.

---

## Edge cases

- Unsupported class (`frlm`, `p_dispersion`, regionalization): an informative error naming the supported functions.
- No valid assignments (infeasible `p_center`, zero-coverage `lscp`/`mclp`): return a 0-row `sf` with the full column schema and CRS, plus a warning.
- `huff` with `"all"` and a tiny `min_share` on large inputs: warn when the output exceeds about 1e6 rows.
- `cost_matrix` with wrong dimensions, or `min_share` outside `[0, 1)`: error.
- Index out of bounds (e.g. a result edited after solving): treated as unassigned, plus a warning suggesting the row order was changed.

---

## Testing (`tests/testthat/test-spider-lines.R`)

1. **Per solver** (p_median, p_center, cflp, lscp, mclp, huff) on small point data: row counts, indices match `.facility`/`.primary_store`, and endpoints equal the input coordinates.
2. **Invalid assignments**: infeasible `p_center` (zeros), and NA for uncovered MCLP/LSCP. Neither produces lines; `n_unassigned` is correct.
3. **CFLP splits**: build a fixture that actually forces a split, since the current regression fixture has none (`test-solvers-regression.R:173`). Check that `"all"` gives more lines than demand, shares are unrenormalized and sum to about 1 per row, and residue below `tol` is excluded.
4. **Share filtering on primary links**: a split CFLP primary below `min_share` is dropped and counted in `n_filtered`, not in `n_unassigned`.
5. **Huff**: the primary share equals the matrix probability; the `"all"` count equals `sum(probability_matrix > min_share)`; zero-probability (`Inf` cost) stores are excluded.
6. **Polygons**: centroid anchor equals `st_centroid()`; `point_on_surface` falls inside; MULTIPOINT is reduced; Z/M input is accepted.
7. **CRS**: lon/lat and projected inputs; demand/facility CRS mismatch is transformed; one CRS missing errors; both missing gives straight lines with numeric `length`, and `route_fun` errors.
8. **Coincident endpoints** produce zero-length lines.
9. **Cost**:
   - `.cost` on each solver equals the original matrix value for the assigned pair, and is `NA` for unassigned rows.
   - With NA/Inf entries in the supplied matrix, `.cost` reports the original value, not the internal penalty.
   - Automatic (`.cost`) and explicit (`cost_matrix`) costs agree for primary links.
   - `"all"` with extra links and no matrix: `cost` is omitted and a message is emitted.
   - A wrong-dimension matrix errors.
   - A pre-change result without `.cost` omits `cost`.
10. **`weight`**: matches `as.numeric(weight) * share`; a pre-change result without `weight_col` omits the column.
11. **`route_fun`** with mock routers: dogleg geometry; an empty-geometry failure; duplicate pairs (checked with a closure call counter: one call, unique pairs only, and a failed duplicated pair falls back for every affected link); malformed returns (wrong length, no CRS, POINT geometry) error; output cast to MULTILINESTRING.
12. **Zero-row output** has the correct schema.

---

## Documentation

- Roxygen page with a runnable straight-line example and `\dontrun{}` router recipes. Explain that these are the "allocation lines" ArcGIS refers to.
- Facility-location vignette: add a short "Mapping allocations" section using the vignette's existing data and solver output, with spider lines in mapgl colored by facility. Use real data; no synthetic examples in vignettes.
- NEWS entry: `spider_lines()`; new `.cost` column on solver demand output; `weight_col` now recorded in result metadata.
- Document `.cost` in each solver's `@return`.
- Fix the `huff()` roxygen to say `.primary_store` is a 1-based index.

---

## Review decisions (Codex review, 2026-10-08; pending author sign-off)

1. **Name**: `spider_lines()`; mention "allocation lines" in the docs.
2. **Cost**: hybrid. Solvers store `.cost` for the primary assignment, and `spider_lines()` uses it automatically. `cost_matrix` is an optional override and the only source of costs for additional CFLP/Huff links. (This replaces the v2 explicit-only decision.)
3. **Huff**: in scope, using matrix probabilities even for primary links.
4. **Router CRS**: a fixed EPSG:4326 contract; require a known input CRS. Adapters handle any other router needs.
5. **Plot defaults**: keep separate. There is no registered `plot.spopt_locate` method today.

---

## Related finding (out of scope)

`tests/testthat/test-p_median.R` and `tests/testthat/test-lscp.R` start with an unconditional `skip("Rust compilation required")`, so they never run. These solvers now use HiGHS, not Rust, so the skips look stale. Re-enable them in a separate change.
