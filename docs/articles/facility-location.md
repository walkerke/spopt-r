# Solving facility location problems with spopt

*Facility location* problems ask: given a set of demand points and
candidate facility sites, where should we place facilities to best serve
demand? These problems arise across industries - siting fire stations to
minimize emergency response times, placing retail stores to maximize
customer coverage, or locating warehouses to minimize shipping
distances.

spopt implements several classic facility location algorithms, each
optimizing a different objective. This vignette demonstrates these
algorithms using Census data from Tarrant County, Texas (home of Fort
Worth), showing how to select optimal locations for public services.

## Setting up the problem

Facility location problems have three core components:

1.  **Demand points**: Locations that need to be served (e.g.,
    population centers, customer addresses)
2.  **Candidate facilities**: Potential sites where facilities could be
    built
3.  **Cost matrix**: Distance or travel time between each demand point
    and candidate facility

A key distinction: in practice, you typically have *many* demand points
but *few* candidate facility sites. For example, you might need to serve
500 Census tracts but only have 25 potential building sites to choose
from. This asymmetry is what makes the problem tractable - the solver
selects from a limited set of candidates rather than considering every
possible location.

Let’s set up a realistic scenario: we want to place community health
centers to serve the population of Tarrant County. We’ll use Census
tract centroids as demand points, and sample 30 candidate locations from
across the county to represent potential facility sites.

``` r

library(spopt)
library(tidycensus)
library(tidyverse)
library(sf)
library(mapgl)

# Get tract-level population data
tarrant <- get_acs(
  geography = "tract",
  variables = "B01003_001",
  state = "TX",
  county = "Tarrant",
  geometry = TRUE,
  year = 2023
) |>
  filter(estimate > 0) |>
  rename(population = estimate)

# Demand points: all tract centroids
demand_pts <- tarrant |>
  st_centroid()

# Candidate facilities: sample 30 locations across the county
# In practice, these might be specific parcels, existing buildings, or zoned commercial sites
set.seed(1983)
n_candidates <- 30

county_boundary <- tarrant |> st_union()
candidate_pts <- st_sample(county_boundary, n_candidates) |>
  st_as_sf() |>
  mutate(id = row_number())
```

We now have 448 demand points (tract centroids) and 30 candidate
facility locations. This setup mirrors real-world planning where you’re
evaluating a shortlist of potential sites.

``` r

# Visualize the setup
maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_circle_layer(
    id = "demand",
    source = demand_pts,
    circle_color = "steelblue",
    circle_radius = 3,
    circle_opacity = 0.5
  ) |>
  add_circle_layer(
    id = "candidates",
    source = candidate_pts,
    circle_color = "black",
    circle_radius = 6,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

## P-Median: Minimizing total distance

The *P-Median* problem ([Hakimi 1964](#ref-hakimi1964)) minimizes the
total weighted distance from demand points to their assigned facilities.
This is the classic efficiency-focused location model - it finds
locations that minimize how far people, on average, must travel.

``` r

result_pmedian <- p_median(
  demand = demand_pts,
  facilities = candidate_pts,
  n_facilities = 5,
  weight_col = "population"
)
```

The solver runs quickly with 30 candidates - the optimization scales
with the number of candidate sites, not demand points. Each demand point
now has a `.facility` column (the row of its assigned facility) and a
`.cost` column (the distance to it, here in meters).

To visualize the allocation,
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
draws a line from each demand point to the facility that serves it.
These are sometimes called *allocation lines* or *desire lines*.

``` r

# One line per tract, from the tract to its assigned facility
pmedian_lines <- spider_lines(result_pmedian) |>
  mutate(facility = as.character(facility_index))

# Get selected facility locations
selected <- result_pmedian$facilities |>
  filter(.selected) |>
  mutate(id = as.character(id))

facility_colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")

# Map the results
maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_line_layer(
    id = "allocations",
    source = pmedian_lines,
    line_color = match_expr(
      column = "facility",
      values = selected$id,
      stops = facility_colors
    ),
    line_width = 1,
    line_opacity = 0.6
  ) |>
  add_circle_layer(
    id = "facilities",
    source = selected,
    circle_color = match_expr(
      column = "id",
      values = selected$id,
      stops = facility_colors
    ),
    circle_radius = 10,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Each line connects a tract to its assigned facility, colored by
facility. The solution minimizes the total population-weighted distance.
Along with the geometry,
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
returns each link’s `weight` (the tract’s population) and `cost` (its
distance), so you can summarize the allocation directly:

``` r

pmedian_lines |>
  st_drop_geometry() |>
  group_by(facility) |>
  summarize(
    tracts = n(),
    population = sum(weight),
    mean_km = weighted.mean(cost, weight) / 1000
  )
```

    # A tibble: 5 × 4
      facility tracts population mean_km
      <chr>     <int>      <dbl>   <dbl>
    1 14           59     330421    6.85
    2 18          108     430427    7.72
    3 27           81     429709    8.78
    4 28          101     496283    8.39
    5 7            99     448903    7.01

The `facility_index` and `demand_index` columns refer to rows of
`result_pmedian$facilities` and `result_pmedian$demand`, so you can join
any other attributes back on.

You can access solution metadata through the `spopt` attribute:

``` r

attr(result_pmedian, "spopt")
```

    $algorithm
    [1] "p_median"

    $weight_col
    [1] "population"

    $n_selected
    [1] 5

    $n_facilities
    [1] 5

    $objective
    [1] 16669528322

    $mean_distance
    [1] 7805.025

    $n_fixed
    [1] 0

    $fixed_col
    NULL

    $fixed_facilities
    NULL

    $solve_time
    [1] 0.514914

The `objective` value represents the total weighted distance - lower is
better.

## P-Center: Minimizing maximum distance

While P-Median optimizes *overall* accessibility, the *P-Center* problem
([Hakimi 1965](#ref-hakimi1965)) focuses on *equity* - it minimizes the
maximum distance any demand point must travel. This is critical for
emergency services where we need to guarantee that *everyone* is within
a reasonable distance.

``` r

result_pcenter <- p_center(
  demand = demand_pts,
  facilities = candidate_pts,
  n_facilities = 5
)

selected_pcenter <- result_pcenter$facilities |>
  filter(.selected)

# Compare to P-Median locations
maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_circle_layer(
    id = "pmedian",
    source = selected,
    circle_color = "#3498db",
    circle_radius = 10,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  ) |>
  add_circle_layer(
    id = "pcenter",
    source = selected_pcenter,
    circle_color = "#e74c3c",
    circle_radius = 5,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Notice how the P-Center solution (red) pushes facilities toward the
edges of the county to ensure no one is too far away, while P-Median
(blue) concentrates facilities where population is densest. Despite
this, two common facilities are selected between the solutions.

## MCLP: Maximum coverage with limited facilities

The *Maximum Coverage Location Problem* (MCLP) ([Church and ReVelle
1974](#ref-church1974)) maximizes the demand covered within a service
radius when you can only build a fixed number of facilities. This is
useful when you have budget constraints but want to cover as many people
as possible.

``` r

result_mclp <- mclp(
  demand = demand_pts,
  facilities = candidate_pts,
  n_facilities = 5,
  service_radius = 5000,  # 5 km
  weight_col = "population"
)

# Calculate coverage
covered_pop <- result_mclp$demand |>
  filter(.covered) |>
  pull(population) |>
  sum()

total_pop <- sum(demand_pts$population, na.rm = TRUE)

cat(sprintf("Coverage: %s of %s (%.1f%%)",
            format(covered_pop, big.mark = ","),
            format(total_pop, big.mark = ","),
            100 * covered_pop / total_pop))
```

    Coverage: 520,623 of 2,135,743 (24.4%)

The `service_radius` parameter defines what “covered” means - any demand
point within this distance of a selected facility is considered covered.
The algorithm then selects facilities to maximize the total covered
population.

Covered demand points are assigned to their nearest selected facility,
so
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
draws lines only for covered tracts. Uncovered tracts have no assignment
and get no line:

``` r

mclp_lines <- spider_lines(result_mclp)

uncovered <- result_mclp$demand |>
  filter(!.covered)

maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_line_layer(
    id = "covered",
    source = mclp_lines,
    line_color = "#2c7fb8",
    line_width = 1
  ) |>
  add_circle_layer(
    id = "uncovered",
    source = uncovered,
    circle_color = "#d7301f",
    circle_radius = 3,
    circle_opacity = 0.7
  ) |>
  add_circle_layer(
    id = "facilities",
    source = filter(result_mclp$facilities, .selected),
    circle_color = "#2c7fb8",
    circle_radius = 10,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Red points are tracts outside the 5 km service radius of every selected
facility.

## LSCP: Minimum facilities for full coverage

The *Location Set Covering Problem* (LSCP) ([Toregas et al.
1971](#ref-toregas1971)) asks the opposite question from MCLP: what’s
the minimum number of facilities needed to cover *all* demand within a
service radius?

``` r

result_lscp <- lscp(
  demand = demand_pts,
  facilities = candidate_pts,
  service_radius = 8000  # 8 km
)

n_selected <- sum(result_lscp$facilities$.selected)
cat(sprintf("Minimum facilities needed for full coverage: %d", n_selected))
```

    Minimum facilities needed for full coverage: 19

The solver finds 19 facilities in this example. LSCP is particularly
useful for planning emergency services where coverage is mandatory -
every resident must be within a certain response time of a fire station
or hospital.

## P-Dispersion: Spreading facilities apart

Most facility location problems assume demand points want to be *close*
to facilities. But some facilities are *obnoxious* - landfills, prisons,
or polluting industries that communities want far away. The
*P-Dispersion* problem ([Kuby 1987](#ref-kuby1987)) maximizes the
minimum distance between facilities.

P-Dispersion is also useful for environmental monitoring networks or
cell tower placement where you want sensors spread across a region.

``` r

result_pdispersion <- p_dispersion(
  facilities = candidate_pts,
  n_facilities = 10
)

selected_disp <- result_pdispersion |>
  filter(.selected)

maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_circle_layer(
    id = "facilities",
    source = selected_disp,
    circle_color = "#2ecc71",
    circle_radius = 10,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Notice how the facilities are spread around the county’s perimeter and
interior - maximizing the minimum inter-facility distance.

## CFLP: Capacitated facility location

Real facilities have capacity limits - a clinic can only see so many
patients per day, a warehouse can only store so much inventory. The
*Capacitated Facility Location Problem* (CFLP) ([Daskin
2013](#ref-daskin2013); [Sridharan 1995](#ref-sridharan1995)) adds
capacity constraints and allows demand to be split across multiple
facilities.

### Varying capacity limits

In practice, candidate sites often have different capacities based on
lot size, zoning, or building constraints. Let’s simulate a realistic
scenario with small, medium, and large sites:

``` r

# Create candidates with varying capacities
set.seed(1983)
candidate_facilities <- candidate_pts |>
  mutate(
    # Assign site sizes: small (200k), medium (400k), large (800k)
    site_type = sample(c("small", "medium", "large"), n(), replace = TRUE),
    capacity = case_when(
      site_type == "small" ~ 200000,
      site_type == "medium" ~ 400000,
      site_type == "large" ~ 800000
    )
  )

result_cflp <- cflp(
  demand = demand_pts,
  facilities = candidate_facilities,
  n_facilities = 5,
  weight_col = "population",
  capacity_col = "capacity"
)

# Check which sites were selected and their utilization
result_cflp$facilities |>
  filter(.selected) |>
  st_drop_geometry() |>
  select(id, site_type, capacity, .utilization)
```

      id site_type capacity .utilization
    1  5    medium    4e+05    0.9712625
    2  7    medium    4e+05    1.0000000
    3 12    medium    4e+05    1.0000000
    4 21     large    8e+05    0.6974662
    5 27    medium    4e+05    0.9731625

The solver selects facilities that together can serve all demand while
minimizing total distance. Notice that `.utilization` shows what
fraction of each facility’s capacity is used.

When demand exceeds capacity at the nearest facility, the solver splits
demand across multiple facilities:

``` r

# How many demand points are split?
n_split <- sum(result_cflp$demand$.split)
cat(sprintf("%d of %d demand points are served by multiple facilities",
            n_split, nrow(demand_pts)))
```

    2 of 448 demand points are served by multiple facilities

By default,
[`spider_lines()`](https://walker-data.com/spopt-r/reference/spider_lines.md)
draws one line per demand point, to its *primary* facility (the one
serving the largest share). Use `allocations = "all"` to draw every
allocation, so a split demand point gets one line per facility. The
`share` column gives the fraction of that point’s demand sent along each
line, and `weight` gives the population it carries:

``` r

cflp_lines <- spider_lines(result_cflp, allocations = "all") |>
  mutate(split = if_else(demand_index %in% which(result_cflp$demand$.split),
                         "split", "single"))

maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_line_layer(
    id = "single",
    source = filter(cflp_lines, split == "single"),
    line_color = "#969696",
    line_width = 0.75,
    line_opacity = 0.5
  ) |>
  add_line_layer(
    id = "split",
    source = filter(cflp_lines, split == "split"),
    line_color = "#d7301f",
    line_width = interpolate(
      column = "share",
      values = c(0, 1),
      stops = c(1.5, 5)
    )
  ) |>
  add_circle_layer(
    id = "facilities",
    source = filter(result_cflp$facilities, .selected),
    circle_color = "black",
    circle_radius = 8,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Red lines belong to demand points split across facilities, with line
width showing the share of demand sent along each line. Gray lines are
demand points served entirely by one facility.

### Incorporating facility costs

What if different sites have different costs - perhaps due to real
estate prices, construction costs, or lease rates? CFLP can incorporate
these costs to find the economically optimal solution.

When you provide a `facility_cost_col` and set `n_facilities = 0`, the
solver determines the optimal *number* of facilities by balancing fixed
costs (opening facilities) against variable costs (transportation
distance). The key is scaling costs appropriately - fixed costs should
be comparable to total transportation costs:

``` r

# Add costs based on site size
# Scale to be comparable with total transport costs (population * distance)
candidate_with_costs <- candidate_facilities |>
  mutate(
    # Fixed cost to open each facility (scaled to match transport cost units)
    fixed_cost = case_when(
      site_type == "small" ~ 5e8,   # Higher cost per unit capacity
      site_type == "medium" ~ 8e8,
      site_type == "large" ~ 1e9
    )
  )

result_with_costs <- cflp(
  demand = demand_pts,
  facilities = candidate_with_costs,
  n_facilities = 0,  # Let solver determine optimal number
  weight_col = "population",
  capacity_col = "capacity",
  facility_cost_col = "fixed_cost"
)

# How many facilities does the cost-optimized solution select?
cost_meta <- attr(result_with_costs, "spopt")
cat(sprintf("Optimal number of facilities: %d\n", cost_meta$n_selected))
```

    Optimal number of facilities: 9

This mode is powerful for budget planning: instead of arbitrarily
choosing 5 facilities, let the optimization determine the
cost-minimizing configuration.

### Real estate costs example

Real estate prices typically vary spatially - central locations cost
more than peripheral ones. Let’s create a realistic scenario where costs
increase toward the county centroid:

``` r

# Calculate distance from county centroid (proxy for "centrality")
county_centroid <- st_centroid(county_boundary)

# Compute distances to center
dist_to_center <- as.numeric(st_distance(candidate_facilities, county_centroid))

candidate_with_realestate <- candidate_facilities |>
  mutate(
    dist_to_center = dist_to_center,
    # Costs higher near center, lower at periphery
    # Normalize to 0-1 range and invert (closer = higher cost)
    centrality = 1 - (dist_to_center - min(dist_to_center)) /
                     (max(dist_to_center) - min(dist_to_center)),
    # Cost ranges from $200-$500 per sqft based on location
    cost_per_sqft = 200 + centrality * 300,
    sqft = capacity / 10,  # Assume 10 people per sqft capacity
    # Scale fixed cost to be comparable with transport costs
    fixed_cost = sqft * cost_per_sqft * 50
  )
```

Let’s visualize how costs vary across the county:

``` r

maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_circle_layer(
    id = "candidates",
    source = candidate_with_realestate,
    circle_color = interpolate(
      column = "cost_per_sqft",
      values = c(200, 350, 500),
      stops = c("#2166ac", "#f7f7f7", "#b2182b")
    ),
    circle_radius = interpolate(
      column = "capacity",
      values = c(200000, 500000, 800000),
      stops = c(6, 10, 14)
    ),
    circle_stroke_color = "white",
    circle_stroke_width = 1
  )
```

Red sites are expensive (central), blue sites are cheaper (peripheral).
Larger circles indicate higher capacity.

Now let’s compare two solutions: one ignoring costs (P-Median) and one
incorporating real estate costs:

``` r

# Solution ignoring costs (just minimize distance)
result_no_cost <- p_median(
  demand = demand_pts,
  facilities = candidate_with_realestate,
  n_facilities = 5,
  weight_col = "population"
)

# Solution with real estate costs
result_with_realestate <- cflp(
  demand = demand_pts,
  facilities = candidate_with_realestate,
  n_facilities = 5,
  weight_col = "population",
  capacity_col = "capacity",
  facility_cost_col = "fixed_cost"
)

# Compare selections
selected_no_cost <- result_no_cost$facilities |>
  filter(.selected) |>
  mutate(method = "Distance only")

selected_with_cost <- result_with_realestate$facilities |>
  filter(.selected) |>
  mutate(method = "With real estate costs")
```

``` r

# Side-by-side comparison
maplibre(bounds = tarrant) |>
  add_fill_layer(
    id = "tracts",
    source = tarrant,
    fill_color = "lightgray",
    fill_opacity = 0.3
  ) |>
  add_circle_layer(
    id = "no_cost",
    source = selected_no_cost,
    circle_color = "#e41a1c",
    circle_radius = 12,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  ) |>
  add_circle_layer(
    id = "with_cost",
    source = selected_with_cost,
    circle_color = "#377eb8",
    circle_radius = 7,
    circle_stroke_color = "white",
    circle_stroke_width = 2
  )
```

Red circles show the distance-only solution (P-Median), blue circles
show the cost-aware solution. Notice how the cost-aware solution may
shift toward cheaper peripheral sites, trading some accessibility for
lower real estate costs.

## Comparing algorithms

Here’s a quick reference for choosing the right algorithm:

| Algorithm | Objective | Use case |
|----|----|----|
| **P-Median** | Minimize total weighted distance | General efficiency; warehouses, retail |
| **P-Center** | Minimize maximum distance | Equity-focused; emergency services |
| **MCLP** | Maximize covered demand | Fixed budget; expand coverage |
| **LSCP** | Minimize facilities for full coverage | Mandatory coverage; emergency planning |
| **P-Dispersion** | Maximize minimum inter-facility dist. | Obnoxious facilities; monitoring networks |
| **CFLP** | Minimize cost with capacity limits | Logistics; realistic capacity constraints |

For most public service planning, start with **P-Median** for an
efficiency baseline, then compare to **P-Center** to see the equity
trade-offs. If you have a hard budget constraint, **MCLP** helps
understand what coverage is achievable with limited facilities.

## Candidate site selection strategies

The examples above used randomly sampled candidate locations. In
practice, you’d generate candidates more thoughtfully:

- **Existing infrastructure**: Current facility locations that could be
  expanded
- **Zoning-based**: Sites zoned for commercial or institutional use
- **Network nodes**: Major intersections or highway interchanges
- **Population centers**: Centroids of high-density areas
- **Grid sampling**: Regular grid across the study area for exploratory
  analysis

You can also combine strategies - start with a coarse grid of candidates
for initial analysis, then refine to specific parcels for detailed
planning.

## Using custom cost matrices

By default, spopt calculates Euclidean distances between points. But
real-world accessibility depends on road networks and travel times, not
straight-line distance. All facility location functions accept a
`cost_matrix` parameter for custom distances.

See the [Travel-Time Cost
Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md)
vignette for how to generate travel time matrices using r5r, and then
pass them to these functions.

``` r

# Example with custom cost matrix
cost_mat <- my_travel_time_matrix  # Generated from r5r or similar

result <- p_median(
  demand = demand_pts,
  facilities = candidate_pts,
  n_facilities = 5,
  weight_col = "population",
  cost_matrix = cost_mat
)
```

## Next steps

- [Regionalization](https://walker-data.com/spopt-r/articles/regionalization.md) -
  Build spatially-contiguous regions
- [Huff Model](https://walker-data.com/spopt-r/articles/huff-model.md) -
  Model market share and retail competition
- [Travel-Time Cost
  Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md) -
  Use real-world travel times with r5r

## References

Church, R., and C. ReVelle. 1974. “The Maximal Covering Location
Problem.” *Papers in Regional Science* 32 (1): 101–18.
<https://doi.org/10.1007/BF01942293>.

Daskin, M. S. 2013. *Network and Discrete Location: Models, Algorithms,
and Applications*. 2nd ed. John Wiley & Sons.
<https://doi.org/10.1002/9781118537015>.

Hakimi, S. L. 1964. “Optimum Locations of Switching Centers and the
Absolute Centers and Medians of a Graph.” *Operations Research* 12 (3):
450–59. <https://doi.org/10.1287/opre.12.3.450>.

Hakimi, S. L. 1965. “Optimum Distribution of Switching Centers in a
Communication Network and Some Related Graph Theoretic Problems.”
*Operations Research* 13 (3): 462–75.
<https://doi.org/10.1287/opre.13.3.462>.

Kuby, M. J. 1987. “Programming Models for Facility Dispersion: The
p-Dispersion and Maxisum Dispersion Problems.” *Geographical Analysis*
19 (4): 315–29. <https://doi.org/10.1111/j.1538-4632.1987.tb00133.x>.

Sridharan, R. 1995. “The Capacitated Plant Location Problem.” *European
Journal of Operational Research* 87 (2): 203–13.
<https://doi.org/10.1016/0377-2217(95)00042-O>.

Toregas, C., R. Swain, C. ReVelle, and L. Bergman. 1971. “The Location
of Emergency Service Facilities.” *Operations Research* 19 (6): 1363–73.
<https://doi.org/10.1287/opre.19.6.1363>.
