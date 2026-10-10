# Regionalization with spopt

*Regionalization* refers to the process of grouping smaller geographic
units into larger, spatially contiguous regions. Unlike standard
clustering methods, regionalization algorithms enforce that resulting
regions form connected geographic areas - you can’t have a region with
disconnected pieces scattered across the map.

This vignette walks through spopt’s regionalization algorithms using
Census tract data from Dallas, Texas. We’ll explore how to build regions
that minimize internal heterogeneity while maintaining spatial
contiguity and meeting population thresholds.

## When would you use regionalization?

Regionalization solves problems across many fields:

- **Political redistricting**: Building compact, contiguous districts
  that balance population
- **Market segmentation**: Creating sales territories with similar
  customer characteristics
- **Health planning**: Aggregating small-area data while preserving
  spatial relationships
- **Urban planning**: Delineating neighborhoods based on socioeconomic
  similarity
- **Census data analysis**: Addressing differential privacy concerns by
  aggregating blocks into larger areas

## Getting Census data

Let’s start by pulling some demographic data for Census tracts in Dallas
County, Texas. We’ll use the tidycensus package to get population,
median household income, and percentage with a bachelor’s degree -
variables that might define meaningful neighborhood clusters.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spopt`](https://walker-data.com/spopt-r/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidycensus`](https://walker-data.com/tidycensus/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`mapgl`](https://walker-data.com/mapgl/)`)`\
\
`dallas`` ``<-`` `[`get_acs`](https://walker-data.com/tidycensus/reference/get_acs.html)`(`\
`  geography ``=`` ``"tract"``,`\
`  variables ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    pop ``=`` ``"B01003_001"``,`\
`    income ``=`` ``"B19013_001"``,`\
`    bachelors ``=`` ``"DP02_0068P"`\
`  ``)``,`\
`  state ``=`` ``"TX"``,`\
`  county ``=`` ``"Dallas"``,`\
`  geometry ``=`` ``TRUE``,`\
`  year ``=`` ``2023``,`\
`  output ``=`` ``"wide"`\
`)`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``incomeE``)``, ``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``bachelorsE``)``)`

We now have 642 Census tracts with population, income, and education
data. Let’s take a quick look at the geographic distribution of median
household income:

\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``dallas``, column ``=`` ``"incomeE"``)`

The map reveals the familiar spatial pattern of income inequality in
Dallas - higher incomes concentrated in the Park Cities north of
downtown, with lower incomes in the southern part of the county.

## Max-P regionalization

The *Max-P* algorithm ([Duque et al. 2012](#ref-duque2012)) finds the
maximum number of regions such that each region exceeds a specified
threshold while minimizing within-region heterogeneity. This is
particularly useful when you need regions that meet minimum population
requirements for statistical reliability. Recent extensions support
compactness constraints ([Feng et al. 2022](#ref-feng2022)) and improved
efficiency ([Wei et al. 2021](#ref-wei2021)).

Let’s create regions where each must contain at least 50,000 people:

\
`maxp_result`` ``<-`` `[`max_p_regions`](https://walker-data.com/spopt-r/reference/max_p_regions.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  threshold_var ``=`` ``"popE"``,`\
`  threshold ``=`` ``50000``,`\
`  n_iterations ``=`` ``100``,`\
`  seed ``=`` ``1983`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``maxp_result``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

Let’s step through the key parameters:

- `attrs`: The variables used to measure similarity. Tracts with similar
  income and education levels will be grouped together.
- `threshold_var`: The variable that must meet the minimum threshold
  (population in this case).
- `threshold`: Each region must have at least this many people.
- `n_iterations`: The algorithm uses a tabu search heuristic; more
  iterations generally yield better solutions.
- `seed`: For reproducibility, since the algorithm has stochastic
  elements.

The result is an sf object with a new `.region` column indicating each
tract’s assigned region. The algorithm found 41 regions, each with at
least 50,000 residents.

You can access metadata about the solution through the `spopt`
attribute:

\
[`attr`](https://rdrr.io/r/base/attr.html)`(``maxp_result``, ``"spopt"``)`

    $algorithm
    [1] "max_p"

    $n_regions
    [1] 41

    $objective
    [1] 671.3679

    $threshold_var
    [1] "popE"

    $threshold
    [1] 50000

    $region_stats
       region n_areas threshold_sum meets_threshold
    1      39      18         60104            TRUE
    2      21      22         68044            TRUE
    3      19      22         60327            TRUE
    4       1      18         59620            TRUE
    5       6      12         53251            TRUE
    6       4      19         79974            TRUE
    7      11      16         62948            TRUE
    8      13      14         56700            TRUE
    9      15      19         80211            TRUE
    10     12      11         50975            TRUE
    11     35      12         54458            TRUE
    12     37      24        103358            TRUE
    13     33      13         51101            TRUE
    14     25      18         55656            TRUE
    15     34      22         62657            TRUE
    16     38      15         52941            TRUE
    17     18      15         56183            TRUE
    18     27      12         56522            TRUE
    19     24      13         68227            TRUE
    20     17      13         56948            TRUE
    21     40      17         74129            TRUE
    22     10      17         83011            TRUE
    23     26      13         58473            TRUE
    24      3      12         51346            TRUE
    25     28      15         53597            TRUE
    26     16      13         55192            TRUE
    27     41      15         58047            TRUE
    28     31      25         83231            TRUE
    29     20      19         69435            TRUE
    30      5      18         69699            TRUE
    31      7      12         52141            TRUE
    32     22      15         62023            TRUE
    33     32      11         57099            TRUE
    34      2      13         68651            TRUE
    35     29      16         68612            TRUE
    36      8      12         60990            TRUE
    37     23      13         68022            TRUE
    38     30      12         57042            TRUE
    39     36      11         52370            TRUE
    40     14      13         55764            TRUE
    41      9      22         89777            TRUE

    $solve_time
    [1] 0.04264092

    $scaled
    [1] TRUE

    $n_iterations
    [1] 100

    $n_sa_iterations
    [1] 100

    $compact
    [1] FALSE

    $compact_weight
    [1] 0.5

    $homogeneous
    [1] TRUE

    $mean_compactness
    NULL

    $region_compactness
    NULL

### Spatial weights

By default, all regionalization functions use **queen contiguity** - two
tracts are neighbors if they share any boundary point (including
corners). You can also use **rook contiguity**, where tracts must share
an edge to be neighbors:

\
`maxp_rook`` ``<-`` `[`max_p_regions`](https://walker-data.com/spopt-r/reference/max_p_regions.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  threshold_var ``=`` ``"popE"``,`\
`  threshold ``=`` ``50000``,`\
`  weights ``=`` ``"rook"``,`\
`  n_iterations ``=`` ``100``,`\
`  seed ``=`` ``1983`\
`)`

For more control, you can specify weights as a list:

- `list(type = "knn", k = 6)`: K-nearest neighbors (useful for point
  data or ensuring connectivity)
- `list(type = "distance", d = 5000)`: Distance-based weights (units
  match your CRS)

You can also pass an `nb` object created with spdep or spopt’s
[`sp_weights()`](https://walker-data.com/spopt-r/reference/sp_weights.md)
function.

### Compact regions

For applications like sales territories or electoral districts, you may
want regions with compact, regular shapes. The `compact` parameter
optimizes for compactness in addition to attribute homogeneity:

\
`maxp_compact`` ``<-`` `[`max_p_regions`](https://walker-data.com/spopt-r/reference/max_p_regions.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  threshold_var ``=`` ``"popE"``,`\
`  threshold ``=`` ``50000``,`\
`  weights ``=`` ``"rook"``,`\
`  compact ``=`` ``TRUE``,`\
`  compact_weight ``=`` ``0.5``,`\
`  n_iterations ``=`` ``100``,`\
`  seed ``=`` ``1983`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``maxp_compact``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

The `compact_weight` parameter (0 to 1) controls the trade-off between
attribute homogeneity and geometric compactness. Higher values
prioritize compact shapes. The parameter `compact_metric` provides a
choice between two compactness metrics. The default, “centroid
dispersion”, is appropriate for both polygons (e.g., state borders) and
point geometries (e.g., store locations). The alternative option, “NMI”
(normalized moment of inertia), is the original metric proposed by Feng
et al. ([2022](#ref-feng2022)), and is appropriate only for polygons.

## SKATER algorithm

*SKATER* (Spatial K’luster Analysis by Tree Edge Removal) ([Assunção et
al. 2006](#ref-assuncao2006)) takes a different approach. It first
builds a minimum spanning tree connecting all tracts based on their
attribute similarity, then iteratively removes edges to create clusters.
The algorithm is fast and produces spatially coherent regions.

\
`skater_result`` ``<-`` `[`skater`](https://walker-data.com/spopt-r/reference/skater.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``6``,`\
`  seed ``=`` ``1983`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``skater_result``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

SKATER supports a `floor` and `floor_value` parameter if you need
minimum population constraints:

\
`skater_constrained`` ``<-`` `[`skater`](https://walker-data.com/spopt-r/reference/skater.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``6``,`\
`  floor ``=`` ``"popE"``,`\
`  floor_value ``=`` ``150000``,`\
`  seed ``=`` ``1983`\
`)`

## AZP: Automatic Zoning Procedure

The *Automatic Zoning Procedure* (AZP) ([Openshaw
1977](#ref-openshaw1977); [Openshaw and Rao 1995](#ref-openshaw1995))
uses local search optimization with three algorithm variants: basic
(greedy), tabu search, and simulated annealing.

\
`azp_result`` ``<-`` `[`azp`](https://walker-data.com/spopt-r/reference/azp.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``20``,`\
`  method ``=`` ``"tabu"``,`\
`  tabu_length ``=`` ``10``,`\
`  max_iterations ``=`` ``100``,`\
`  seed ``=`` ``1983`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``azp_result``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

The `method` parameter controls which algorithm variant to use:

- `"basic"`: Simple greedy local search (fastest)
- `"tabu"`: Tabu search, which maintains a list of recent moves to avoid
  getting stuck in local optima
- `"sa"`: Simulated annealing, which accepts some worse solutions early
  to explore more of the solution space

For large problems, you may also want to use the simulated annealing
variant:

\
`azp_sa`` ``<-`` `[`azp`](https://walker-data.com/spopt-r/reference/azp.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``20``,`\
`  method ``=`` ``"sa"``,`\
`  cooling_rate ``=`` ``0.85``,`\
`  max_iterations ``=`` ``100``,`\
`  seed ``=`` ``1983`\
`)`

## SPENC: Spatially-Encouraged Spectral Clustering

*SPENC* ([Wolf 2021](#ref-wolf2021)) combines spectral clustering with
spatial constraints. It uses a radial basis function (RBF) kernel to
measure attribute similarity and incorporates spatial connectivity into
the spectral embedding. This approach can find clusters with complex,
non-convex shapes that other methods might miss.

\
`spenc_result`` ``<-`` `[`spenc`](https://walker-data.com/spopt-r/reference/spenc.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``15``,`\
`  gamma ``=`` ``1.0``,`\
`  seed ``=`` ``1983`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``spenc_result``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

The `gamma` parameter controls the RBF kernel bandwidth - higher values
create “tighter” clusters in attribute space.

## Ward spatial clustering

Spatially-constrained *Ward* clustering is a hierarchical method that
only allows merging adjacent clusters. At each step, it merges the pair
of adjacent clusters that minimizes the increase in total within-cluster
variance.

\
`ward_result`` ``<-`` `[`ward_spatial`](https://walker-data.com/spopt-r/reference/ward_spatial.md)`(`\
`  ``dallas``,`\
`  attrs ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"incomeE"``, ``"bachelorsE"``)``,`\
`  n_regions ``=`` ``15`\
`)`\
\
[`maplibre_view`](https://walker-data.com/mapgl/reference/maplibre_view.html)`(``ward_result``, column ``=`` ``".region"``, legend ``=`` ``FALSE``)`

Ward clustering is deterministic (no random seed needed) and tends to
produce compact, roughly equal-sized regions.

## Choosing an algorithm

Each regionalization algorithm has strengths for different scenarios:

| Algorithm | Best for | Key features |
|----|----|----|
| **Max-P** | Population thresholds | Maximizes number of regions meeting constraints |
| **SKATER** | Fast, interpretable results | Tree-based, good for large datasets |
| **AZP** | High-quality solutions | Multiple optimization variants |
| **SPENC** | Complex cluster shapes | Spectral embedding with spatial constraints |
| **Ward** | Deterministic, balanced regions | Hierarchical, no tuning required |

For most applications, I’d recommend starting with **Max-P** if you have
population constraints, or **SKATER** for a quick first pass. If you
want to explore the solution space more thoroughly, try **AZP** with
tabu search or simulated annealing.

## Next steps

- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Solve location-allocation problems
- [Huff Model](https://walker-data.com/spopt-r/articles/huff-model.md) -
  Model market share and retail competition
- [Travel-Time Cost
  Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md) -
  Use real-world travel times

## References

Assunção, R. M., M. C. Neves, G. Câmara, and C. da Costa Freitas. 2006.
“Efficient Regionalization Techniques for Socio-Economic Geographical
Units Using Minimum Spanning Trees.” *International Journal of
Geographical Information Science* 20 (7): 797–811.
<https://doi.org/10.1080/13658810600665111>.

Duque, J. C., L. Anselin, and S. J. Rey. 2012. “The Max-p-Regions
Problem.” *Journal of Regional Science* 52 (3): 397–419.
<https://doi.org/10.1111/j.1467-9787.2011.00743.x>.

Feng, X., S. Rey, and R. Wei. 2022. “The Max-p-Compact-Regions Problem.”
*Transactions in GIS* 26: 717–34. <https://doi.org/10.1111/tgis.12874>.

Openshaw, S. 1977. “A Geographical Solution to Scale and Aggregation
Problems in Region-Building, Partitioning and Spatial Modelling.”
*Transactions of the Institute of British Geographers* 2 (4): 459–72.
<https://doi.org/10.2307/622300>.

Openshaw, S., and L. Rao. 1995. “Algorithms for Reengineering 1991
Census Geography.” *Environment and Planning A* 27 (3): 425–46.
<https://doi.org/10.1068/a270425>.

Wei, R., S. Rey, and E. Knaap. 2021. “Efficient Regionalization for
Spatially Explicit Neighborhood Delineation.” *International Journal of
Geographical Information Science* 35 (1): 135–51.
<https://doi.org/10.1080/13658816.2020.1759806>.

Wolf, L. J. 2021. “Spatially-Encouraged Spectral Clustering: A Technique
for Blending Map Typologies and Regionalization.” *International Journal
of Geographical Information Science* 35 (11): 2356–73.
<https://doi.org/10.1080/13658816.2021.1934475>.
