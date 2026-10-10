# Market analysis with the Huff model

The *Huff model* ([Huff 1963](#ref-huff1963), [1964](#ref-huff1964)) is
a classic approach to modeling consumer spatial behavior and market
share. It predicts the probability that a consumer will visit a
particular store based on two factors: the store’s attractiveness (often
measured by size) and the distance the consumer must travel.

This probabilistic model is widely used in retail site selection, market
analysis, and competitive assessment. In this vignette, we’ll use
spopt’s [`huff()`](https://walker-data.com/spopt-r/reference/huff.md)
function to analyze grocery store competition in the Austin, Texas
metropolitan area.

## The Huff model formula

The probability that a consumer at location $`i`$ will visit store $`j`$
is:

``` math
P_{ij} = \frac{A_j^\alpha \cdot D_{ij}^\beta}{\sum_{k=1}^{n} A_k^\alpha \cdot D_{ik}^\beta}
```

Where:

- $`A_j`$ is the attractiveness of store $`j`$ (e.g., square footage)
- $`D_{ij}`$ is the distance from consumer $`i`$ to store $`j`$
- $`\alpha`$ is the attractiveness exponent (typically positive)
- $`\beta`$ is the distance decay exponent (typically negative)

The model assumes consumers trade off attractiveness and distance - a
larger store might attract customers from farther away, but consumers
also prefer closer stores, all else being equal.

## Setting up the analysis

Let’s analyze grocery store competition in the Austin metro area. We’ll
use Census block group data to represent consumer demand locations, and
a set of real store locations.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spopt`](https://walker-data.com/spopt-r/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidycensus`](https://walker-data.com/tidycensus/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`tidyverse`](https://tidyverse.tidyverse.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sf`](https://r-spatial.github.io/sf/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`mapgl`](https://walker-data.com/mapgl/)`)`\
\
`# Get block group data for Travis, Williamson, and Hays counties`\
`austin_counties`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Travis"``, ``"Williamson"``, ``"Hays"``)`\
\
`demand`` ``<-`` `[`get_acs`](https://walker-data.com/tidycensus/reference/get_acs.html)`(`\
`  geography ``=`` ``"block group"``,`\
`  variables ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``pop ``=`` ``"B01003_001"``, income ``=`` ``"B19013_001"``)``,`\
`  state ``=`` ``"TX"``,`\
`  county ``=`` ``austin_counties``,`\
`  geometry ``=`` ``TRUE``,`\
`  year ``=`` ``2023``,`\
`  output ``=`` ``"wide"`\
`)`` ``|>`\
`  `[`st_transform`](https://r-spatial.github.io/sf/reference/st_transform.html)`(``4326``)`` ``|>`\
`  `[`filter`](https://dplyr.tidyverse.org/reference/filter.html)`(``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``incomeE``)``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(`\
`    ``# Spending potential: population weighted by relative income`\
`    spending ``=`` ``popE`` ``*`` ``(``incomeE`` ``/`` `[`median`](https://rspatial.github.io/terra/reference/summarize-generics.html)`(``incomeE``, na.rm ``=`` ``TRUE``)``)`\
`  ``)`

We’ve created a “spending potential” variable that combines population
with income - block groups with higher incomes have more spending
potential per capita. This will serve as our demand weight.

Now let’s define some store locations. For this example, I’ll create a
simulated set of HEB and Whole Foods locations in the Austin area:

\
`stores`` ``<-`` `[`tibble`](https://tibble.tidyverse.org/reference/tibble.html)`(`\
`  id ``=`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"Store_"``, ``1``:``8``)``,`\
`  name ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"HEB Mueller"``, ``"HEB Tech Ridge"``, ``"HEB Hancock"``, ``"HEB South Congress"``,`\
`    ``"Whole Foods Downtown"``, ``"Whole Foods Domain"``, ``"HEB Round Rock"``, ``"HEB Cedar Park"`\
`  ``)``,`\
`  chain ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"HEB"``, ``4``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``"Whole Foods"``, ``2``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``"HEB"``, ``2``)``)``,`\
`  sqft ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``80000``, ``75000``, ``55000``, ``70000``, ``40000``, ``35000``, ``85000``, ``72000``)``,`\
`  lon ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``-``97.7025``, ``-``97.6920``, ``-``97.7215``, ``-``97.7830``,`\
`          ``-``97.7495``, ``-``97.7235``, ``-``97.6790``, ``-``97.8200``)``,`\
`  lat ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``30.2950``, ``30.4420``, ``30.3030``, ``30.2280``,`\
`          ``30.2690``, ``30.4020``, ``30.5080``, ``30.5100``)`\
`)`` ``|>`\
`  `[`st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html)`(``coords ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"lon"``, ``"lat"``)``, crs ``=`` ``4326``)`

## Running the Huff model

The [`huff()`](https://walker-data.com/spopt-r/reference/huff.md)
function calculates market probabilities and expected sales for each
store. We need to provide:

- `demand`: Consumer locations - can be polygons (like block groups) or
  points. For polygons, distances are computed from centroids
  automatically.
- `stores`: Store locations with attractiveness data
- `attractiveness_col`: The column measuring store attractiveness
- `distance_exponent`: Controls how quickly distance reduces
  attractiveness (negative values)

\
`result`` ``<-`` `[`huff`](https://walker-data.com/spopt-r/reference/huff.md)`(`\
`  demand ``=`` ``demand``,`\
`  stores ``=`` ``stores``,`\
`  attractiveness_col ``=`` ``"sqft"``,`\
`  attractiveness_exponent ``=`` ``1.0``,`\
`  distance_exponent ``=`` ``-``1.5``,`\
`  sales_potential_col ``=`` ``"spending"`\
`)`

The function returns a list with two **sf objects** that include all
results:

- `$demand`: The original demand sf with added columns:
  - `.primary_store`: ID of the highest-probability store for each
    location
  - `.entropy`: Competition measure (higher = more competition between
    stores)
  - `.prob_<store_id>`: Probability column for each store
- `$stores`: The original stores sf with added columns:
  - `.market_share`: Proportion of total market captured
  - `.expected_sales`: Expected sales (sum of probability × sales
    potential)

## Visualizing market areas

Let’s map the *primary store* for each block group - the store with the
highest probability of being visited. Since `result$demand` is already
an sf object with `.primary_store` included, we can map it directly:

\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``result``$``demand``)`` ``|>`\
`  `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`    id ``=`` ``"market_areas"``,`\
`    source ``=`` ``result``$``demand``,`\
`    fill_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`      column ``=`` ``".primary_store"``,`\
`      values ``=`` ``1``:``8``,`\
`      stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#e41a1c"``, ``"#377eb8"``, ``"#4daf4a"``, ``"#984ea3"``,`\
`                ``"#ff7f00"``, ``"#a65628"``, ``"#f781bf"``, ``"#999999"``)`\
`    ``)``,`\
`    fill_opacity ``=`` ``0.6`\
`  ``)`` ``|>`\
`  `[`add_line_layer`](https://walker-data.com/mapgl/reference/add_line_layer.html)`(`\
`    id ``=`` ``"borders"``,`\
`    source ``=`` ``result``$``demand``,`\
`    line_color ``=`` ``"white"``,`\
`    line_width ``=`` ``0.3`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"stores"``,`\
`    source ``=`` ``result``$``stores``,`\
`    circle_color ``=`` ``"black"``,`\
`    circle_radius ``=`` ``8``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2`\
`  ``)`

The map shows each block group colored by its primary store. This gives
us a quick sense of each store’s *trade area* - the geographic region
from which it draws most customers. Note how assignments are not an
exclusive function of distance; these probabilities will also be
influenced by store size, which we are using as our attractiveness
metric. In a real-world example, other retail metrics (traffic counts,
foot traffic, etc.) might be used here.

## Understanding market share

The `stores` result contains market share and expected sales for each
location:

\
`result``$``stores`` ``|>`\
`  `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``name``, ``chain``, ``sqft``, ``.expected_sales``, ``.market_share``)`` ``|>`\
`  `[`arrange`](https://dplyr.tidyverse.org/reference/arrange.html)`(`[`desc`](https://dplyr.tidyverse.org/reference/desc.html)`(``.market_share``)``)`

    # A tibble: 8 × 5
      name                 chain        sqft .expected_sales .market_share
      <chr>                <chr>       <dbl>           <dbl>         <dbl>
    1 HEB Round Rock       HEB         85000         397953.        0.168
    2 HEB Cedar Park       HEB         72000         358571.        0.151
    3 HEB South Congress   HEB         70000         353076.        0.149
    4 HEB Tech Ridge       HEB         75000         335597.        0.142
    5 HEB Mueller          HEB         80000         328209.        0.139
    6 HEB Hancock          HEB         55000         234390.        0.0990
    7 Whole Foods Downtown Whole Foods 40000         194714.        0.0822
    8 Whole Foods Domain   Whole Foods 35000         165036.        0.0697

The market share column (`.market_share`) shows each store’s share of
total expected sales across the study area. Larger stores in less
competitive areas tend to capture more market share.

## Competition and entropy

The Huff model also calculates *entropy* for each demand location - a
measure of competition intensity. Entropy is high when consumers have
multiple similarly-attractive options, and low when one store dominates.
This is included directly in `result$demand`:

\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``result``$``demand``)`` ``|>`\
`  `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`    id ``=`` ``"entropy"``,`\
`    source ``=`` ``result``$``demand``,`\
`    fill_color ``=`` ``interpolate``(`\
`      column ``=`` ``".entropy"``,`\
`      values ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``1.5``, ``2``)``,`\
`      stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#2166ac"``, ``"#f7f7f7"``, ``"#b2182b"``)`\
`    ``)``,`\
`    fill_opacity ``=`` ``0.7``,`\
`    tooltip ``=`` ``".entropy"`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"stores"``,`\
`    source ``=`` ``result``$``stores``,`\
`    circle_color ``=`` ``"black"``,`\
`    circle_radius ``=`` ``8``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2`\
`  ``)`

Areas with high entropy (red) have intense competition - consumers there
have multiple good options. Low-entropy areas (blue) are dominated by a
single store; these will often be areas immediately adjacent to a store
location.

## Evaluating new store locations

One of the most powerful applications of the Huff model is *what-if
analysis* - evaluating how a new store would affect market share. Let’s
test a potential new location in south Austin:

\
`# Add a hypothetical new store`\
`new_store`` ``<-`` `[`tibble`](https://tibble.tidyverse.org/reference/tibble.html)`(`\
`  id ``=`` ``"Store_9"``,`\
`  name ``=`` ``"New HEB South Austin"``,`\
`  chain ``=`` ``"HEB"``,`\
`  sqft ``=`` ``70000``,`\
`  lon ``=`` ``-``97.78``,`\
`  lat ``=`` ``30.20`\
`)`` ``|>`\
`  `[`st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html)`(``coords ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"lon"``, ``"lat"``)``, crs ``=`` ``4326``)`\
\
`stores_with_new`` ``<-`` `[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)`(``stores``, ``new_store``)`\
\
`# Re-run the model`\
`result_new`` ``<-`` `[`huff`](https://walker-data.com/spopt-r/reference/huff.md)`(`\
`  demand ``=`` ``demand``,`\
`  stores ``=`` ``stores_with_new``,`\
`  attractiveness_col ``=`` ``"sqft"``,`\
`  attractiveness_exponent ``=`` ``1.0``,`\
`  distance_exponent ``=`` ``-``1.5``,`\
`  sales_potential_col ``=`` ``"spending"`\
`)`\
\
`# Compare market shares`\
`comparison`` ``<-`` ``result``$``stores`` ``|>`\
`  `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
`  `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``name``, original_share ``=`` ``.market_share``)`` ``|>`\
`  `[`left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html)`(`\
`    ``result_new``$``stores`` ``|>`\
`      `[`st_drop_geometry`](https://r-spatial.github.io/sf/reference/st_geometry.html)`(``)`` ``|>`\
`      `[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``name``, new_share ``=`` ``.market_share``)``,`\
`    by ``=`` ``"name"`\
`  ``)`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``change ``=`` ``new_share`` ``-`` ``original_share``)`\
\
`comparison`

    # A tibble: 8 × 4
      name                 original_share new_share   change
      <chr>                         <dbl>     <dbl>    <dbl>
    1 HEB Mueller                  0.139     0.121  -0.0176
    2 HEB Tech Ridge               0.142     0.131  -0.0105
    3 HEB Hancock                  0.0990    0.0865 -0.0125
    4 HEB South Congress           0.149     0.114  -0.0349
    5 Whole Foods Downtown         0.0822    0.0697 -0.0125
    6 Whole Foods Domain           0.0697    0.0639 -0.00585
    7 HEB Round Rock               0.168     0.158  -0.0105
    8 HEB Cedar Park               0.151     0.141  -0.01000

This analysis shows how the new store would cannibalize sales from
existing locations. Stores near the new location see their market share
decline, while distant stores are not as affected.

Let’s visualize how the new store reshapes market areas:

\
[`maplibre`](https://walker-data.com/mapgl/reference/maplibre.html)`(``bounds ``=`` ``result_new``$``demand``)`` ``|>`\
`  `[`add_fill_layer`](https://walker-data.com/mapgl/reference/add_fill_layer.html)`(`\
`    id ``=`` ``"market_areas"``,`\
`    source ``=`` ``result_new``$``demand``,`\
`    fill_color ``=`` `[`match_expr`](https://walker-data.com/mapgl/reference/match_expr.html)`(`\
`      column ``=`` ``".primary_store"``,`\
`      values ``=`` ``1``:``9``,`\
`      stops ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#e41a1c"``, ``"#377eb8"``, ``"#4daf4a"``, ``"#984ea3"``,`\
`                ``"#ff7f00"``, ``"#a65628"``, ``"#f781bf"``, ``"#999999"``, ``"#66c2a5"``)`\
`    ``)``,`\
`    fill_opacity ``=`` ``0.6`\
`  ``)`` ``|>`\
`  `[`add_line_layer`](https://walker-data.com/mapgl/reference/add_line_layer.html)`(`\
`    id ``=`` ``"borders"``,`\
`    source ``=`` ``result_new``$``demand``,`\
`    line_color ``=`` ``"white"``,`\
`    line_width ``=`` ``0.3`\
`  ``)`` ``|>`\
`  `[`add_circle_layer`](https://walker-data.com/mapgl/reference/add_circle_layer.html)`(`\
`    id ``=`` ``"stores"``,`\
`    source ``=`` ``result_new``$``stores``,`\
`    circle_color ``=`` ``"black"``,`\
`    circle_radius ``=`` ``8``,`\
`    circle_stroke_color ``=`` ``"white"``,`\
`    circle_stroke_width ``=`` ``2`\
`  ``)`

Compare this to the original market area map above - the new store
(teal) carves out its own trade area in south Austin, primarily at the
expense of the nearby HEB South Congress location.

## Adjusting model parameters

The `distance_exponent` and `attractiveness_exponent` parameters
significantly affect results:

**Distance exponent** (typically -0.5 to -3.0):

- Values closer to 0: Consumers willing to travel far for attractive
  stores
- More negative values: Strong preference for nearby stores

**Attractiveness exponent** (typically 0.5 to 2.0):

- Values near 1: Linear relationship between size and attractiveness
- Values \> 1: Larger stores disproportionately attractive
- Values \< 1: Diminishing returns to size

\
`# Test with stronger distance decay`\
`result_steep`` ``<-`` `[`huff`](https://walker-data.com/spopt-r/reference/huff.md)`(`\
`  demand ``=`` ``demand``,`\
`  stores ``=`` ``stores``,`\
`  attractiveness_col ``=`` ``"sqft"``,`\
`  attractiveness_exponent ``=`` ``1.0``,`\
`  distance_exponent ``=`` ``-``2.5``,  ``# Much steeper decay`\
`  sales_potential_col ``=`` ``"spending"`\
`)`

With a steeper distance decay (-2.5 vs -1.5), trade areas become smaller
and more localized. This might be more realistic for
convenience-oriented shopping where consumers prioritize proximity.

## Using composite attractiveness

Store attractiveness isn’t always just about size. You might want to
combine multiple factors - size, brand perception, parking availability,
product selection. spopt supports multiple attractiveness variables
directly:

\
`# Add parking data to stores`\
`stores_extended`` ``<-`` ``stores`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(``parking ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``300``, ``250``, ``150``, ``200``, ``120``, ``100``, ``350``, ``275``)``)`\
\
`# Use both sqft and parking as attractiveness factors`\
`result_multi`` ``<-`` `[`huff`](https://walker-data.com/spopt-r/reference/huff.md)`(`\
`  demand ``=`` ``demand``,`\
`  stores ``=`` ``stores_extended``,`\
`  attractiveness_col ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"sqft"``, ``"parking"``)``,`\
`  attractiveness_exponent ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1.0``, ``0.5``)``,  ``# parking has diminishing returns`\
`  distance_exponent ``=`` ``-``1.5``,`\
`  sales_potential_col ``=`` ``"spending"`\
`)`

The composite attractiveness is computed as
$`A_j = \text{sqft}^{1.0} \times \text{parking}^{0.5}`$. You can also
pre-compute a single attractiveness column if you prefer:

\
`stores_composite`` ``<-`` ``stores`` ``|>`\
`  `[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html)`(`\
`    ``# HEB stores get a 20% brand premium`\
`    brand_factor ``=`` `[`if_else`](https://dplyr.tidyverse.org/reference/if_else.html)`(``chain`` ``==`` ``"HEB"``, ``1.2``, ``1.0``)``,`\
`    attractiveness ``=`` ``sqft`` ``*`` ``brand_factor`\
`  ``)`\
\
`result_composite`` ``<-`` `[`huff`](https://walker-data.com/spopt-r/reference/huff.md)`(`\
`  demand ``=`` ``demand``,`\
`  stores ``=`` ``stores_composite``,`\
`  attractiveness_col ``=`` ``"attractiveness"``,`\
`  distance_exponent ``=`` ``-``1.5``,`\
`  sales_potential_col ``=`` ``"spending"`\
`)`

## Next steps

- [Facility
  Location](https://walker-data.com/spopt-r/articles/facility-location.md) -
  Optimize facility placement
- [Regionalization](https://walker-data.com/spopt-r/articles/regionalization.md) -
  Build custom regions from Census data
- [Travel-Time Cost
  Matrices](https://walker-data.com/spopt-r/articles/travel-time-matrices.md) -
  Use driving times instead of straight-line distances

## References

Huff, D. L. 1963. “A Probabilistic Analysis of Shopping Center Trade
Areas.” *Land Economics* 39 (1): 81–90.
<https://doi.org/10.2307/3144521>.

Huff, D. L. 1964. “Defining and Estimating a Trading Area.” *Journal of
Marketing* 28 (3): 34–38. <https://doi.org/10.1177/002224296402800307>.
