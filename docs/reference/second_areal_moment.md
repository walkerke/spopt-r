# Second areal moment (i.e., second moment of inertia)

Computes the second moment of area (also known as the second moment of
inertia) for polygon geometries. This is a measure of how the area of a
shape is distributed relative to its centroid.

## Usage

``` r
second_areal_moment(x, project = TRUE)
```

## Arguments

- x:

  An sf object, sfc geometry column, or sfg geometry.

- project:

  Logical. If the geometries have geodetic coordinates, then they will
  be projected using an Albers Equal Area Conic projection centered on
  the data.

## Value

Numeric vector of second areal moments

## Details

The second moment of area is the sum of the inertia across the x and y
axes:

The inertia for the x axis is: \$\$I_x = \frac{1}{12}\sum\_{i=1}^{N}
(x_i y\_{i+1} - x\_{i+1}y_i) (x_i^2 + x_ix\_{i+1} + x\_{i+1}^2)\$\$

While the y axis is in a similar form: \$\$I_y =
\frac{1}{12}\sum\_{i=1}^{N} (x_i y\_{i+1} - x\_{i+1}y_i) (y_i^2 +
y_iy\_{i+1} + y\_{i+1}^2)\$\$

where \\(x_i, y_i\\) is the current point and \\(x\_{i+1}, y\_{i+1}\\)
is the next point, and where \\(x\_{n+1} = x_1, y\_{n+1} = y_1\\).

For multipart polygons with holes, all parts are treated as separate
contributions to the overall centroid, which provides the same result as
if all parts with holes are separately computed, and then merged
together using the parallel axis theorem.

The code and documentation are adapted from the PySAL Python package
(Ray and Anselin, 2007). See Hally (1987) and Li et al. (2013) for
additional details.

## References

Hally, D. 1987. "The calculations of the moments of polygons." Canadian
National Defense Research and Development Technical Memorandum 87/209.
<https://apps.dtic.mil/sti/tr/pdf/ADA183444.pdf>

Li, W., Goodchild, M.F., and Church, R.L. 2013. "An Efficient Measure of
Compactness for Two-Dimensional Shapes and Its Application in
Regionalization Problems." International Journal of Geographical
Information Science 27 (6): 1227–50.
[doi:10.1080/13658816.2012.752093](https://doi.org/10.1080/13658816.2012.752093)
.

Rey, Sergio J., and Luc Anselin. 2007. "PySAL: A Python Library of
Spatial Analytical Methods." Review of Regional Studies 37 (1): 5–27.
[doi:10.52324/001c.8285](https://doi.org/10.52324/001c.8285) .

## See also

[`nmi()`](https://walker-data.com/spopt-r/reference/nmi.md), which
computes the normalized moment of inertia.

## Examples

``` r
library(sf)
poly <- st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)))
second_areal_moment(poly)
#> [1] 0.1666667
```
