# Normalized moment of inertia (NMI)

Computes the normalized moment of inertia (NMI), a compactness measure
for polygon geometries. The NMI ranges between 0 and 1, where 1 is the
most compact shape (a circle) and 0 is an infinitely extending shape
(Feng et al. 2022).

## Usage

``` r
nmi(x)
```

## Arguments

- x:

  An sf object, sfc geometry column, or sfg geometry

## Value

Numeric vector of normalized moments of inertia.

## Details

The NMI is defined as follows, where \\(A\\) is the area of a geometry,
and \\(I\\) is the second moment of inertia (i.e., the second areal
moment): \$\$\frac{A^2}{2 \pi I}\$\$ See Li et al. (2013, 2014) for
additional details.

## References

Feng, X., Rey, S., and Wei, R. (2022). "The max-p-compact-regions
problem." Transactions in GIS, 26, 717–734.
[doi:10.1111/tgis.12874](https://doi.org/10.1111/tgis.12874) .

Li, W., Goodchild, M.F., and Church, R.L. 2013. "An Efficient Measure of
Compactness for Two-Dimensional Shapes and Its Application in
Regionalization Problems." International Journal of Geographical
Information Science 27 (6): 1227–50.
[doi:10.1080/13658816.2012.752093](https://doi.org/10.1080/13658816.2012.752093)
.

Li, W., Church, R.L. and Goodchild, M.F. 2014. "The p-Compact-regions
Problem." Geogr Anal, 46: 250-273.
[doi:10.1111/gean.12038](https://doi.org/10.1111/gean.12038) .

## See also

[`second_areal_moment()`](https://walker-data.com/spopt-r/reference/second_areal_moment.md)
