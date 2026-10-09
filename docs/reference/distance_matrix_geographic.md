# Compute distance matrix for geographic CRS using sf::st_distance

Internal function that uses sf::st_distance() for accurate great circle
distances on geographic (longlat) coordinate systems.

## Usage

``` r
distance_matrix_geographic(x, y = NULL, use_centroids = NULL)
```

## Arguments

- x:

  An sf object

- y:

  An sf object or NULL

- use_centroids:

  Logical. If TRUE, use geometry centroids.

## Value

A numeric matrix of distances in meters.
