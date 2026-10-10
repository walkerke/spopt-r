# Compute Huff Model probabilities

Computes probability surface based on distance decay and attractiveness.
Formula: \\(P\_{ij} = (A_j \times D\_{ij}^\beta) / \Sigma_k(A_k \times
D\_{ik}^\beta)\\)

## Usage

``` r
rust_huff(cost_matrix, attractiveness, distance_exponent, sales_potential)
```

## Arguments

- cost_matrix:

  Cost/distance matrix (demand x stores)

- attractiveness:

  Attractiveness values for each store (pre-computed with exponents)

- distance_exponent:

  Distance decay exponent (typically negative, e.g., -1.5)

- sales_potential:

  Optional sales potential for each demand point

## Value

List with probabilities, market shares, expected sales
