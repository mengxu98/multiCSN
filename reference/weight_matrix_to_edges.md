# Convert a square weight matrix to an edge table

Convert a square weight matrix to an edge table

## Usage

``` r
weight_matrix_to_edges(weight_matrix, remove_self = TRUE)
```

## Arguments

- weight_matrix:

  Target-by-regulator weight matrix (rows are targets, columns are
  regulators).

- remove_self:

  Drop self-edges.

## Value

A data frame with `regulator`, `target` and `weight`, sorted by
decreasing absolute weight.

## Examples

``` r
m <- matrix(c(0, 0.2, -0.7, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
weight_matrix_to_edges(m)
#>   regulator target weight
#> 1         B      A   -0.7
#> 2         A      B    0.2
```
