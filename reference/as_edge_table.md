# Standardize an edge table

Accepts the common layouts of predicted networks and ground truths
(`regulator`/`target`, `Gene1`/`Gene2`, or the first two columns) and
returns `regulator`, `target` and, for predictions, a numeric `weight`.
The weight column is the first of `weight`, `score`, `cor_value`,
`importance`, `corr` and `value`, then the third column, otherwise 1.

## Usage

``` r
as_edge_table(table, weight = TRUE)
```

## Arguments

- table:

  A data frame of edges.

- weight:

  Return a `weight` column.

## Value

A data frame.

## Examples

``` r
as_edge_table(data.frame(Gene1 = "A", Gene2 = "B", EdgeWeight = 0.4))
#>   regulator target weight
#> 1         A      B    0.4
```
