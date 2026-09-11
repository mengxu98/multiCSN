# Function to return shortest path from 1 regulator to 1 target in a static network

Function to return shortest path from 1 regulator to 1 target in a
static network

## Usage

``` r
static_shortest_path(
  network_table,
  regulator,
  target,
  weight_column = "weight",
  compare_to_average = FALSE
)
```

## Arguments

- network_table:

  The weight data table of network.

- regulator:

  The starting gene.

- target:

  The end gene.

- weight_column:

  column name in network_table with edge weights that will be converted
  to distances

- compare_to_average:

  if TRUE will compute normalized against average path length

## Value

shortest path, distance, normalized distance, and action

## Examples

``` r
data("example_matrix", package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix)
#> ℹ [2026-09-11 13:23:15] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:23:15] Checking parameters...
#> ✔ [2026-09-11 13:23:15] Inferring network done
#> ℹ [2026-09-11 13:23:15] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
static_shortest_path(
  network_table,
  regulator = "g1",
  target = "g2"
)
#> $path
#> [1] "g1" "g6" "g2"
#> 
#> $distance
#> [1] 2.869565
#> 
#> $action
#> [1] 1
#> 
```
