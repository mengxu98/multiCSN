# rough_hierarchy

returns rough roots in the network, rough roots selected as those
connected to most number of nodes

## Usage

``` r
rough_hierarchy(network_table, abs_weight = TRUE, directed = TRUE)
```

## Arguments

- network_table:

  The weight data table of network.

- abs_weight:

  Whether to use absolute weight values.

- directed:

  Whether the network is directed.

## Value

list

## Examples

``` r
data("example_matrix", package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix)
#> ℹ [2026-09-11 12:07:25] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 12:07:25] Checking parameters...
#> ✔ [2026-09-11 12:07:25] Inferring network done
#> ℹ [2026-09-11 12:07:25] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
rough_hierarchy(network_table)
#> $roots
#> [1] "g6" "g2" "g4" "g5" "g3" "g1"
#> 
#> $num_paths
#> [1] 6
#> 
```
