# Plot dynamic networks

Plot dynamic networks

## Usage

``` r
plot_static_networks(
  network_table,
  regulators = NULL,
  targets = NULL,
  legend_position = "right"
)
```

## Arguments

- network_table:

  The weight data table of network.

- regulators:

  A character vector of regulators to include.

- targets:

  A character vector of targets to include.

- legend_position:

  The position of legend.

## Value

A ggplot2 object

## Examples

``` r
data(example_matrix, package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix)
#> ℹ [2026-09-11 13:23:11] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:23:11] Checking parameters...
#> ✔ [2026-09-11 13:23:11] Inferring network done
#> ℹ [2026-09-11 13:23:11] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
plot_static_networks(
  network_table,
  regulators = "g1"
)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the multiCSN package.
#>   Please report the issue at <https://github.com/mengxu98/multiCSN/issues>.

plot_static_networks(
  network_table,
  targets = "g1"
)

plot_static_networks(
  network_table,
  regulators = "g2",
  targets = "g3"
)
```
