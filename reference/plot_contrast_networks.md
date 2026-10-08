# Plot contrast networks

Plot contrast networks

## Usage

``` r
plot_contrast_networks(
  network_table,
  degree_value = 0,
  weight_value = 0,
  cols = NULL,
  regulator_color = "#8C4985",
  target_color = "#B0C4DE",
  label.size = 3.5,
  base_family = "",
  legend_position = "bottom"
)
```

## Arguments

- network_table:

  The weight data table of network.

- degree_value:

  Degree value to filter nodes. Default is `0`.

- weight_value:

  Weight value to filter edges. Default is `0`.

- cols:

  Named vector of colors for edge interactions. Default uses
  `c("Activation" = "#1F78B4", "Repression" = "#E31A1C")`.

- regulator_color:

  Color of regulator nodes. Default is `"#8C4985"`.

- target_color:

  Color of target nodes. Default is `"#B0C4DE"`.

- label.size:

  Node label font size. Default is `3.5`.

- base_family:

  Font family for graph theme. Default is `""`.

- legend_position:

  The position of legend. Default is `"bottom"`.

## Value

A ggplot2 object.

## Examples

``` r
data(example_matrix, package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix)
#> ℹ [2026-10-08 08:48:43] Inferring network for <matrix/array>...
#> ✔ [2026-10-08 08:48:43] Inferring network done
#> ℹ [2026-10-08 08:48:43] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6
plot_contrast_networks(network_table[1:50, ])
```
