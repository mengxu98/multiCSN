# Plot coefficients for multiple targets

Plot coefficients for multiple targets

## Usage

``` r
plot_coefficients(data, targets = NULL, nrow = NULL, ...)
```

## Arguments

- data:

  Input data.

- targets:

  Targets to plot. Default is \`NULL\`.

- nrow:

  Number of rows for the plot. Default is \`NULL\`.

- ...:

  Other arguments passed to \[plot_coefficient\].

## Value

A patchwork of ggplot objects

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_matrix, package = "inferCSN")
network_table <- inferCSN::inferCSN(
  example_matrix,
  targets = c("g1", "g2", "g3")
)
plot_coefficients(network_table, show_values = FALSE)
plot_coefficients(network_table, targets = "g1")
} # }
```
