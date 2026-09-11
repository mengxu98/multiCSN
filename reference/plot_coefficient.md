# Plot coefficients

Plot coefficients

## Usage

``` r
plot_coefficient(
  data,
  style = c("continuous", "binary"),
  positive_color = "#3d67a2",
  negative_color = "#c82926",
  neutral_color = "#cccccc",
  bar_width = 0.7,
  text_size = 3,
  show_values = TRUE,
  ...
)
```

## Arguments

- data:

  Input data with `regulator` and `weight` columns.

- style:

  Plotting style: `"binary"` or `"continuous"`. Both use
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html)
  signed value bars; `"binary"` uses the two endpoint colors and
  `"continuous"` includes a neutral midpoint color.

- positive_color:

  Color for positive weights. Default is `"#3d67a2"`.

- negative_color:

  Color for negative weights. Default is `"#c82926"`.

- neutral_color:

  Color for weights near zero (used in `"continuous"`). Default is
  `"#cccccc"`.

- bar_width:

  Width of the bars. Default is `0.7`.

- text_size:

  Size of the text for weight values. Default is `3`.

- show_values:

  Whether to show weight values on bars. Default is `TRUE`.

- ...:

  Additional arguments passed to
  [thisplot::StatPlot](https://mengxu98.github.io/thisplot/reference/StatPlot.html).

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_matrix, package = "inferCSN")
network_table <- inferCSN::inferCSN(example_matrix, targets = "g1")
plot_coefficient(network_table)
plot_coefficient(network_table, style = "binary")
} # }
```
