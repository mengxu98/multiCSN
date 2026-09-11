# Useful plotting function to plot heatmap of module expression across time with pre-split matrix

Useful plotting function to plot heatmap of module expression across
time with pre-split matrix

## Usage

``` r
heatmap_by_treatment_group(
  expList,
  meta_data,
  pseudotime_column = "pseudotime",
  toScale = T,
  limits = c(0, 5),
  smooth = TRUE,
  order_by = "WAG",
  thresh_on = 0.02,
  fontSize = 8,
  anno_colors = NA
)
```

## Arguments

- expList:

  list of expression matrices

- meta_data:

  sample table

- pseudotime_column:

  column in sample table containing pseudotime annotation

- toScale:

  whether or not to scale expression

- limits:

  limits on expression

- smooth:

  whether or not to smooth expression across pseudotime for cleaner
  plotting

- order_by:

  name in expList that is used to order rows in the heatmap

- thresh_on:

  threshold expression is considered on, used in ordering the rows

- fontSize:

  heatmap font size

- anno_colors:

  annotation colors

## Value

pheatmap
