# Useful plotting function to plot heatmap with pre-split matrix

Useful plotting function to plot heatmap with pre-split matrix

## Usage

``` r
plot_heatmap_by_treatment(
  expList,
  meta_data,
  pseudotime_column = "latent_time",
  toScale = T,
  limits = c(0, 5),
  smooth = TRUE,
  anno_colors = NA,
  fontSize = 8
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

- anno_colors:

  annotation colors

- fontSize:

  heatmap font size

## Value

pheatmap
