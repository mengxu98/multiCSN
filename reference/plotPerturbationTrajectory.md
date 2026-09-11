# Plot perturbation trajectory

Plot perturbation trajectory

Plot perturbation trajectory

## Usage

``` r
plotPerturbationTrajectory(
  object,
  genes = NULL,
  show_heatmap = FALSE,
  n_top_genes = 50,
  heatmap_style = "combined",
  ...
)

plotPerturbationTrajectory(
  object,
  genes = NULL,
  show_heatmap = FALSE,
  n_top_genes = 50,
  heatmap_style = "combined",
  ...
)

# S4 method for class 'Network'
plotPerturbationTrajectory(
  object,
  genes = NULL,
  show_heatmap = FALSE,
  n_top_genes = 50,
  heatmap_style = "combined",
  ...
)

# S4 method for class 'Seurat'
plotPerturbationTrajectory(
  object,
  genes = NULL,
  show_heatmap = FALSE,
  n_top_genes = 50,
  heatmap_style = "combined",
  ...
)
```

## Arguments

- object:

  Network or CSNObject

- genes:

  Genes to plot

- show_heatmap:

  Logical, whether to show expression heatmap

- n_top_genes:

  Number of top variable genes to show in heatmap

- heatmap_style:

  Style of heatmap visualization ("combined" or "separate")

- ...:

  Additional arguments

## Value

ggplot object or list of plots
