# Plot gene ranks and network properties

Plot gene ranks and network properties

## Usage

``` r
plot_gene_rank(object, ...)

# S4 method for class 'Network'
plot_gene_rank(
  object,
  method = c("page_rank", "degree_distribution"),
  weight_cutoff = 0.1,
  compare_random = TRUE,
  top_n = 30,
  ...
)

# S4 method for class 'Seurat'
plot_gene_rank(
  object,
  network = NULL,
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  weight_cutoff = 0.1,
  compare_random = TRUE,
  combine = TRUE,
  top_n = 30,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  ...
)

# S4 method for class 'CSNObject'
plot_gene_rank(object, ...)

# S4 method for class 'data.frame'
plot_gene_rank(
  object,
  method = c("page_rank", "degree_distribution"),
  weight_cutoff = 0.1,
  compare_random = TRUE,
  top_n = 30
)
```

## Arguments

- object:

  Network object

- ...:

  Other params

- method:

  Character, ranking method: "page_rank" or "degree_distribution"

- weight_cutoff:

  Numeric, threshold for edge weight filtering

- compare_random:

  Logical, whether to compare with randomized network

- top_n:

  Integer, number of top genes to show in centrality plot.

- network:

  Name of the network to use; defaults to the active network.

- celltypes:

  Character vector specifying which states/cell types to include; NULL
  for all.

- combine:

  Logical. When `TRUE`, combine per-state plots with
  [`patchwork::wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html);
  otherwise return a list of plots.

- nrow, ncol:

  Integer or `NULL`. Rows/columns for
  [`patchwork::wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
  when `combine=TRUE`.

- byrow:

  Logical. Fill plots by row in `wrap_plots`; default `TRUE`.

## Value

Combined ggplot object
