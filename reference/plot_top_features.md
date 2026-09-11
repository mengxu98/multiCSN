# quick plot of top regulators in dynamic networks

quick plot of top regulators in dynamic networks

## Usage

``` r
plot_top_features(
  network_list,
  gene_ranks = NULL,
  regulators = NULL,
  targets = NULL,
  regulators_num = NULL,
  targets_num = NULL,
  network_order = NULL,
  method = "weight"
)
```

## Arguments

- network_list:

  the result of running epochGRN

- gene_ranks:

  the result of running compute_pagerank

- regulators:

  regulators

- targets:

  targets

- regulators_num:

  number of top regulators to plot

- targets_num:

  number of top targets to plot for each regulator

- network_order:

  which epochs or transitions to plot

- method:

  Choose method to rank regulators for plot. Defaule set to \`weight\`,
  mean choose regulators rely weights in network infer by
  [inferCSN](https://mengxu98.github.io/inferCSN/reference/inferCSN.html).
  Also can choose \`page_rank\`, then will ranking regulators rely the
  result of [page_rank](https://r.igraph.org/reference/page_rank.html).

## Value

plot
