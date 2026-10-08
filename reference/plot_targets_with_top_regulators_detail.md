# quick plot of top regulators given targets in dynamic networks based on reconstruction weight, colored by expression and interaction type

quick plot of top regulators given targets in dynamic networks based on
reconstruction weight, colored by expression and interaction type

## Usage

``` r
plot_targets_with_top_regulators_detail(
  network,
  targets,
  states_list,
  weight_column = "weight",
  gene_ranks = NULL,
  regulators_num = 5,
  network_order = NULL,
  fixed_layout = TRUE,
  layout_alg = "fr",
  declutter = TRUE
)
```

## Arguments

- network:

  the result of running split_network_by_states

- targets:

  targets

- states_list:

  result of running assign_genes_to_states

- weight_column:

  column name containing reconstruction weights to use

- gene_ranks:

  gene_ranks

- regulators_num:

  number of top regulators to plot

- network_order:

  which states or transitions to plot

- fixed_layout:

  whether or not to fix node positions across state networks

- layout_alg:

  layout algorithm if fixed_layout. Defaults to FR. Ignored if
  fixed_layout==FALSE.

- declutter:

  if TRUE, will only label nodes with active interactions in given
  network

## Value

plot
