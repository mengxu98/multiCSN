# Plot the dynamic differential network but colored by communities and optionally faded by igraph::betweenness

Plot the dynamic differential network but colored by communities and
optionally faded by igraph::betweenness

## Usage

``` r
plot_detail_network(
  network,
  regulators,
  top_edges = NULL,
  only_TFs = TRUE,
  network_order = NULL,
  communities = NULL,
  compute_betweenness = TRUE
)
```

## Arguments

- network:

  the dynamic network

- regulators:

  regulators

- top_edges:

  top_edges

- only_TFs:

  whether or not to only plot regulators and exclude non-regulators

- network_order:

  the network_order in which to plot epochs, or which epochs to plot

- communities:

  community assignments or the result of running find_commumities. The
  names in this object should match the names of the epoch networks in
  network. If NULL, it will be automatically run.

- compute_betweenness:

  whether or not to fade nodes by igraph::betweenness

## Value

plot
