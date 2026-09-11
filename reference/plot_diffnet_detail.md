# Plot the dynamic differential network but colored by communities and optionally faded by betweenness

Plot the dynamic differential network but colored by communities and
optionally faded by betweenness

## Usage

``` r
plot_diffnet_detail(
  grn,
  tfs,
  only_TFs = TRUE,
  order = NULL,
  compute_betweenness = TRUE
)
```

## Arguments

- grn:

  the dynamic network

- tfs:

  TFs

- only_TFs:

  whether or not to only plot TFs and exclude non-regulators

- order:

  the order in which to plot epochs, or which epochs to plot

- compute_betweenness:

  whether or not to fade nodes by betweenness

## Value

plot
