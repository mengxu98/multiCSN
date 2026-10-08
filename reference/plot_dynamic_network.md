# quick plot of dynamic networks

quick plot of dynamic networks

## Usage

``` r
plot_dynamic_network(
  network,
  regulators,
  only_TFs = TRUE,
  network_order = NULL,
  weight_threshold = NULL
)
```

## Arguments

- network:

  the result of running split_network_by_states

- regulators:

  regulators

- only_TFs:

  plot only regulator network

- network_order:

  which states or transitions to plot

- weight_threshold:

  weight_threshold

## Value

plot
