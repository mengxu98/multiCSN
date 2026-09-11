# Computes betweenness and degree of each regulator for each network in a list of networks

Computes betweenness and degree of each regulator for each network in a
list of networks

## Usage

``` r
biglist_compute_betweenness_degree(
  netlist,
  weight_column = "zscore",
  normalized = TRUE
)
```

## Arguments

- netlist:

  netlist

- weight_column:

  column name in grnDFs containing edge weights

- normalized:

  whether or not to normalize degree and betweenness

## Value

dataframe listing network, betweenness of each regulator, degree of each
regulator
