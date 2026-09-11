# Computes frobenius distance in a pairwise manner between two sets of networks

Computes frobenius distance in a pairwise manner between two sets of
networks

## Usage

``` r
compute_frobenius_distance(
  netlist1,
  netlist2,
  weight_column = "weight",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
)
```

## Arguments

- netlist1:

  list of grnDFs

- netlist2:

  list of grnDFs

- weight_column:

  column name containing edge weights

- compare_within_netlist1:

  whether or not to do pairwise comparisons between networks in netlist1

- compare_within_netlist2:

  whether or not to do pairwise comparisons between networks in netlist2

## Value

dataframe of frobenius distances
