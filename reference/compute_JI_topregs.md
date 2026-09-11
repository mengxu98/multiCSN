# Computes Jaccard similarity between top regulators in two sets of networks

Computes Jaccard similarity between top regulators in two sets of
networks

## Usage

``` r
compute_JI_topregs(
  netlist1,
  netlist2,
  n_regs = 15,
  method = "pagerank",
  compare_within_netlist1 = TRUE,
  compare_within_netlist2 = TRUE
)
```

## Arguments

- netlist1:

  list of grnDFs

- netlist2:

  list of grnDFs

- n_regs:

  the number of top regulators to compare from each network

- method:

  method to find top regulators. Currently only supports "pagerank"

- compare_within_netlist1:

  whether or not to do pairwise comparisons between networks in netlist1

- compare_within_netlist2:

  whether or not to do pairwise comparisons between networks in netlist2

## Value

dataframe of Jaccard similarities of top regulators
