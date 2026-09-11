# Computes Jaccard similarity between top regulators in two sets of networks across a range of top X regulators

Computes Jaccard similarity between top regulators in two sets of
networks across a range of top X regulators

## Usage

``` r
JI_across_topregs(
  netlist1,
  netlist2,
  n_regs = 3:15,
  func = "mean",
  method = "pagerank",
  weight_column = "zscore",
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

  a vector indicating which values of top regulators to scan across

- func:

  func

- method:

  method to find top regulators. Currently only supports "pagerank"

- weight_column:

  column name in grnDFs containing edge weights

- compare_within_netlist1:

  whether or not to do pairwise comparisons between networks in netlist1

- compare_within_netlist2:

  whether or not to do pairwise comparisons between networks in netlist2

## Value

dataframe of Jaccard similarities of top regulators
