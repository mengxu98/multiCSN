# Calculate Early Precision Ratio

Calculates the Early Precision Ratio (EPR) on the fixed candidate
universe. EPR compares the precision among the top-ranked predicted
edges to the precision expected from a random predictor over the same
candidate edge pool. All non-zero predictions tied at the reference
edge-count cutoff are retained; the precision denominator is the actual
selected count. Empty or all-zero predictions have EPR zero. Saved score
precision is retained.

## Usage

``` r
calculate_epr(network_table, ground_truth, tf_edges = FALSE)
```

## Arguments

- network_table:

  A data frame of predicted network structure

- ground_truth:

  A data frame of ground truth network

- tf_edges:

  Whether to restrict the candidate universe to regulator-to-gene edges

## Value

A list containing the metric
