# BEELINE-style network evaluation

Scores a predicted network against a ground truth on the complete
directed gene-pair universe of the ground-truth genes (self-pairs
excluded). Predictions are ranked by absolute weight and pairs without a
prediction share score 0. Ranking metrics are tie-aware. The binary
summaries (accuracy, precision, recall, F1, Jaccard index and signal
index) use the Youden-optimal score threshold, and early precision keeps
every non-zero prediction tied at the ground-truth edge count, dividing
by the number actually selected.

## Usage

``` r
evaluate_network(network_table, ground_truth)
```

## Arguments

- network_table:

  Predicted edges, see
  [`as_edge_table`](https://mengxu98.github.io/multiCSN/reference/as_edge_table.md).

- ground_truth:

  Ground-truth edges with non-missing, non-empty endpoints.

## Value

A named list with `AUROC`, `AUPRC`, `AP`, `ACC`, `Precision`, `Recall`,
`F1`, `JI`, `SI`, `EarlyPrecision`, `EPR`, `TopK`, `EarlySelectedEdges`,
`EarlyTruePositives`, `PositiveRate` and `GoldEdges`.

## Examples

``` r
truth <- data.frame(regulator = c("A", "A", "B"), target = c("B", "C", "C"))
prediction <- data.frame(
  regulator = c("A", "B", "C", "A"),
  target = c("B", "C", "A", "C"),
  weight = c(0.9, 0.5, 0.4, 0.1)
)
str(evaluate_network(prediction, truth))
#> List of 16
#>  $ AUROC             : num 0.889
#>  $ AUPRC             : num 0.903
#>  $ AP                : num 0.917
#>  $ ACC               : num 0.833
#>  $ Precision         : num 1
#>  $ Recall            : num 0.667
#>  $ F1                : num 0.8
#>  $ JI                : num 0.667
#>  $ SI                : int 2
#>  $ EarlyPrecision    : num 0.667
#>  $ EPR               : num 1.33
#>  $ TopK              : int 3
#>  $ EarlySelectedEdges: int 3
#>  $ EarlyTruePositives: int 2
#>  $ PositiveRate      : num 0.5
#>  $ GoldEdges         : int 3
```
