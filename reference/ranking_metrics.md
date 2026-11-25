# Tie-aware ranking metrics

Ranking metrics for a score vector against binary labels. Scores tied at
a threshold are treated as one block, so every metric is invariant to
the order of tied predictions.

`ranking_metrics()` returns the threshold average precision
(`average_precision`), the trapezoidal area under the precision-recall
curve (`auprc`), the AUROC with average ranks for ties and expected
precision, recall and F1 for each `top_k`, where a top-k boundary that
falls inside a tie block takes the block's expected positives.

## Usage

``` r
ranking_metrics(scores, labels, top_k = c(25L, 50L, 100L, 200L, 500L, 1000L))

tie_aware_average_precision(labels, scores)

tie_aware_auroc(labels, scores)

expected_precision_at_k(labels, scores, k)
```

## Arguments

- scores:

  Numeric vector of prediction scores; larger is more confident.

- labels:

  Logical or numeric 0/1 vector of truth labels, with the same length as
  `scores`.

- top_k:

  Non-empty numeric vector of positive whole-number edge budgets no
  larger than `.Machine$integer.max`.

- k:

  Edge budget for `expected_precision_at_k()`.

## Value

`ranking_metrics()` returns a list with `auroc`, `average_precision`,
`auprc`, a `topk` data frame and the descending score `order`. If either
class is absent, the metrics and cutoff scores are `NA` and the
requested budgets are retained. The other functions return one number
and drop pairs with missing labels or non-finite scores. Average
precision requires a positive label, AUROC requires both classes, and
expected precision at `k` requires a usable pair and a positive budget;
otherwise they return `NA`.

## Details

`ranking_metrics()` requires at least one finite score and rejects
non-finite scores and missing labels. With both classes present, budgets
above the number of scores are reduced to that number.

## Examples

``` r
scores <- c(0.9, 0.8, 0.8, 0.3, 0.1)
labels <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
ranking_metrics(scores, labels, top_k = 2)$topk
#>   top_k cutoff_score precision recall   f1
#> 1     2          0.8      0.75   0.75 0.75
tie_aware_auroc(labels = labels, scores = scores)
#> [1] 0.9166667
tie_aware_average_precision(labels = labels, scores = scores)
#> [1] 0.8333333
expected_precision_at_k(labels = labels, scores = scores, k = 2)
#> [1] 0.75
```
