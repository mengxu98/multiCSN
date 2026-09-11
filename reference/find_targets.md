# Finds binding targets given list of dataframes containing binding info for effectors

Finds binding targets given list of dataframes containing binding info
for effectors

## Usage

``` r
find_targets(
  aList,
  column = "max_score",
  by_rank = FALSE,
  n_targets = 2000,
  threshold = NULL
)
```

## Arguments

- aList:

  the result of running score_targets

- column:

  column name used in either ranking or thresholding the possible
  targets

- by_rank:

  if TRUE will return top n_targets; if FALSE will use thresholds

- n_targets:

  the number of targets to return for each effector if by_rank=TRUE

- threshold:

  threshold to use if by_rank=FALSE

## Value

list of binding targets for list of effectors
