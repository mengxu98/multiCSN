# Function to score targets of effectors

Adds columns for mean binding score, maximum binding score, and percent
frequency target is a hit based on threshold

## Usage

``` r
score_targets(aList, threshold = 50)
```

## Arguments

- aList:

  list of dataframes containing binding data for each effector

- threshold:

  binding score threshold to compute hit frequency

## Value

updated list of dataframes
