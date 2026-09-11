# Splits data into epochs manually

Splits data into epochs given group assignment

## Usage

``` r
split_epochs_by_group(dynamic_object, assignment)
```

## Arguments

- dynamic_object:

  result of running findDynGenes or compileDynGenes

- assignment:

  a list of vectors where names(assignment) are epoch names, and vectors
  contain groups belonging to corresponding epoch

## Value

updated dynamic_object with epoch column included in
dynamic_object\$cells
