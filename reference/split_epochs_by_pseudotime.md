# Splits data into epochs

Splits data into epochs by assigning cells to epochs

## Usage

``` r
split_epochs_by_pseudotime(dynamic_object, cuts, epoch_names = NULL)
```

## Arguments

- dynamic_object:

  result of running findDynGenes or compileDynGenes

- cuts:

  vector of pseudotime cutoffs

- epoch_names:

  names of resulting epochs, must have length of length(cuts)+1

## Value

updated dynamic_object with epoch column included in
dynamic_object\$cells
