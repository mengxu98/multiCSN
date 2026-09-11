# find_cuts_by_similarity

Returns cuts to define epochs via sliding window comparison

## Usage

``` r
find_cuts_by_similarity(
  matrix,
  dynamic_object,
  winSize = 2,
  limit_to = NULL,
  p_value = 0.05
)
```

## Arguments

- matrix:

  genes-by-cells expression matrix

- dynamic_object:

  result of running findDynGenes or define_epochs

- winSize:

  number of cells to each side to compare each cell to

- limit_to:

  vector of genes on which to base epoch cuts, for example, limiting to
  TFs

- p_value:

  pval threshold if gene is dynamically expressed

## Value

vector of pseudotimes at which to cut data into epochs
