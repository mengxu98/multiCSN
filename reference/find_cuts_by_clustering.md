# Returns cuts to define epochs

Returns cuts to define epochs via clustering

## Usage

``` r
find_cuts_by_clustering(
  matrix,
  dynamic_object,
  num_epochs,
  limit_to = NULL,
  method = "kmeans",
  p_value = 0.05
)
```

## Arguments

- matrix:

  genes-by-cells expression matrix

- dynamic_object:

  result of running findDynGenes or define_epochs

- num_epochs:

  the number of epochs

- limit_to:

  vector of genes on which to base epoch cuts, for example, limiting to
  TFs

- method:

  what clustering method to use, either 'kmeans' or 'hierarchical'

- p_value:

  pval threshold if gene is dynamically expressed

## Value

vector of pseudotimes at which to cut data into epochs
