# Returns cuts to define states

Returns cuts to define states via clustering

## Usage

``` r
find_cuts_by_clustering(
  matrix,
  dynamic_object,
  num_states,
  limit_to = NULL,
  method = "kmeans",
  p_value = 0.05
)
```

## Arguments

- matrix:

  genes-by-cells expression matrix

- dynamic_object:

  result of running findDynGenes or define_states

- num_states:

  the number of states

- limit_to:

  vector of genes on which to base state cuts, for example, limiting to
  TFs

- method:

  what clustering method to use, either 'kmeans' or 'hierarchical'

- p_value:

  pval threshold if gene is dynamically expressed

## Value

vector of pseudotimes at which to cut data into states
