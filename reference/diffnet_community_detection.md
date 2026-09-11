# Perform community detection on a dynamic network

Perform community detection on a dynamic network

## Usage

``` r
diffnet_community_detection(
  diffnet,
  method = "louvain",
  use_weights = FALSE,
  weight_column = NULL
)
```

## Arguments

- diffnet:

  diffnet

- method:

  community detection method, currently only "louvain"

- use_weights:

  whether or not to use edge weights (for weighted graphs)

- weight_column:

  if using weights, name of the column containing edge weights

## Value

community assignments of nodes in the dynamic network
