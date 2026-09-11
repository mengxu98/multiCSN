# Function to compute page rank of TF+target networks

Function to compute page rank of TF+target networks

## Usage

``` r
compute_pagerank(
  networks_list,
  regulators = NULL,
  targets = NULL,
  directed = FALSE
)
```

## Arguments

- networks_list:

  result of dynamic network reconstruction

- regulators:

  list of regulators

- targets:

  list of targets

- directed:

  if network is directed or not

## Value

list of each network, ranked by \`page_rank\`
