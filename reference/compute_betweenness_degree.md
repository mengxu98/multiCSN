# Function to compute betweenness and degree

Function to compute betweenness and degree

## Usage

``` r
compute_betweenness_degree(
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

list of each network, ranked by \`betweenness\*degree\`
