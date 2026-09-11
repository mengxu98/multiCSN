# find_communities

Find communities in a static or dynamic network

## Usage

``` r
find_communities(network_table, use_weights = FALSE, directed = FALSE)
```

## Arguments

- network_table:

  network_table

- use_weights:

  whether or not to use edge weights (for weighted graphs)

- directed:

  Default set to \`FALSE\`. If directed = TRUE, remember to flip target
  and regulator in diffnet dfs

## Value

community assignments of nodes in the dynamic network
