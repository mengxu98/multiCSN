# Get stored shortest-path results from Network/CSNObject

Get stored shortest-path results from Network/CSNObject

## Usage

``` r
ShortestPaths(object, ...)

# S4 method for class 'Seurat'
ShortestPaths(object, network = DefaultNetwork(object), celltypes = NULL, ...)

# S4 method for class 'Network'
ShortestPaths(object, ...)
```

## Arguments

- object:

  The input data.

- ...:

  Arguments for other methods

- network:

  Name of the network.

- celltypes:

  Cell types/states to get paths for; NULL for all.
