# Get gene ranks from Network/CSNObject

Get gene ranks from Network/CSNObject

## Usage

``` r
GeneRanks(object, ...)

# S4 method for class 'Seurat'
GeneRanks(object, network = DefaultNetwork(object), celltypes = NULL, ...)

# S4 method for class 'Network'
GeneRanks(object, ...)
```

## Arguments

- object:

  The input data.

- ...:

  Arguments for other methods

- network:

  Name of the network.

- celltypes:

  Cell types to get ranks for; NULL for all.
