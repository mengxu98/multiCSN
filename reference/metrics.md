# Metrics

Metrics

Get metrics

## Usage

``` r
metrics(object, ...)

metrics(object, ...)

# S4 method for class 'Seurat'
metrics(object, network = DefaultNetwork(object), celltypes = NULL, ...)

# S4 method for class 'Network'
metrics(object, celltypes = NULL, ...)
```

## Arguments

- object:

  The input data, a csn object.

- ...:

  Arguments for other methods

- network:

  network

- celltypes:

  cell types to analyze, NULL for all cell types
