# Get stored extract_genes_transition results from CSNObject

Get stored extract_genes_transition results from CSNObject

## Usage

``` r
ExtractGenesTransition(object, ...)

# S4 method for class 'Seurat'
ExtractGenesTransition(object, network = DefaultNetwork(object), ...)
```

## Arguments

- object:

  The input data.

- ...:

  Arguments for other methods

- network:

  Name of the stored network whose transition results are returned;
  defaults to `DefaultNetwork(object)`.
