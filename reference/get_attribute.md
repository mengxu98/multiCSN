# Get any attribute from a CSNObject object

Get any attribute from a CSNObject object

Get any attribute from a CSNObject

## Usage

``` r
get_attribute(object, ...)

get_attribute(object, ...)

# S4 method for class 'Seurat'
get_attribute(
  object,
  celltypes = NULL,
  active_network = NULL,
  attribute = c("genes", "tfs", "peaks", "regulators", "targets", "celltypes", "cells",
    "modules", "coefficients"),
  ...
)

# S4 method for class 'CSNObject'
get_attribute(object, ...)
```

## Arguments

- object:

  The input data, a csn object.

- ...:

  Arguments for other methods

- celltypes:

  A character vector specifying the celltypes to get attributes for.

- active_network:

  A character string specifying the active network.

- attribute:

  Attribute to get: genes, tfs, peaks, regulators, targets, celltypes,
  cells, modules, coefficients.

## Value

A character vector of regulatory genes

A character vector or list of regulatory genes
