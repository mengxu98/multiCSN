# Export CSN

Export network data from a CSN object with optional filtering

## Usage

``` r
export_csn(object, ...)

export_csn(object, ...)

# S4 method for class 'Seurat'
export_csn(
  object,
  active_network = NULL,
  celltypes = NULL,
  weight_cutoff = NULL,
  ...
)

# S4 method for class 'Network'
export_csn(object, weight_cutoff = NULL, ...)
```

## Arguments

- object:

  A CSNObject object

- ...:

  Additional arguments

- active_network:

  Character string specifying which network to export

- celltypes:

  Character vector of cell types to export

- weight_cutoff:

  Numeric threshold for filtering edges by absolute weight

## Value

A list of network data frames by cell type

## Examples

``` r
if (FALSE) { # \dontrun{
# Supply a multiome Seurat object with a completed CSN.

networks <- export_csn(object)
head(networks[[1]])
} # }
```
