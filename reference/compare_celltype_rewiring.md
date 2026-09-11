# Compare rewiring between static cell type networks

Compare rewiring between static cell type networks

## Usage

``` r
compare_celltype_rewiring(
  object,
  network = DefaultNetwork(object),
  celltypes = NULL,
  weight_column = "weight"
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Static network name. Defaults to `DefaultNetwork(object)`.

- celltypes:

  Optional subset of cell types.

- weight_column:

  Weight column used when exporting networks.

## Value

A list with pair summaries and per-TF rewiring tables.
