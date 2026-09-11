# Edge uniqueness across multiple GRNs

Edge uniqueness across multiple GRNs

## Usage

``` r
edge_uniqueness(object, ...)

# S4 method for class 'list'
edge_uniqueness(object, tfs, weight_column)

# S4 method for class 'Seurat'
edge_uniqueness(
  object,
  network = NULL,
  celltypes = NULL,
  tfs = NULL,
  weight_column = "weight",
  weight_cutoff = NULL
)

# S4 method for class 'CSNObject'
edge_uniqueness(object, ...)
```

## Arguments

- object:

  Named list of GRN data.frames (regulator, target, weight) or
  CSNObject.

- ...:

  Passed to methods.

- tfs:

  If `NULL`, uses `object@metadata$tfs`.

- weight_column:

  Column name for edge weights.

- network:

  Name of dynamic network in `object@networks`.

- celltypes:

  States to include; `NULL` for all.

- weight_cutoff:

  Passed to `export_csn`.
