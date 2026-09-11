# Aggregate Seurat assay over groups

Aggregate Seurat assay over groups

## Usage

``` r
aggregate_assay(
  object,
  group_name,
  fun = "mean",
  assay = "RNA",
  layer = "data"
)
```

## Arguments

- object:

  The csn object.

- group_name:

  A character vector indicating the metadata column to aggregate over.

- fun:

  The summary function to be applied to each group.

- assay:

  The assay to summarize.

- layer:

  The layer to summarize.

## Value

A Seurat object.
