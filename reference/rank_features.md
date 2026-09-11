# Rank network features across ordered states

Rank network features across ordered states

## Usage

``` r
rank_features(
  object,
  type = c("regulators", "targets"),
  network = DefaultNetwork(object),
  celltypes = NULL,
  method = c("page_rank", "degree_distribution"),
  top_n = 20,
  weight_column = "weight"
)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Feature type, `"regulators"` or `"targets"`.

- network:

  Dynamic network name. Defaults to `DefaultNetwork(object)`.

- celltypes:

  Optional subset of states.

- method:

  Gene-rank method, `"page_rank"` or `"degree_distribution"`.

- top_n:

  Number of rows to keep per state. Use `NULL` for all.

- weight_column:

  Weight column used when `type = "targets"`.

## Value

A data.frame combining ranked features for the selected states.
