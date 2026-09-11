# Select top prioritized targets for chosen TFs

Select top prioritized targets for chosen TFs

## Usage

``` r
select_key_targets(
  object,
  tfs = NULL,
  tf_scores = NULL,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  top_n = 10
)
```

## Arguments

- object:

  A `Seurat` object.

- tfs:

  Optional TF vector. If `NULL`, derived from `tf_scores`.

- tf_scores:

  Optional output of
  `prioritize_features(type = "regulators", strategy = "integrated")`.

- network:

  Dynamic or static network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- assay:

  RNA assay used for target-expression support.

- peak_assay:

  Chromatin assay used for peak support.

- top_n:

  Number of targets kept per TF.

## Value

A data.frame of top targets per TF.
