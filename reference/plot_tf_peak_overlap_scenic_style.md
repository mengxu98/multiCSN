# Plot TF-supported peak overlap across ordered states

Plot TF-supported peak overlap across ordered states

## Usage

``` r
plot_tf_peak_overlap_scenic_style(
  object,
  key_target_tbl,
  tfs,
  pseudotime_column,
  peak_assay = NULL,
  top_targets_per_tf = 6
)
```

## Arguments

- object:

  A `Seurat` object.

- key_target_tbl:

  Output of
  [`prioritize_features`](https://mengxu98.github.io/multiCSN/reference/prioritize_features.md)
  with `type = "targets"` and `strategy = "integrated"`.

- tfs:

  TFs to visualize.

- pseudotime_column:

  Pseudotime column used for state assignment.

- peak_assay:

  Peak assay name.

- top_targets_per_tf:

  Number of target genes retained per TF.

## Value

A `ggplot` object.
