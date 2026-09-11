# Plot a target-gene locus using Signac when available, with robust fallbacks

Plot a target-gene locus using Signac when available, with robust
fallbacks

## Usage

``` r
plot_target_locus_signac(
  object,
  target_gene,
  key_target_tbl,
  tfs = NULL,
  pseudotime_column,
  assay = "RNA",
  peak_assay = NULL,
  extend_upstream = 5000,
  extend_downstream = 5000,
  auto_extend = TRUE,
  include_unassigned = FALSE,
  top_links = 8,
  show_tf_labels = TRUE,
  normalize_tracks = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- target_gene:

  Target gene symbol.

- key_target_tbl:

  Output of
  [`prioritize_features`](https://mengxu98.github.io/multiCSN/reference/prioritize_features.md)
  with `type = "targets"` and `strategy = "integrated"`.

- tfs:

  Optional TF subset used to build links.

- pseudotime_column:

  Pseudotime column used for state assignment.

- assay:

  RNA assay used for expression tracks when `CoveragePlot()` is
  available.

- peak_assay:

  Peak assay name.

- extend_upstream, extend_downstream:

  Genomic extension around the gene.

- auto_extend:

  Whether to expand the plotting window to include all inferred linked
  peaks for the target gene.

- include_unassigned:

  Whether to keep unassigned cells in grouped tracks.

- top_links:

  Maximum number of links shown in fallback link plots.

- show_tf_labels:

  Whether to display TF labels above linked peaks in the fallback link
  panel.

- normalize_tracks:

  Whether to scale each grouped pseudo-coverage track to a maximum of 1
  in fallback mode.

## Value

A patchwork/ggplot object.
