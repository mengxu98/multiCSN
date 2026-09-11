# Select top prioritized TFs

Select top prioritized TFs

## Usage

``` r
select_key_tfs(
  object,
  network = DefaultNetwork(object),
  pseudotime_column = NULL,
  assay = NULL,
  peak_assay = NULL,
  method = c("page_rank", "degree_distribution"),
  require_all_states = TRUE,
  min_state_fraction = 1,
  min_targets = 10,
  top_n = 6
)
```

## Arguments

- object:

  A `Seurat` object.

- network:

  Dynamic or static network name. Defaults to `DefaultNetwork(object)`.

- pseudotime_column:

  Optional pseudotime column.

- assay:

  RNA assay used for TF expression summaries.

- peak_assay:

  Chromatin assay used for motif-accessibility summaries.

- method:

  Regulator ranking method.

- require_all_states:

  Require TFs to be present in all inferable states.

- min_state_fraction:

  Minimum fraction of states where TF must be present when
  `require_all_states = FALSE`.

- min_targets:

  Minimum number of targets in at least one state.

- top_n:

  Number of TFs returned.

## Value

Character vector of TFs.
