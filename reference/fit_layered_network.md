# Fit a layered multi-omic gene regulatory network

Fits the layered multiCSN model on single-cell multiome data: a
transcription-factor-to-gene layer, a transcription-factor-to-region
layer and a region-to-gene layer, each solved as independent per-target
Greedy-l0 problems, followed by the factorized chain projection that
links the three layers. By default only the selected cells are used and
the expression layers stored in the object are taken as-is.

## Usage

``` r
fit_layered_network(
  object,
  cell_group = NULL,
  cells = NULL,
  regulators = NULL,
  targets = NULL,
  min_detected = 0L,
  min_peak_cells = 3L,
  upstream = 1e+05,
  downstream = 0L,
  tf_region_scope = c("motif_peaks", "domain_peaks"),
  exclude_response_alias = FALSE,
  renormalize = FALSE,
  response_chunk = 1024L,
  max_support_size = NULL,
  min_improvement = 1e-10,
  sort_regulators = FALSE,
  cores = 1L,
  checkpoint_dir = NULL,
  store = FALSE,
  network_name = "multiCSN",
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with chromatin accessibility and gene expression
  assays plus the multiCSN state (see
  [`initiate_object`](https://mengxu98.github.io/multiCSN/reference/initiate_object.md)).

- cell_group:

  Catalogue key (cell type, state or timepoint) used to select cells,
  peaks and target genes.

- cells:

  Optional character vector of exact cells to fit; when supplied the fit
  uses exactly these cells and validates them against `cell_group`.

- regulators, targets:

  Optional character vectors restricting the transcription factors and
  target genes.

- min_detected:

  Minimum number of cells in which a target gene or a transcription
  factor must be detected; `0` disables the filter.

- min_peak_cells:

  Minimum number of cells in which a candidate regulatory region must be
  detected.

- upstream, downstream:

  Window in base pairs used to assign candidate regions to genes.

- tf_region_scope:

  Candidate regions of the TF-region layer: `"motif_peaks"` fits every
  motif-carrying accessible region, while `"domain_peaks"` restricts the
  layer to regions inside gene domains.

- exclude_response_alias:

  Drop candidates whose standardized cross-product with the response is
  an affine alias (\|r\| \>= 1 - 1e-10).

- renormalize:

  Re-normalize RNA counts and re-run TF-IDF on the selected cells before
  fitting.

- response_chunk:

  Number of targets fitted per chunk.

- max_support_size, min_improvement:

  Greedy-l0 stopping rules forwarded to
  [`inferCSN::fit_greedy_l0_batch()`](https://mengxu98.github.io/inferCSN/reference/fit_greedy_l0_batch.html).

- sort_regulators:

  Sort the transcription factors before fitting.

- cores:

  Number of forked workers used per chunk (ignored on Windows).

- checkpoint_dir:

  Optional directory storing per-chunk checkpoints so an interrupted fit
  can resume; `NULL` keeps everything in memory.

- store:

  Store the three layers as `Network` entries in the returned object.

- network_name:

  Name of the stored network entries.

- verbose:

  Print progress messages.

## Value

A list with `tf_gene`, `tf_region`, `region_gene`, `triplets`,
`mediated_tf_gene`, `accounting`, `summary`, `settings`, `cells` and,
when `store = TRUE`, `object`.
