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
  target_chunk = 256L,
  max_support_size = NULL,
  min_improvement = 1e-10,
  sort_regulators = FALSE,
  cores = 1L,
  checkpoint_dir = NULL,
  store = FALSE,
  network_name = "multiCSN",
  verbose = TRUE,
  preset = NULL,
  output_dir = NULL,
  format = c("tsv", "csv", "tsv.gz", "csv.gz")
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

- preset:

  Default settings for `min_detected`, `tf_region_scope`,
  `exclude_response_alias`, `renormalize` and `sort_regulators`:
  `"state"` or `"timepoint"`. `NULL` keeps the individual defaults.

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

  Number of responses fitted per task in the shared-design layers.

- target_chunk:

  Number of region-gene targets fitted per task; larger values amortise
  the fork cost of the parallel backend.

- max_support_size, min_improvement:

  Greedy-l0 stopping rules forwarded to
  [`inferCSN::fit_greedy_l0_batch()`](https://mengxu98.github.io/inferCSN/reference/fit_greedy_l0_batch.html).

- sort_regulators:

  Sort the transcription factors before fitting.

- cores:

  Number of workers used per chunk (fork on Unix, PSOCK on Windows).

- checkpoint_dir:

  Optional directory storing per-chunk checkpoints so an interrupted fit
  can resume. Existing checkpoints are reused only when the selected
  inputs and fitting settings match. A SHA-256 identity is computed when
  checkpointing is enabled; `NULL` keeps everything in memory.

- store:

  Store the three layers as `Network` entries in the returned object.
  The default `FALSE` returns only the fitted result and does not
  include a copy of the Seurat object.

- network_name:

  Name of the stored network entries.

- output_dir:

  Optional directory receiving the endpoint tables (see
  [`write_network_endpoints`](https://mengxu98.github.io/multiCSN/reference/write_network_endpoints.md)).
  The tables are written to a temporary sibling directory and published
  only after the fit completes. Use `checkpoint_dir` to retain resumable
  checkpoints. An existing `output_dir` is never overwritten.

- format:

  File format of `output_dir`.

- verbose:

  Print progress messages.

## Value

A list with `tf_gene`, `tf_region`, `region_gene`, `triplets`,
`mediated_tf_gene`, `accounting`, `summary`, `settings`, `cells`,
`targets`, `regulators` and, when `store = TRUE`, `object`. With
`output_dir` the list is returned invisibly and also holds the written
`paths`.

## Details

A typical workflow is
[`prepare_candidates()`](https://mengxu98.github.io/multiCSN/reference/prepare_candidates.md)
followed by
`fit_layered_network(object, cell_group = "state", output_dir = "out")`.
Two presets cover the common designs: `"state"` fits cell states that
share one normalization, and `"timepoint"` re-normalizes each group,
requires every target and TF to be detected in 20 cells and drops affine
response aliases. Arguments supplied explicitly override the preset.

Each assay must have one layer named `data`; split layers require
[`SeuratObject::JoinLayers()`](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
before fitting. Per-group normalization also requires a single `counts`
layer in each assay. With `store = FALSE`, fitting leaves the supplied
object unchanged.

Each edge layer retains its own `standardized_beta`,
`deletion_delta_bic` and signed ordinal `weight`. The coefficient sign
gives that layer's fitted direction. Ordinal weights rank deletion
evidence within an endpoint and retain coefficient signs. They are not
coefficients or probabilities.

`triplets` links selected TF-region and region-gene edges. Its
`chain_delta_bic` is the minimum of the two endpoint deletion scores;
`chain_direction` is the product of their signs. `backbone_selected`
records whether the direct TF-gene edge exists, and
`backbone_sign_concordant` records agreement with the direct sign.
`mediated_tf_gene` keeps the strongest sign-consistent path per TF-gene
pair. These projections retain their own scores and support flags.

`settings` records the effective settings after applying the preset and
explicit arguments. `cells`, `targets` and `regulators` record the
entities actually used. The five endpoint tables, accounting tables and
summary can be saved independently with
[`write_network_endpoints()`](https://mengxu98.github.io/multiCSN/reference/write_network_endpoints.md)
and read with
[`read_network_endpoints()`](https://mengxu98.github.io/multiCSN/reference/read_network_endpoints.md).

Valid inputs can yield empty endpoint tables. If no TF has variable
expression, the TF layers have no selected edges and their accounting
status is `no_variable_candidate`; the region-gene layer is still fitted
independently.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- prepare_candidates(object, group.by = "cell_type", genome = genome)
fit <- fit_layered_network(
  object,
  cell_group = "Erythroblast",
  preset = "state",
  output_dir = "networks/erythroblast"
)
networks <- read_network_endpoints("networks/erythroblast")
} # }
```
