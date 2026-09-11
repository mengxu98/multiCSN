# inferring cell-type specific gene regulatory network

Infers cell-specific gene regulatory networks. Matrix methods come from
inferCSN; `Network` and `Seurat` methods are implemented here.

## Usage

``` r
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)

# S4 method for class 'Network'
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  r_squared_threshold = 0,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  gene_cor_threshold = 0,
  ...
)

# S4 method for class 'Seurat'
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  penalty = "L0",
  cross_validation = FALSE,
  seed = 1,
  n_folds = 5,
  r_squared_threshold = 0,
  celltypes = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  alpha = 0.5,
  family = "gaussian",
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  pseudotime_column = NULL,
  dynamic_features = NULL,
  ...
)

# S4 method for class 'CSNObject'
inferCSN(
  object,
  pseudotime = NULL,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  Expression matrix, `Network`, `Seurat`, or `CSNObject`. Matrix input
  is handled by
  [`inferCSN::inferCSN()`](https://mengxu98.github.io/inferCSN/reference/inferCSN.html).

- pseudotime:

  Optional pseudotime vector or branch matrix for matrix input. For
  `Seurat` objects, a single column name is treated as
  `pseudotime_column`.

- regulators, targets:

  Optional gene subsets.

- max_support_size:

  Optional support-size limit for matrix inference.

- lag_fraction:

  Fractional lag used when `lag_steps` is `NULL`.

- lag_steps:

  Optional integer lag for matrix inference.

- cores:

  Number of workers.

- verbose:

  Whether to report progress.

- ...:

  Additional method arguments.

- penalty, cross_validation, seed, n_folds:

  Stored on `Network` objects for provenance. Fitting uses `method`.

- r_squared_threshold:

  Threshold of \\R^2\\ used when formatting `Network` / `Seurat`
  results.

- method:

  Fitting backend for `Network` / `Seurat` methods: `"greedy_l0"`
  (greedy L0 via
  [`inferCSN::fit_greedy_l0()`](https://mengxu98.github.io/inferCSN/reference/fit_greedy_l0.html)),
  `"glm"`, `"glmnet"`, `"cv.glmnet"`, `"xgb"`, or `"susie"`.

- gene_cor_threshold:

  Threshold for TF-target gene correlation.

- celltypes:

  Character vector of cell types to infer networks for.

- network_name:

  Name of the stored network.

- peak_to_gene_method:

  Method to link peaks to genes. One of `Signac` or `GREAT`.

- upstream:

  Distance upstream of the gene to consider as a potential regulatory
  region.

- downstream:

  Distance downstream of the gene to consider as a potential regulatory
  region.

- extend:

  Distance added to the basal regulatory region when
  `peak_to_gene_method = "GREAT"`.

- only_tss:

  Measure distance from the TSS (`TRUE`) or from the gene body
  (`FALSE`).

- peak_to_gene_domains:

  `GenomicRanges` object with regulatory regions for each gene.

- peak_cor_threshold:

  Threshold for binding peak-target gene correlation.

- aggregate_rna_col:

  Metadata column used to aggregate RNA.

- aggregate_peaks_col:

  Metadata column used to aggregate peaks.

- alpha:

  Elastic-net mixing parameter. See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html).

- family:

  Error distribution and link function. See
  [`family`](https://rdrr.io/r/stats/family.html).

- interaction_term:

  Interaction between TF and binding site: `"+"`, `":"`, or `"*"`.

- adjust_method:

  Method for adjusting p-values.

- scale:

  Whether to z-transform expression and accessibility matrices.

- pseudotime_column:

  Metadata column for pseudotime. If set, runs state-dynamic inference.

- dynamic_features:

  Optional dynamic-gene filter when `pseudotime_column` is set.

## Value

A regulator-target edge table for matrix input, a `Network` object, or a
`Seurat` object.

## Details

The Seurat method with chromatin peaks fits TF-by-region interaction
predictors through `fit_models`; projecting those coefficients does not
fit independent TF-gene, TF-region and region-gene models. For
factorized analyses, fit the three endpoint responses independently and
use `factorized_chain_support` to assemble their selected supports.

## Examples

``` r
data("example_matrix", package = "inferCSN")
network_table_1 <- inferCSN(
  example_matrix
)
#> ℹ [2026-09-11 13:27:18] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:27:18] Checking parameters...
#> ✔ [2026-09-11 13:27:18] Inferring network done
#> ℹ [2026-09-11 13:27:18] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6

network_table_2 <- inferCSN(
  example_matrix,
  cores = 2
)
#> ℹ [2026-09-11 13:27:18] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:27:18] Checking parameters...
#> ✔ [2026-09-11 13:27:18] Inferring network done
#> ℹ [2026-09-11 13:27:18] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1    12          6       6

head(network_table_1)
#>   regulator target     weight
#> 1        g6     g2 -0.9583333
#> 2        g2     g6 -0.8750000
#> 3        g6     g1  0.7916667
#> 4        g4     g5  0.7083333
#> 5        g5     g4  0.6250000
#> 6        g2     g3 -0.5416667

identical(
  network_table_1,
  network_table_2
)
#> [1] TRUE

inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g4")
)
#> ℹ [2026-09-11 13:27:18] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:27:18] Checking parameters...
#> ✔ [2026-09-11 13:27:18] Inferring network done
#> ℹ [2026-09-11 13:27:18] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     4          2       2
#>   regulator target weight
#> 1        g2     g3 -0.875
#> 2        g1     g4 -0.625
#> 3        g2     g4 -0.375
#> 4        g1     g3 -0.125
inferCSN(
  example_matrix,
  regulators = c("g1", "g2"),
  targets = c("g3", "g0")
)
#> ℹ [2026-09-11 13:27:18] Inferring network for <matrix/array>...
#> ◌ [2026-09-11 13:27:18] Checking parameters...
#> Warning: Ignoring 1 requested targets absent from `object`: g0
#> ✔ [2026-09-11 13:27:18] Inferring network done
#> ℹ [2026-09-11 13:27:18] Network information:
#> ℹ                         Edges Regulators Targets
#> ℹ                       1     2          2       1
#>   regulator target weight
#> 1        g2     g3  -0.75
#> 2        g1     g3  -0.25

if (FALSE) { # \dontrun{
data("example_ground_truth", package = "inferCSN")
calculate_metrics(
  network_table_1,
  example_ground_truth,
  return_plot = TRUE
)
} # }
if (FALSE) { # \dontrun{
data("example_matrix", package = "inferCSN")
object <- initiate_object(example_matrix)
object <- inferCSN(object)
} # }
if (FALSE) { # \dontrun{
# multiome input: motif scanning comes first
data(pbmcmultiome_sub, package = "scop")
data(motifs)
data(motif2tf)

object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)
object <- initiate_object(
  object,
  group.by = "CellType",
  rna_assay = "RNA",
  peak_assay = "peaks",
  verbose = FALSE
)

genome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
object <- find_motifs(
  object,
  pfm = motifs,
  motif_tfs = motif2tf,
  genome = genome,
  backend = "motifmatchr",
  verbose = FALSE
)

object <- inferCSN(object, cores = 2, verbose = FALSE)
networks <- export_csn(object)
head(networks[[1]])
} # }
```
