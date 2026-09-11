# Fit models for gene expression

Fit statistical models to infer gene regulatory relationships

## Usage

``` r
fit_models(
  object,
  regulators = NULL,
  targets = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0.1,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  verbose = TRUE,
  cores = 1,
  celltype = NULL,
  ...
)

# S4 method for class 'Seurat'
fit_models(
  object,
  regulators = NULL,
  targets = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0.1,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  verbose = TRUE,
  cores = 1,
  celltype = NULL,
  ...
)

# S4 method for class 'CSNObject'
fit_models(
  object,
  regulators = NULL,
  targets = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0.1,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  verbose = TRUE,
  cores = 1,
  celltype = NULL,
  ...
)

# S4 method for class 'matrix'
fit_models(
  object,
  regulators = NULL,
  targets = NULL,
  network_name = paste0(method, "_network"),
  peak_to_gene_method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  peak_to_gene_domains = NULL,
  gene_cor_threshold = 0.1,
  peak_cor_threshold = 0.1,
  aggregate_rna_col = NULL,
  aggregate_peaks_col = NULL,
  method = c("greedy_l0", "glm", "glmnet", "cv.glmnet", "xgb", "susie"),
  interaction_term = ":",
  adjust_method = "fdr",
  scale = FALSE,
  verbose = TRUE,
  cores = 1,
  celltype = NULL,
  ...
)
```

## Arguments

- object:

  The input data for model fitting

- regulators:

  A character vector of regulators

- targets:

  A character vector of targets

- network_name:

  Name for the inferred network

- peak_to_gene_method:

  Method to link peaks to genes ("Signac" or "GREAT")

- upstream:

  Distance upstream of TSS to consider

- downstream:

  Distance downstream of TSS to consider

- extend:

  Extension distance for GREAT method

- only_tss:

  Whether to measure distance from TSS only

- peak_to_gene_domains:

  Custom regulatory domains

- gene_cor_threshold:

  Correlation threshold for TFs

- peak_cor_threshold:

  Correlation threshold for peaks

- aggregate_rna_col:

  Column name for RNA aggregation

- aggregate_peaks_col:

  Column name for peaks aggregation

- method:

  Statistical method for model fitting

- interaction_term:

  Type of interaction variable

- adjust_method:

  Method for p-value adjustment

- scale:

  Whether to scale the data

- verbose:

  Whether to show progress messages

- cores:

  Number of cores for parallel processing

- celltype:

  Current cell type being processed

- ...:

  Additional arguments passed to methods

## Value

A list containing model coefficients and fit statistics
