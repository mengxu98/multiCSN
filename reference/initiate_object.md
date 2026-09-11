# Initiate object

Initiate object

Initiate the `CSNObject` object

## Usage

``` r
initiate_object(object, ...)

initiate_object(object, ...)

# S4 method for class 'matrix'
initiate_object(object, regulators = character(0), targets = character(0), ...)

# S4 method for class 'Seurat'
initiate_object(
  object,
  celltype = NULL,
  regulators = NULL,
  targets = NULL,
  group.by = NULL,
  filter_mode = c("celltype_specific", "variable", "celltype", "unfiltered"),
  filter_by = c("celltype", "aggregate"),
  rna_assay = "RNA",
  rna_min_pct = 0.1,
  rna_logfc_threshold = 0.25,
  rna_test_method = "wilcox",
  n_variable_genes = 2000,
  peak_assay = NULL,
  peak_min_pct = 0.05,
  peak_logfc_threshold = 0.1,
  peak_test_method = "LR",
  n_variable_peaks = "q5",
  regions = NULL,
  regions_extend = 0,
  exclude_exons = TRUE,
  only_pos = TRUE,
  verbose = TRUE,
  p_value = 0.05,
  ...
)

# S4 method for class 'CSNObject'
initiate_object(object, ...)
```

## Arguments

- object:

  A Seurat object containing gene expression and/or chromatin
  accessibility data.

- ...:

  Additional arguments passed to marker detection functions.

- regulators:

  Regulators to consider for CSN inference. Can be a character vector
  applied to all cell types, or a named list (names must match cell
  types) of character vectors per cell type.

- targets:

  Targets to consider for CSN inference. Same format rules as
  `regulators`.

- celltype:

  Selected celltype(s) to analyze. If NULL, all cells will be used.

- group.by:

  Column name in Seurat meta.data for cell type information.

- filter_mode:

  Filter mode for identifying cell-type specific features.

- filter_by:

  When filter_mode is "variable", specify "aggregate" or "celltype".

- rna_assay:

  Name of the gene expression assay.

- rna_min_pct, rna_logfc_threshold, rna_test_method, n_variable_genes:

  RNA params.

- peak_assay:

  Name of the chromatin accessibility assay. If NULL, only RNA.

- peak_min_pct, peak_logfc_threshold, peak_test_method,
  n_variable_peaks:

  Peak params.

- regions:

  Candidate regions (GRanges or data frame). If NULL, all peaks.

- regions_extend:

  Base pairs to extend regions.

- exclude_exons:

  Whether to consider exons for binding site inference.

- only_pos:

  Only return positive markers.

- verbose:

  Print progress messages.

- p_value:

  Significance threshold for filtering features.

## Value

CSNObject object.

## Examples

``` r
data(pbmcmultiome_sub, package = "scop")

object <- Seurat::NormalizeData(pbmcmultiome_sub, assay = "RNA", verbose = FALSE)
object <- Signac::RunTFIDF(object, assay = "peaks", verbose = FALSE)

object <- initiate_object(
  object,
  group.by = "CellType",
  rna_assay = "RNA",
  peak_assay = "peaks",
  verbose = FALSE
)
object
#> An object of class Seurat 
#> 24405 features across 500 samples within 2 assays 
#> Active assay: peaks (12000 features, 0 variable features)
#>  2 layers present: counts, data
#>  1 other assay present: RNA
```
