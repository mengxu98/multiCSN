# Layered multi-omic networks

multiCSN fits three regression layers on paired single-cell RNA and
chromatin data, transcription factor (TF) to gene, TF to region and
region to gene, and links them into factorized TF-region-gene chains.
This vignette runs the workflow on simulated data.

## Input

The input is a Seurat object with a normalized RNA assay, a TF-IDF
normalized chromatin assay with a gene annotation, and a metadata column
with the cell groups to fit.

``` r

library(multiCSN)
set.seed(1)
simulated <- simulate_csn_data(n_cells = 300, n_genes = 80, n_tfs = 30, n_peaks = 3000)
object <- suppressWarnings(suppressMessages(create_seurat_object(simulated, variable_features = 50)))
object <- Signac::RunTFIDF(object, assay = "ATAC", verbose = FALSE)
object
#> An object of class Seurat 
#> 3080 features across 300 samples within 2 assays 
#> Active assay: RNA (80 features, 50 variable features)
#>  3 layers present: counts, data, scale.data
#>  1 other assay present: ATAC
```

## Candidate catalogue

[`prepare_candidates()`](https://mengxu98.github.io/multiCSN/reference/prepare_candidates.md)
keeps every RNA gene as a candidate target, scans the peaks for motifs
with the bundled `motifs` and `motif2tf` resources and records the TFs
with a motif.

``` r

genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
object <- prepare_candidates(object, group.by = "celltype", genome = genome, verbose = FALSE)
get_attribute(object, attribute = "celltypes")
#> [1] "celltype_3" "celltype_1" "celltype_2"
length(get_attribute(object, attribute = "tfs"))
#> [1] 80
```

## Fitting

Fit one network per cell group. `preset = "state"` suits cell states
sharing one normalization; `preset = "timepoint"` re-normalizes the
group and requires targets and TFs to be detected in 20 cells. Any
argument you pass overrides the preset.

``` r

output <- file.path(tempdir(), "celltype_1")
fit <- fit_layered_network(
  object,
  cell_group = "celltype_1",
  preset = "state",
  upstream = 5e7,
  output_dir = output,
  verbose = FALSE
)
fit$summary
#>           endpoint selected
#> 1          TF-gene      100
#> 2        TF-region      310
#> 3      region-gene       50
#> 4 factorized-chain       16
#> 5 mediated-TF-gene       16
```

With `output_dir`, the tables are written only after the fit completes,
and a re-run never overwrites an existing result. Each table has signed
ordinal weights in `weight`.

## Reading results

``` r

validate_network_endpoints(output)
networks <- read_network_endpoints(output)
names(networks)
#> [1] "TF-gene"          "TF-region"        "region-gene"      "factorized-chain"
#> [5] "mediated-TF-gene"
head(networks$`TF-gene`)
#>    regulator target standardized_beta deletion_delta_bic weight
#>       <char> <char>             <num>              <num>  <num>
#> 1:     SOX15  DMBX1         0.2432997           1.041878  0.215
#> 2:     FOXA3   DRGX         0.3358523           8.708487  0.795
#> 3:      GLI3   DRGX         0.1966337           0.234426  0.065
#> 4:     LM187   DRGX         0.2264841           1.719762  0.335
#> 5:    TFAP2A   DRGX         0.2308855           1.898556  0.375
#> 6:    ZNF101   DRGX         0.2909058           5.354368  0.715
```

## Comparing states

``` r

other <- fit_layered_network(
  object, cell_group = "celltype_2", preset = "state", upstream = 5e7, verbose = FALSE
)
first <- read_network_endpoints(output, "TF-gene", direction = TRUE, simplify = TRUE)
second <- data.frame(
  regulator = other$tf_gene$regulator,
  target = other$tf_gene$target,
  direction = sign(other$tf_gene$weight)
)
transitions <- classify_edge_transitions(first, second, from = "celltype_1", to = "celltype_2")
table(transitions$topology_class)
#> 
#>             gained               lost shared_sign_stable 
#>                 73                 99                  1
```
