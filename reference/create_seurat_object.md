# Create a Seurat object from simulated data

Create a Seurat object from simulated data

## Usage

``` r
create_seurat_object(
  sim_data,
  min_cells = 3,
  min_features = 0,
  variable_features = 2000
)
```

## Arguments

- sim_data:

  Output from
  [`simulate_csn_data`](https://mengxu98.github.io/multiCSN/reference/simulate_csn_data.md),
  or a list with `expression_matrix`, `peak_matrix`, `tfs`, `targets`,
  `cell_metadata`, and `peak_annotations`.

- min_cells:

  Include genes/peaks detected in at least this many cells

- min_features:

  Include cells where at least this many genes are detected

- variable_features:

  Number of variable features to select

## Value

A Seurat object with RNA and ATAC assays

## Examples

``` r
if (FALSE) { # \dontrun{
sim_data <- simulate_csn_data()
seurat_obj <- create_seurat_object(sim_data)
} # }
```
