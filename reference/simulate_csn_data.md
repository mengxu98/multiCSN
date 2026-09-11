# Simulate data for testing inferCSN

Simulate data for testing inferCSN

## Usage

``` r
simulate_csn_data(
  n_cells = 300,
  n_genes = 100,
  n_tfs = 50,
  n_peaks = 1000,
  n_celltypes = 3,
  sparsity = 0.9,
  atac_sparsity = 0.9
)
```

## Arguments

- n_cells:

  Number of cells to simulate

- n_genes:

  Number of genes to simulate

- n_tfs:

  Number of transcription factors

- n_peaks:

  Number of ATAC peaks to simulate

- n_celltypes:

  Number of cell types to simulate

- sparsity:

  Sparsity level of the expression matrix (0-1)

- atac_sparsity:

  Sparsity level of the ATAC matrix (0-1)

## Value

A list containing simulated data
