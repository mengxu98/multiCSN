# Subsample an expression matrix

Subsample an expression matrix

## Usage

``` r
subsampling(
  matrix,
  method = c("sample", "meta_cells", "pseudobulk"),
  ratio = 1,
  seed = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- matrix:

  Input matrix.

- method:

  Subsampling strategy.

- ratio:

  Fraction retained or aggregated.

- seed:

  Random seed.

- verbose:

  Whether to report progress.

- ...:

  Additional arguments.

## Value

The subsampled matrix.
