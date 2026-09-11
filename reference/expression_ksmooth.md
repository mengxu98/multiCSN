# expression_ksmooth

smooths expression across cells in path.

## Usage

``` r
expression_ksmooth(
  matrix,
  meta_data,
  pseudotime_column = "pseudotime",
  bandwith = 0.25
)
```

## Arguments

- matrix:

  expression matrix.

- meta_data:

  meta_data includes "cell_name", "pseudotime", "group".

- pseudotime_column:

  pseudotime_column.

- bandwith:

  bandwith for kernel smoother.

## Value

smoothed matrix
