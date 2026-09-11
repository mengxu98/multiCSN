# Function that orders genes based on peak expression

Function that orders genes based on peak expression

## Usage

``` r
order_genes(matrix, meta_data, pseudotime_column, smooth = TRUE)
```

## Arguments

- matrix:

  expression matrix

- meta_data:

  sample table

- pseudotime_column:

  column in sample table containing pseudotime annotation

- smooth:

  whether or not to smooth expression across pseudotime

## Value

ordered genes
