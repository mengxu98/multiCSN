# Find peaks or regions near gene body or TSS

Find peaks or regions near gene body or TSS

## Usage

``` r
find_peaks_near_genes(
  peaks,
  genes,
  sep = c("-", "-"),
  method = c("Signac", "GREAT"),
  upstream = 1e+05,
  downstream = 0,
  extend = 1e+06,
  only_tss = FALSE,
  verbose = TRUE
)
```

## Arguments

- peaks:

  A `GRanges` object with peak regions.

- genes:

  A `GRanges` object with gene coordinates.

- sep:

  Vector of separators to use for genomic string.

- method:

  Character specifying the method to link peak overlapping motif regions
  to nearby genes. One of 'Signac' or 'GREAT'.

- upstream:

  Integer defining the distance upstream of the gene/TSS to consider.

- downstream:

  Integer defining the distance downstream of the gene/TSS to consider.

- extend:

  Integer defining the distance from the upstream and downstream of the
  basal regulatory region. Used only when method is 'GREAT'.

- only_tss:

  Logical value. Measure distance from the TSS (`TRUE`) or from the
  entire gene body (`FALSE`).

- verbose:

  Logical value. Display messages

## Value

A sparse binary Matrix with gene/peak matches.
