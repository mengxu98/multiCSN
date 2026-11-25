# Write the endpoint tables of a layered network

Write the endpoint tables of a layered network

## Usage

``` r
write_network_endpoints(
  fit,
  root,
  format = c("tsv", "csv", "tsv.gz", "csv.gz"),
  overwrite = FALSE
)
```

## Arguments

- fit:

  A result of
  [`fit_layered_network`](https://mengxu98.github.io/multiCSN/reference/fit_layered_network.md).

- root:

  Output directory.

- format:

  File format: tsv, csv, tsv.gz or csv.gz.

- overwrite:

  Replace existing files.

## Value

The written paths, invisibly.
