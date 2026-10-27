# Write layered network result tables

Writes the five network endpoint tables, three per-target accounting
tables and the endpoint summary. TF-region accounting uses `target` as
the response column name. No Seurat object or checkpoint files are
written.

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

  The result returned by
  [`fit_layered_network()`](https://mengxu98.github.io/multiCSN/reference/fit_layered_network.md).

- root:

  Output directory.

- format:

  File format: tsv, csv, tsv.gz or csv.gz.

- overwrite:

  Replace existing output files. Defaults to `FALSE`.

## Value

An invisible named character vector of output paths.

## Details

Use the same `format` with
[`read_network_endpoints()`](https://mengxu98.github.io/multiCSN/reference/read_network_endpoints.md)
to read the five network tables. Accounting and summary tables can be
read as ordinary delimited tables. Empty networks retain their column
headers.
