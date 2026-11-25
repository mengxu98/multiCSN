# Validate a layered network result directory

Checks that every endpoint and accounting file is present and non-empty,
that tables can be read, that endpoint keys are complete and unique, and
that finite nonnegative integer summary counts match the endpoint row
counts. Accounting tables retain their endpoint, target and
selected-count columns. Header-only tables are accepted for complete
empty networks. With `error = FALSE`, malformed files and invalid result
metadata return `FALSE` without signaling a parsing or validation error.

## Usage

``` r
validate_network_endpoints(root, format = NULL, error = TRUE)
```

## Arguments

- root:

  Directory holding the endpoint files of one fitted network.

- format:

  File format: tsv, csv, tsv.gz or csv.gz. Detected from the files in
  `root` when `NULL`.

- error:

  Signal an error on the first problem; otherwise return `FALSE`.

## Value

`TRUE` (invisibly) for a complete result.
