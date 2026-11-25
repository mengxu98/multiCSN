# Read layered network endpoints from a result directory

Read layered network endpoints from a result directory

## Usage

``` r
read_network_endpoints(
  root,
  endpoints = NULL,
  direction = FALSE,
  weight_column = "weight",
  validate = TRUE,
  format = NULL,
  simplify = FALSE
)
```

## Arguments

- root:

  Directory holding the endpoint files of one fitted network.

- endpoints:

  Endpoints to read; defaults to all five.

- direction:

  Return only the key columns plus a `direction` column holding the sign
  of `weight_column` (finite, non-zero).

- weight_column:

  Name of the signed weight column.

- validate:

  Check that key columns are present, complete and unique.

- format:

  File format: tsv, csv, tsv.gz or csv.gz. Detected from the files in
  `root` when `NULL`.

- simplify:

  Return the table itself when one endpoint is requested.

## Value

A named list of `data.table`s, one per endpoint. When a single endpoint
is requested with `simplify = TRUE`, that table.
