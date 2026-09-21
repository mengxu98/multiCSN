# Read layered network endpoints from a result directory

Read layered network endpoints from a result directory

## Usage

``` r
read_network_endpoints(
  root,
  endpoints = NULL,
  direction = FALSE,
  weight_column = "weight",
  validate = TRUE
)
```

## Arguments

- root:

  Directory holding the endpoint files of one fitted network.

- endpoints:

  Endpoints to read; defaults to all five.

- direction:

  Replace the signed weight column with a `direction` column holding its
  sign (finite, non-zero), which is what the transition and stability
  helpers consume.

- weight_column:

  Name of the signed weight column.

- validate:

  Check that key columns are present, complete and unique.

## Value

A named list of `data.table`s, one per endpoint.
