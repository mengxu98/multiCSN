# Endpoint tables of a layered multi-omic network

The five endpoint tables written by
[`fit_layered_network`](https://mengxu98.github.io/multiCSN/reference/fit_layered_network.md)
(and by the manuscript pipeline) share one contract: a file name, a key
column set and a signed weight. This function returns that contract so
writers, readers and validators cannot drift apart.

## Usage

``` r
network_endpoints()
```

## Value

A data frame with columns `endpoint` and `file` plus a list column
`keys` holding the key columns of each endpoint.
