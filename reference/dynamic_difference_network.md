# Compute a dynamic difference network

Compute a dynamic difference network

## Usage

``` r
dynamic_difference_network(
  edgeDF,
  epochs,
  condition,
  type,
  diff_thresh = 3,
  condition_thresh = 6
)
```

## Arguments

- edgeDF:

  the result of running edge_uniqueness

- epochs:

  list of epoch gene assignments

- condition:

  condition of interest, should be one of the treatment (aka network
  name) or column names in edgeDF

- type:

  "on" or "off" specifies either finding edges that are active in the
  condition network but off others or inactive in condition network but
  active in others

- diff_thresh:

  edge difference threshold to determine if edge is uniquely on or off

- condition_thresh:

  edge weight theshold in condition network to keep or filter out in
  difference network

## Value

list of dataframes representing the dynamic difference network
