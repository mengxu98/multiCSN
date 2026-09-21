# Classify edges between two network states

Merges two endpoint tables on their key columns and labels every edge as
gained, lost, sign-stable or sign-flipped. This is the shared reduction
behind the state-transition, donor-deletion and cell-budget summaries.

## Usage

``` r
classify_edge_transitions(
  previous,
  next_state,
  keys,
  direction_column = "direction",
  labels = c(gained = "gained", lost = "lost", stable = "shared_sign_stable", flip =
    "shared_sign_flip")
)
```

## Arguments

- previous, next_state:

  Endpoint tables with the key columns and a signed direction column.

- keys:

  Key columns identifying an edge.

- direction_column:

  Column holding the signed direction in both tables.

- labels:

  Named character vector with the class names for `gained`, `lost`,
  `stable` and `flip`.

## Value

A `data.table` with the key columns, `previous_direction`,
`next_direction` and `topology_class`; edges present in only one state
have `NA` on the other side.
