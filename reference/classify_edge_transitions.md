# Classify edges between two network states

Merges two endpoint tables on their key columns and labels every edge as
gained, lost, sign-stable or sign-flipped. This is the shared reduction
behind the state-transition, donor-deletion and cell-budget summaries.

## Usage

``` r
classify_edge_transitions(
  previous,
  next_state,
  keys = c("regulator", "target"),
  direction_column = "direction",
  labels = c(gained = "gained", lost = "lost", stable = "shared_sign_stable", flip =
    "shared_sign_flip"),
  from = NULL,
  to = NULL
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

- from, to:

  Optional state names added as `previous_state`, `next_state` and
  `transition` (`"from -> to"`) columns.

## Value

A `data.table` with the key columns, `previous_direction`,
`next_direction` and `topology_class`; edges present in only one state
have `NA` on the other side.

## Examples

``` r
early <- data.frame(regulator = c("A", "B"), target = "X", direction = c(1, -1))
late <- data.frame(regulator = c("B", "C"), target = "X", direction = c(1, 1))
classify_edge_transitions(
  early, late, c("regulator", "target"),
  from = "early", to = "late"
)
#>    regulator target previous_direction next_direction   topology_class
#>       <char> <char>              <num>          <num>           <char>
#> 1:         A      X                  1             NA             lost
#> 2:         B      X                 -1              1 shared_sign_flip
#> 3:         C      X                 NA              1           gained
#>    previous_state next_state    transition
#>            <char>     <char>        <char>
#> 1:          early       late early -> late
#> 2:          early       late early -> late
#> 3:          early       late early -> late
```
