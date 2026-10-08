# Define states

Define states

## Usage

``` r
define_states(
  dynamic_object,
  matrix,
  method = "pseudotime",
  num_states = 2,
  pseudotime_cuts = NULL,
  group_assignments = NULL,
  p_value = 0.05,
  winSize = 2
)
```

## Arguments

- dynamic_object:

  A trajectory result with a named vector of gene P values in `genes`
  and a cell table in `cells`, or a named list of such results. The cell
  table contains `pseudotime` and cell identifiers; group-based
  partitioning also requires `group`.

- matrix:

  genes-by-cells expression matrix, or a list of expression matrices per
  path. If list, names should match names of dynamic_object.

- method:

  method to define states. Either "pseudotime", "cell_order", "group",
  "con_similarity", "kmeans", "hierarchical"

- num_states:

  number of states to define. Ignored when explicit pseudotime cuts or
  group assignments define the states.

- pseudotime_cuts:

  vector of pseudotime cutoffs. If NULL, cuts are set to
  max(pseudotime)/num_states.

- group_assignments:

  a list of vectors where names(assignment) are state names, and vectors
  contain groups belonging to corresponding state

- p_value:

  p_value

- winSize:

  winSize

## Value

updated list of dynamic_object with state column included in
dynamic_object\$cells
