# Extract transition-specific regulators and targets

Extract transition-specific regulators and targets

## Usage

``` r
extract_genes_transition(object, ...)

# S4 method for class 'list'
extract_genes_transition(
  object,
  regulators = NULL,
  targets = NULL,
  transition_networks = NULL,
  top_regulators_num = 5,
  top_targets_num = 5,
  common_regulators = TRUE,
  common_targets = TRUE
)

# S4 method for class 'Seurat'
extract_genes_transition(
  object,
  network = NULL,
  celltypes = NULL,
  weight_cutoff = NULL,
  regulators = NULL,
  targets = NULL,
  transition_networks = NULL,
  top_regulators_num = 5,
  top_targets_num = 5,
  common_regulators = TRUE,
  common_targets = TRUE
)

# S4 method for class 'CSNObject'
extract_genes_transition(object, ...)
```

## Arguments

- object:

  Named list of network tables or CSNObject.

- ...:

  Passed to methods.

- regulators:

  Optional regulator subset; kept for backward compatibility and
  currently not used in the ranking.

- targets:

  Optional target subset; kept for backward compatibility and currently
  not used in the ranking.

- transition_networks:

  Optional list of state triplets defining transitions; `NULL` derives
  triplets from consecutive states.

- top_regulators_num:

  Number of top regulators kept per state.

- top_targets_num:

  Number of top targets kept per regulator.

- common_regulators:

  Logical; when `TRUE`, restrict regulators to those shared across the
  states of a transition.

- common_targets:

  Logical; when `TRUE`, split targets into transition-specific sets and
  a shared `intersect` set.

- network:

  Name of dynamic network in `object@networks`.

- celltypes:

  States to include; `NULL` for all (sorted by name for transition
  triplets).

- weight_cutoff:

  Passed to `export_csn`.
