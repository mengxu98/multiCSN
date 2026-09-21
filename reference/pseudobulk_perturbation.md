# Pseudobulk a perturbation screen into a paired design

Sums single-cell counts of one or more modalities per replicate x
condition sample, which is the input of the paired effect model used to
validate an inferred network against a perturbation screen. The control
label keeps its own sample in every replicate, so the resulting design
is paired and complete by construction.

## Usage

``` r
pseudobulk_perturbation(
  cells,
  sample_info,
  targets = NULL,
  control_label = "NT",
  roles = c(control = "network_training_control", perturbed = "heldout_perturbation"),
  require_complete = TRUE,
  verbose = TRUE
)
```

## Arguments

- cells:

  A matrix (features x cells) or a named list of matrices, one per
  modality (for example RNA and ATAC). All modalities must share the
  cells.

- sample_info:

  Data frame with `cell_id`, `replicate`, `timepoint` and `target`
  columns.

- targets:

  Perturbation targets to aggregate; defaults to every target of
  `sample_info` except the control label, in order of appearance.

- control_label:

  Label identifying the unperturbed control cells.

- roles:

  Named character vector mapping `control` and `perturbed` to the role
  labels stored in the sample metadata.

- require_complete:

  Fail when a replicate lacks a control or target sample.

- verbose:

  Print progress messages.

## Value

A list with `sample_metadata` and `pseudobulk` (one matrix per modality,
columns ordered exactly like `sample_metadata$sample_id`).
