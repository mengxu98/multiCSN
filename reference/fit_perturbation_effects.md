# Fit paired perturbation effects on pseudobulk counts

Fits one quasi-likelihood negative-binomial model per timepoint and
target with a paired design (`~ replicate + condition`) and reports
per-feature log fold changes and FDR values for the
perturbed-minus-control contrast.

## Usage

``` r
fit_perturbation_effects(
  pseudobulk,
  sample_metadata,
  targets = NULL,
  timepoints = NULL,
  control_label = "NT",
  effect_definition = "heldout_perturbed_minus_NT_control",
  min_replicates = 3L,
  output_dir = NULL,
  output_subdirs = NULL,
  verbose = TRUE
)
```

## Arguments

- pseudobulk:

  Named list of pseudobulk matrices (see
  [`pseudobulk_perturbation`](https://mengxu98.github.io/multiCSN/reference/pseudobulk_perturbation.md)).

- sample_metadata:

  Sample metadata with `sample_id`, `replicate`, `timepoint`, `target`,
  `cells`.

- targets, timepoints:

  Subsets to fit; defaults to every target and timepoint present.

- control_label:

  Label of the control samples.

- effect_definition:

  Label stored in the effect tables.

- min_replicates:

  Minimum number of paired replicates per contrast.

- output_dir:

  Optional directory receiving one `<timepoint>__<target>.tsv.gz` table
  per contrast; when `NULL` the tables are returned in memory.

- output_subdirs:

  Optional named vector mapping modality names to the subdirectory
  created below `output_dir` (defaults to the modality name).

- verbose:

  Print progress messages.

## Value

A list with `effects` (per modality, one table or path per contrast) and
`summary` (one row per modality x timepoint x target).
