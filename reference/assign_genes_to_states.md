# Assigns genes to states

Assigns genes to states

## Usage

``` r
assign_genes_to_states(
  matrix,
  dynamic_object,
  method = "active_expression",
  p_value = 0.05,
  pThresh_DE = 0.05,
  active_thresh = 0.33,
  toScale = FALSE,
  forceGenes = TRUE
)
```

## Arguments

- matrix:

  genes-by-cells expression matrix

- dynamic_object:

  individual path result of running define_states

- method:

  method of assigning state genes, either "active_expression" (looks for
  active expression in state) or "DE" (looks for differentially
  expressed genes per state)

- p_value:

  pval threshold if gene is dynamically expressed

- pThresh_DE:

  pval if gene is differentially expressed. Ignored if method is
  active_expression.

- active_thresh:

  value between 0 and 1. Percent threshold to define activity

- toScale:

  whether or not to scale the data

- forceGenes:

  whether or not to rescue orphan dyanmic genes, forcing assignment into
  state with max expression.

## Value

states a list detailing genes active in each state
