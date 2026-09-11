# Assigns genes to epochs

Assigns genes to epochs

## Usage

``` r
assign_epochs_new1(
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

  individual path result of running define_epochs

- method:

  method of assigning epoch genes, either "active_expression" (looks for
  active expression in epoch) or "DE" (looks for differentially
  expressed genes per epoch)

- p_value:

  p value threshold if gene is dynamically expressed

- pThresh_DE:

  p value if gene is differentially expressed. Ignored if method is
  active_expression.

- active_thresh:

  value between 0 and 1. Percent threshold to define activity

- toScale:

  whether or not to scale the data

- forceGenes:

  whether or not to rescue orphan dyanmic genes, forcing assignment into
  epoch with max expression.

## Value

epochs a list detailing genes active in each epoch
