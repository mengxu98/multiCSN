# Assemble factorized TF-region-target support

Assemble factorized TF-region-target support

## Usage

``` r
factorized_chain_support(tf_region, region_gene, tf_gene = NULL)
```

## Arguments

- tf_region:

  TF-region edges with regulator, region, standardized_beta and
  deletion_delta_bic columns.

- region_gene:

  Region-gene edges with region, target, standardized_beta and
  deletion_delta_bic columns.

- tf_gene:

  Optional TF-gene backbone with regulator, target and the same evidence
  columns.

## Value

A chain table ranked by minimum endpoint deletion evidence, with
backbone concordance flags.
